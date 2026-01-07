use super::input::Search;
use super::output::SageResults;
use super::telemetry;
use anyhow::Context;
use csv::ByteRecord;
use log::info;
use rayon::prelude::*;
use sage_cloudpath::{CloudPath, FileFormat};
use sage_core::database::{IndexedDatabase, Parameters, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::ion_series::Kind;
use sage_core::lfq::{Peak, PrecursorId};
use sage_core::mass::Tolerance;
use sage_core::peptide::Peptide;
use sage_core::scoring::Fragments;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectrum::{MS1Spectra, ProcessedSpectrum, RawSpectrum, SpectrumProcessor};
use sage_core::tmt::TmtQuant;
use std::collections::{HashMap, HashSet};
use std::time::Instant;
// HTML report specific imports
use maud::{html, PreEscaped};
use report_builder::{
    plots::{plot_boxplot, plot_pp, plot_scatter, plot_score_histogram},
    Report, ReportSection,
};

pub struct Runner {
    pub database: IndexedDatabase,
    pub parameters: Search,
    start: Instant,
}

#[derive(Default)]
struct RawSpectrumAccumulator {
    pub ms1: Vec<RawSpectrum>,
    pub msn: Vec<RawSpectrum>,
}

impl RawSpectrumAccumulator {
    pub fn fold_op(mut self, rhs: RawSpectrum) -> Self {
        if rhs.ms_level == 1 {
            self.ms1.push(rhs);
        } else {
            self.msn.push(rhs);
        }
        self
    }

    pub fn reduce(mut self, other: Self) -> Self {
        self.ms1.extend(other.ms1);
        self.msn.extend(other.msn);
        self
    }
}

impl FromParallelIterator<RawSpectrum> for RawSpectrumAccumulator {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = RawSpectrum>,
    {
        par_iter
            .into_par_iter()
            .fold(
                RawSpectrumAccumulator::default,
                RawSpectrumAccumulator::fold_op,
            )
            .reduce(
                RawSpectrumAccumulator::default,
                RawSpectrumAccumulator::reduce,
            )
    }
}

impl FromIterator<RawSpectrum> for RawSpectrumAccumulator {
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = RawSpectrum>,
    {
        iter.into_iter().fold(
            RawSpectrumAccumulator::default(),
            RawSpectrumAccumulator::fold_op,
        )
    }
}

impl Runner {
    pub fn new(parameters: Search, parallel: usize) -> anyhow::Result<Self> {
        let mut parameters = parameters.clone();
        let start = Instant::now();

        // Check if we should load a pre-built index
        let database = if let Some(ref load_path) = parameters.load_index {
            info!("Loading pre-built index from {}", load_path);
            let path = load_path.parse::<CloudPath>()?;
            let loaded_db = sage_cloudpath::index_parquet::deserialize_index(&path)
                .with_context(|| format!("Failed to load index from `{}`", load_path))?;

            info!(
                "loaded {} fragments, {} peptides in {:#?}",
                loaded_db.fragments.len(),
                loaded_db.peptides.len(),
                start.elapsed()
            );

            // Optional validation: build from FASTA and compare
            if parameters.validate_index {
                info!("Validating loaded index against FASTA build...");
                let fasta = sage_cloudpath::util::read_fasta(
                    &parameters.database.fasta,
                    &parameters.database.decoy_tag,
                    parameters.database.generate_decoys,
                )
                .with_context(|| {
                    format!("Failed to read FASTA from `{}`", parameters.database.fasta)
                })?;
                let built_db = parameters.database.clone().build(fasta);
                sage_cloudpath::index_parquet::validate_index(&built_db, &loaded_db)
                    .with_context(|| "Index validation failed")?;
                info!("Index validation passed!");
            }

            loaded_db
        } else {
            // Build from FASTA as usual
            let fasta = sage_cloudpath::util::read_fasta(
                &parameters.database.fasta,
                &parameters.database.decoy_tag,
                parameters.database.generate_decoys,
            )
            .with_context(|| {
                format!(
                    "Failed to build database from `{}`",
                    parameters.database.fasta
                )
            })?;

            let built_db = match parameters.database.prefilter {
                false => parameters.database.clone().build(fasta),
                true => {
                    parameters
                        .database
                        .auto_calculate_prefilter_chunk_size(&fasta);
                    if parameters.database.prefilter_chunk_size >= fasta.targets.len() {
                        parameters.database.clone().build(fasta)
                    } else {
                        info!(
                            "using {} db chunks of size {}",
                            (fasta.targets.len() + parameters.database.prefilter_chunk_size - 1)
                                / parameters.database.prefilter_chunk_size,
                            parameters.database.prefilter_chunk_size,
                        );
                        let mini_runner = Self {
                            database: IndexedDatabase::default(),
                            parameters: parameters.clone(),
                            start,
                        };
                        let peptides = mini_runner.prefilter_peptides(parallel, fasta);
                        parameters.database.clone().build_from_peptides(peptides)
                    }
                }
            };

            info!(
                "generated {} fragments, {} peptides in {:#?}",
                built_db.fragments.len(),
                built_db.peptides.len(),
                start.elapsed()
            );

            // Save index if requested
            if let Some(ref save_path) = parameters.save_index {
                info!("Saving index to {}", save_path);
                let path = save_path.parse::<CloudPath>()?;
                sage_cloudpath::index_parquet::serialize_index(&built_db, &path)
                    .with_context(|| format!("Failed to save index to `{}`", save_path))?;
                info!("Index saved successfully");
            }

            built_db
        };

        Ok(Self {
            database,
            parameters,
            start,
        })
    }

    pub fn prefilter_peptides(self, parallel: usize, fasta: Fasta) -> Vec<Peptide> {
        let spectra: Option<Vec<ProcessedSpectrum<_>>> =
            match parallel >= self.parameters.mzml_paths.len() {
                true => Some(
                    self.read_processed_spectra(&self.parameters.mzml_paths, 0, 0)
                        .1,
                ),
                false => None,
            };

        let mut db_params = self.parameters.database.clone();
        // TODO: Don't generate decoys for fast searching
        // * if `generate_decoys` is used, we should re-generate at the end
        //  to ensure that picked-peptide conditions are used, otherwise,
        //  if the user supplied decoys in the fasta file, then we should retain them
        //
        // db_params.generate_decoys = false;

        let mut all_peptides: Vec<Peptide> = fasta
            .iter_chunks(self.parameters.database.prefilter_chunk_size)
            .enumerate()
            .flat_map(|(chunk_id, fasta_chunk)| {
                let start = Instant::now();
                info!("pre-filtering fasta chunk {}", chunk_id,);
                let mut db = db_params.clone().build(fasta_chunk);

                info!(
                    "generated {} fragments, {} peptides in {}ms",
                    db.fragments.len(),
                    db.peptides.len(),
                    (Instant::now() - start).as_millis()
                );

                let scorer = Scorer {
                    db: &db,
                    precursor_tol: self.parameters.precursor_tol,
                    fragment_tol: self.parameters.fragment_tol,
                    min_matched_peaks: self.parameters.min_matched_peaks,
                    min_isotope_err: self.parameters.isotope_errors.0,
                    max_isotope_err: self.parameters.isotope_errors.1,
                    min_precursor_charge: self.parameters.precursor_charge.0,
                    max_precursor_charge: self.parameters.precursor_charge.1,
                    override_precursor_charge: self.parameters.override_precursor_charge,
                    max_fragment_charge: self.parameters.max_fragment_charge,
                    chimera: self.parameters.chimera,
                    report_psms: self.parameters.report_psms + 1, // Q: Why is 1 being added here? (JSPP: Feb 2024)
                    wide_window: self.parameters.wide_window,
                    annotate_matches: self.parameters.annotate_matches,
                    score_type: self.parameters.score_type,
                };

                // Allocate an array of booleans indicating whether a peptide was identified in a
                // preliminary pass of the data
                let keep = (0..db.peptides.len())
                    .map(|_| std::sync::atomic::AtomicBool::new(false))
                    .collect::<Vec<_>>();

                match &spectra {
                    Some(spectra) => {
                        self.peptide_filter_processed_spectra(&scorer, &spectra, &keep)
                    }
                    None => self
                        .parameters
                        .mzml_paths
                        .chunks(parallel)
                        .enumerate()
                        .for_each(|(chunk_idx, chunk)| {
                            let spectra_chunk =
                                self.read_processed_spectra(chunk, chunk_idx, parallel).1;
                            self.peptide_filter_processed_spectra(&scorer, &spectra_chunk, &keep)
                        }),
                };

                // Retain only peptides where `keep[ix] = true`
                let peptides = db
                    .peptides
                    .drain(..)
                    .enumerate()
                    .filter_map(|(ix, peptide)| {
                        let val = keep[ix].load(std::sync::atomic::Ordering::Relaxed);
                        if val {
                            Some(peptide)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>();

                info!(
                    "found {} pre-filtered peptides for fasta chunk {}",
                    peptides.len(),
                    chunk_id,
                );
                peptides
            })
            .collect();

        Parameters::reorder_peptides(&mut all_peptides);
        all_peptides
    }

    fn peptide_filter_processed_spectra(
        &self,
        scorer: &Scorer,
        spectra: &Vec<ProcessedSpectrum<sage_core::spectrum::Peak>>,
        keep: &[std::sync::atomic::AtomicBool],
    ) {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= self.parameters.min_peaks && spec.level == 2)
            .for_each(|spectrum| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev % 10_000 == 0 {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                scorer.quick_score(
                    spectrum,
                    self.parameters.database.prefilter_low_memory,
                    keep,
                )
            });

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!(
            "- prefilter search:  {:8} ms ({} spectra/s)",
            duration,
            rate
        );
    }

    fn spectrum_fdr(&self, features: &mut [Feature]) -> usize {
        if sage_core::ml::linear_discriminant::score_psms(features, self.parameters.precursor_tol)
            .is_none()
        {
            log::warn!("linear model fitting failed, falling back to heuristic discriminant score");
            features.par_iter_mut().for_each(|feat| {
                feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
            });
        }
        features.par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));
        sage_core::ml::qvalue::spectrum_q_value(features)
    }

    // Create a path for `file_name` in the specified output directory, if it exists,
    // otherwise, write to current directory
    fn make_path<S: AsRef<str>>(&self, file_name: S) -> CloudPath {
        let mut path = self.parameters.output_directory.clone();
        path.push(file_name);
        path
    }

    fn search_processed_spectra(
        &self,
        scorer: &Scorer,
        msn_spectra: &Vec<ProcessedSpectrum<sage_core::spectrum::Peak>>,
    ) -> Vec<Feature> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        let features: Vec<_> = msn_spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= self.parameters.min_peaks && spec.level == 2)
            .map(|x| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev % 10_000 == 0 {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                x
            })
            .flat_map(|spec| scorer.score(spec))
            .collect();

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!("- search:  {:8} ms ({} spectra/s)", duration, rate);
        features
    }

    fn complete_features(
        &self,
        msn_spectra: Vec<ProcessedSpectrum<sage_core::spectrum::Peak>>,
        ms1_spectra: MS1Spectra,
        features: Vec<Feature>,
    ) -> SageResults {
        let quant = self
            .parameters
            .quant
            .tmt
            .as_ref()
            .map(|isobaric| {
                let level = self.parameters.quant.tmt_settings.level;
                if level != 2 && level != 3 {
                    log::warn!("TMT quant level set at {}, is this correct?", level);
                }
                sage_core::tmt::quantify(&msn_spectra, isobaric, Tolerance::Ppm(-20.0, 20.0), level)
            })
            .unwrap_or_default();

        SageResults {
            features,
            quant,
            ms1: ms1_spectra,
        }
    }

    fn requires_ms1(&self) -> bool {
        self.parameters.quant.lfq
    }

    fn process_chunk(
        &self,
        scorer: &Scorer,
        chunk: &[String],
        chunk_idx: usize,
        batch_size: usize,
    ) -> SageResults {
        let spectra = self.read_processed_spectra(chunk, chunk_idx, batch_size);
        let features = self.search_processed_spectra(scorer, &spectra.1);
        self.complete_features(spectra.1, spectra.0, features)
    }

    fn read_processed_spectra(
        &self,
        chunk: &[String],
        chunk_idx: usize,
        batch_size: usize,
    ) -> (
        MS1Spectra,
        Vec<ProcessedSpectrum<sage_core::spectrum::Peak>>,
    ) {
        // Read all of the spectra at once - this can help prevent memory over-consumption issues
        info!(
            "processing files {} .. {} ",
            batch_size * chunk_idx,
            batch_size * chunk_idx + chunk.len()
        );
        let start = Instant::now();

        let sn = self
            .parameters
            .quant
            .tmt_settings
            .sn
            .then_some(self.parameters.quant.tmt_settings.level);

        let min_deisotope_mz = match &self.parameters.quant.tmt {
            Some(i) => match self.parameters.quant.tmt_settings.level {
                2 => i.reporter_masses().last().map(|x| x * (1.0 + 20E-6)),
                _ => None,
            },
            None => None,
        };

        let sp = SpectrumProcessor::new(
            self.parameters.max_peaks,
            self.parameters.deisotope,
            min_deisotope_mz.unwrap_or(0.0),
        );

        // If the file format supports parallel reading, then we can read
        // then it is faster to read each file in series. (since each spectra
        // will be processed internally in parallel).
        let file_serial_read = chunk
            .iter()
            .all(|path| FileFormat::from(path.as_ref()).within_file_parallel());
        log::trace!("file serial read: {}", file_serial_read);
        let inner_closure = |(idx, path)| {
            let file_id = chunk_idx * batch_size + idx;
            let res = sage_cloudpath::util::read_spectra(
                path,
                file_id,
                sn,
                self.parameters.bruker_config.clone(),
                self.requires_ms1(),
            );

            match res {
                Ok(s) => {
                    log::trace!("- {}: read {} spectra", path, s.len());
                    Ok(s)
                }
                Err(e) => {
                    log::error!("- {}: {}", path, e);
                    Err(e)
                }
            }
        };

        let spectra: RawSpectrumAccumulator = if file_serial_read {
            chunk
                .iter()
                .enumerate()
                .flat_map(inner_closure)
                .flatten()
                .collect()
        } else {
            chunk
                .par_iter()
                .enumerate()
                .flat_map(inner_closure)
                .flatten()
                .collect()
        };

        let msn_spectra = spectra
            .msn
            .into_par_iter()
            .map(|s| sp.process(s))
            .collect::<Vec<_>>();

        // If all the MS1 spectra contain IMS, then we can process them
        // we use the IMS! otherwise we dont.
        // Note: Empty iterators return true.
        let all_contain_ims = spectra.ms1.iter().all(|x| x.mobility.is_some());
        let ms1_empty = spectra.ms1.is_empty();
        let ms1_spectra = if ms1_empty {
            log::trace!("no MS1 spectra found");
            MS1Spectra::Empty
        } else if all_contain_ims {
            log::trace!("Processing MS1 spectra with IMS");
            let spectra = spectra
                .ms1
                .into_iter()
                .map(|x| sp.process_with_mobility(x))
                .collect();
            MS1Spectra::WithMobility(spectra)
        } else {
            log::trace!("Processing MS1 spectra without IMS");
            let spectra = spectra.ms1.into_iter().map(|s| sp.process(s)).collect();
            MS1Spectra::NoMobility(spectra)
        };

        let io_time = Instant::now() - start;
        info!("- file IO: {:8} ms", io_time.as_millis());

        (ms1_spectra, msn_spectra)
    }

    pub fn batch_files(&self, scorer: &Scorer, batch_size: usize) -> SageResults {
        self.parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| self.process_chunk(scorer, chunk, chunk_idx, batch_size))
            .collect::<SageResults>()
    }

    pub fn run(mut self, parallel: usize, parquet: bool) -> anyhow::Result<telemetry::Telemetry> {
        let scorer = Scorer {
            db: &self.database,
            precursor_tol: self.parameters.precursor_tol,
            fragment_tol: self.parameters.fragment_tol,
            min_matched_peaks: self.parameters.min_matched_peaks,
            min_isotope_err: self.parameters.isotope_errors.0,
            max_isotope_err: self.parameters.isotope_errors.1,
            min_precursor_charge: self.parameters.precursor_charge.0,
            max_precursor_charge: self.parameters.precursor_charge.1,
            override_precursor_charge: self.parameters.override_precursor_charge,
            max_fragment_charge: self.parameters.max_fragment_charge,
            chimera: self.parameters.chimera,
            report_psms: self.parameters.report_psms,
            wide_window: self.parameters.wide_window,
            annotate_matches: self.parameters.annotate_matches,
            score_type: self.parameters.score_type,
        };

        //Collect all results into a single container
        let mut outputs = self.batch_files(&scorer, parallel);

        let alignments = if self.parameters.predict_rt {
            // Poisson probability is usually the best single feature for refining FDR.
            // Take our set of 1% FDR filtered PSMs, and use them to train a linear
            // regression model for predicting retention time
            outputs
                .features
                .par_sort_unstable_by(|a, b| a.poisson.total_cmp(&b.poisson));
            sage_core::ml::qvalue::spectrum_q_value(&mut outputs.features);

            let alignments = sage_core::ml::retention_alignment::global_alignment(
                &mut outputs.features,
                self.parameters.mzml_paths.len(),
            );
            let _ = sage_core::ml::retention_model::predict(&self.database, &mut outputs.features);
            let _ = sage_core::ml::mobility_model::predict(&self.database, &mut outputs.features);
            Some(alignments)
        } else {
            None
        };

        let q_spectrum = self.spectrum_fdr(&mut outputs.features);
        let q_peptide = sage_core::fdr::picked_peptide(&self.database, &mut outputs.features);
        let q_protein = sage_core::fdr::picked_protein(&self.database, &mut outputs.features);

        let filenames = self
            .parameters
            .mzml_paths
            .iter()
            .map(|s| {
                s.parse::<CloudPath>()
                    .ok()
                    .and_then(|c| c.filename().map(|s| s.to_string()))
                    .unwrap_or_else(|| s.clone())
            })
            .collect::<Vec<_>>();

        let areas = alignments.and_then(|alignments| {
            if self.parameters.quant.lfq {
                log::trace!("performing LFQ");
                let mut areas = sage_core::lfq::build_feature_map(
                    self.parameters.quant.lfq_settings,
                    self.parameters.precursor_charge,
                    &outputs.features,
                )
                .quantify(&self.database, &outputs.ms1, &alignments);

                let q_precursor = sage_core::fdr::picked_precursor(&mut areas);

                log::info!("discovered {} target MS1 peaks at 5% FDR", q_precursor);
                Some(areas)
            } else {
                None
            }
        });

        log::info!(
            "discovered {} target peptide-spectrum matches at 1% FDR",
            q_spectrum
        );
        log::info!("discovered {} target peptides at 1% FDR", q_peptide);
        log::info!("discovered {} target proteins at 1% FDR", q_protein);
        log::trace!("writing outputs");

        // Write either a single parquet file, or multiple tsv files
        if parquet {
            log::warn!("parquet output format is currently unstable! There may be failures or schema changes!");

            let bytes = sage_cloudpath::parquet::serialize_features(
                &outputs.features,
                &outputs.quant,
                &filenames,
                &self.database,
            )?;

            let path = self.make_path("results.sage.parquet");
            path.write_bytes_sync(bytes)?;
            self.parameters.output_paths.push(path.to_string());

            if self.parameters.annotate_matches {
                let bytes =
                    sage_cloudpath::parquet::serialize_matched_fragments(&outputs.features)?;
                let path = self.make_path("matched_fragments.sage.parquet");
                path.write_bytes_sync(bytes)?;
                self.parameters.output_paths.push(path.to_string());
            }

            if let Some(areas) = &areas {
                let bytes =
                    sage_cloudpath::parquet::serialize_lfq(areas, &filenames, &self.database)?;

                let path = self.make_path("lfq.parquet");
                path.write_bytes_sync(bytes)?;
                self.parameters.output_paths.push(path.to_string());
            }
        } else {
            self.parameters
                .output_paths
                .push(self.write_features(&outputs.features, &filenames)?);

            if self.parameters.annotate_matches {
                self.parameters
                    .output_paths
                    .push(self.write_fragments(&outputs.features)?);
            }

            if !outputs.quant.is_empty() {
                self.parameters
                    .output_paths
                    .push(self.write_tmt(&outputs.quant, &filenames)?);
            }
            if let Some(areas) = areas.clone() {
                self.parameters
                    .output_paths
                    .push(self.write_lfq(areas, &filenames)?);
            }
        }

        // Write percolator input file if requested
        if self.parameters.write_pin {
            self.parameters
                .output_paths
                .push(self.write_pin(&outputs.features, &filenames)?);
        }

        // Write an html report if requested
        if self.parameters.write_report {
            self.parameters.output_paths.push(self.write_report(
                &outputs.features,
                areas,
                &filenames,
            )?);
        }

        let path = self.make_path("results.json");
        self.parameters.output_paths.push(path.to_string());
        println!("{}", serde_json::to_string_pretty(&self.parameters)?);

        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        path.write_bytes_sync(bytes)?;

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);
        info!("cite: \"Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale\" https://doi.org/10.1021/acs.jproteome.3c00486");

        let telemetry = telemetry::Telemetry::new(
            self.parameters,
            self.database.peptides.len(),
            self.database.fragments.len(),
            parquet,
            run_time,
        );

        Ok(telemetry)
    }
    pub fn serialize_feature(&self, feature: &Feature, filenames: &[String]) -> csv::ByteRecord {
        let mut record = csv::ByteRecord::new();

        record.push_field(itoa::Buffer::new().format(feature.psm_id).as_bytes());

        let peptide = &self.database[feature.peptide_idx];
        record.push_field(peptide.to_string().as_bytes());
        record.push_field(
            peptide
                .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(peptide.proteins.len())
                .as_bytes(),
        );
        record.push_field(filenames[feature.file_id].as_bytes());
        record.push_field(feature.spec_id.as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.rank).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.label).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.expmass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.calcmass).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.charge).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.peptide_len).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format(feature.missed_cleavages)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(peptide.semi_enzymatic as u8)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.isotope_error).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_mass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.hyperscore).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_next).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_best).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.aligned_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_rt_model).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ims).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_ims).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_ims_model)
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.matched_peaks).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_b).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_y).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.longest_y_pct).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.matched_intensity_pct)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.scored_candidates)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.poisson).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.discriminant_score)
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.posterior_error)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.spectrum_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.peptide_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.protein_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ms2_intensity).as_bytes());
        record
    }

    pub fn serialize_fragments(
        &self,
        psm_id: usize,
        fragments_: &Option<Fragments>,
    ) -> Vec<ByteRecord> {
        let mut frag_records = vec![];

        if let Some(fragments) = fragments_ {
            for id in 0..fragments.fragment_ordinals.len() {
                let mut record = ByteRecord::new();
                record.push_field(itoa::Buffer::new().format(psm_id).as_bytes());
                let ion_type = match fragments.kinds[id] {
                    Kind::A => "a",
                    Kind::B => "b",
                    Kind::C => "c",
                    Kind::X => "x",
                    Kind::Y => "y",
                    Kind::Z => "z",
                };
                record.push_field(ion_type.as_bytes());
                record.push_field(
                    itoa::Buffer::new()
                        .format(fragments.fragment_ordinals[id])
                        .as_bytes(),
                );
                record.push_field(itoa::Buffer::new().format(fragments.charges[id]).as_bytes());
                record.push_field(
                    ryu::Buffer::new()
                        .format(fragments.mz_calculated[id])
                        .as_bytes(),
                );
                record.push_field(
                    ryu::Buffer::new()
                        .format(fragments.mz_experimental[id])
                        .as_bytes(),
                );
                record.push_field(
                    ryu::Buffer::new()
                        .format(fragments.intensities[id])
                        .as_bytes(),
                );
                frag_records.push(record);
            }
        }

        frag_records
    }

    pub fn write_features(
        &self,
        features: &[Feature],
        filenames: &[String],
    ) -> anyhow::Result<String> {
        let path = self.make_path("results.sage.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        let csv_headers = vec![
            "psm_id",
            "peptide",
            "proteins",
            "num_proteins",
            "filename",
            "scannr",
            "rank",
            "label",
            "expmass",
            "calcmass",
            "charge",
            "peptide_len",
            "missed_cleavages",
            "semi_enzymatic",
            "isotope_error",
            "precursor_ppm",
            "fragment_ppm",
            "hyperscore",
            "delta_next",
            "delta_best",
            "rt",
            "aligned_rt",
            "predicted_rt",
            "delta_rt_model",
            "ion_mobility",
            "predicted_mobility",
            "delta_mobility",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "matched_intensity_pct",
            "scored_candidates",
            "poisson",
            "sage_discriminant_score",
            "posterior_error",
            "spectrum_q",
            "peptide_q",
            "protein_q",
            "ms2_intensity",
        ];

        let headers = csv::ByteRecord::from(csv_headers);

        wtr.write_byte_record(&headers)?;
        for record in features
            .into_par_iter()
            .map(|feat| self.serialize_feature(feat, filenames))
            .collect::<Vec<_>>()
        {
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    pub fn write_fragments(&self, features: &[Feature]) -> anyhow::Result<String> {
        let path = self.make_path("matched_fragments.sage.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        let headers = csv::ByteRecord::from(vec![
            "psm_id",
            "fragment_type",
            "fragment_ordinals",
            "fragment_charge",
            "fragment_mz_calculated",
            "fragment_mz_experimental",
            "fragment_intensity",
        ]);

        wtr.write_byte_record(&headers)?;

        for record in features
            .into_par_iter()
            .map(|feat| self.serialize_fragments(feat.psm_id, &feat.fragments))
            .flatten()
            .collect::<Vec<_>>()
        {
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    fn serialize_pin(
        &self,
        re: &regex::Regex,
        feature: &Feature,
        filenames: &[String],
    ) -> csv::ByteRecord {
        let scannr = re
            .captures_iter(&feature.spec_id)
            .last()
            .and_then(|cap| cap.get(1).map(|cap| cap.as_str()))
            .unwrap_or(&feature.spec_id);

        let mut record = csv::ByteRecord::new();
        let peptide = &self.database[feature.peptide_idx];
        record.push_field(itoa::Buffer::new().format(feature.psm_id).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.label).as_bytes());
        record.push_field(scannr.as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.expmass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.calcmass).as_bytes());
        record.push_field(filenames[feature.file_id].as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ims).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.rank).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 2) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 3) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 4) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 5) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 6) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(if feature.charge < 2 || feature.charge > 6 {
                    feature.charge
                } else {
                    0
                })
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.peptide_len).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format(feature.missed_cleavages)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(peptide.semi_enzymatic as u8)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.isotope_error).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_mass.abs().ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.hyperscore.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_next.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_best.ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.aligned_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_rt_model.clamp(0.001, 1.0).sqrt())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.predicted_ims).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_ims_model)
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.matched_peaks).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_b).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_y).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.longest_y_pct).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.matched_intensity_pct.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.scored_candidates)
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format((-feature.poisson).ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.posterior_error)
                .as_bytes(),
        );
        record.push_field(peptide.to_string().as_bytes());
        record.push_field(
            peptide
                .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                .as_bytes(),
        );
        record
    }

    pub fn write_pin(&self, features: &[Feature], filenames: &[String]) -> anyhow::Result<String> {
        let path = self.make_path("results.sage.pin");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        let headers = csv::ByteRecord::from(vec![
            "SpecId",
            "Label",
            "ScanNr",
            "ExpMass",
            "CalcMass",
            "FileName",
            "retentiontime",
            "ion_mobility",
            "rank",
            "z=2",
            "z=3",
            "z=4",
            "z=5",
            "z=6",
            "z=other",
            "peptide_len",
            "missed_cleavages",
            "semi_enzymatic",
            "isotope_error",
            "ln(precursor_ppm)",
            "fragment_ppm",
            "ln(hyperscore)",
            "ln(delta_next)",
            "ln(delta_best)",
            "aligned_rt",
            "predicted_rt",
            "sqrt(delta_rt_model)",
            "predicted_mobility",
            "sqrt(delta_mobility)",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "ln(matched_intensity_pct)",
            "scored_candidates",
            "ln(-poisson)",
            "posterior_error",
            "Peptide",
            "Proteins",
        ]);

        let re = regex::Regex::new(r"scan=(\d+)").expect("This is valid regex");

        wtr.write_byte_record(&headers)?;
        for record in features
            .into_par_iter()
            .map(|feat| self.serialize_pin(&re, feat, filenames))
            .collect::<Vec<_>>()
        {
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    pub fn write_tmt(&self, quant: &[TmtQuant], filenames: &[String]) -> anyhow::Result<String> {
        let path = self.make_path("tmt.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = csv::ByteRecord::from(vec!["filename", "scannr", "ion_injection_time"]);
        headers.extend(
            self.parameters
                .quant
                .tmt
                .as_ref()
                .map(|tmt| tmt.headers())
                .expect("TMT quant cannot be performed without setting this parameter"),
        );

        wtr.write_byte_record(&headers)?;

        let records = quant
            .into_par_iter()
            .map(|q| {
                let mut record = csv::ByteRecord::new();
                record.push_field(filenames[q.file_id].as_bytes());
                record.push_field(q.spec_id.as_bytes());
                record.push_field(ryu::Buffer::new().format(q.ion_injection_time).as_bytes());
                for peak in &q.peaks {
                    record.push_field(ryu::Buffer::new().format(*peak).as_bytes());
                }
                record
            })
            .collect::<Vec<csv::ByteRecord>>();

        for record in records {
            wtr.write_record(&record)?;
        }
        wtr.flush()?;

        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    pub fn write_lfq(
        &self,
        areas: HashMap<(PrecursorId, bool), (Peak, Vec<f64>), fnv::FnvBuildHasher>,
        filenames: &[String],
    ) -> anyhow::Result<String> {
        let path = self.make_path("lfq.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = csv::ByteRecord::from(vec![
            "peptide",
            "charge",
            "proteins",
            "q_value",
            "score",
            "spectral_angle",
        ]);
        headers.extend(filenames);

        wtr.write_byte_record(&headers)?;

        let records = areas
            .into_par_iter()
            .filter_map(|((id, decoy), (peak, data))| {
                if decoy {
                    return None;
                };
                let mut record = csv::ByteRecord::new();
                let (peptide_ix, charge) = match id {
                    PrecursorId::Combined(x) => (x, None),
                    PrecursorId::Charged((x, charge)) => (x, Some(charge as i32)),
                };
                record.push_field(self.database[peptide_ix].to_string().as_bytes());
                record.push_field(itoa::Buffer::new().format(charge.unwrap_or(-1)).as_bytes());
                record.push_field(
                    self.database[peptide_ix]
                        .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                        .as_bytes(),
                );
                record.push_field(ryu::Buffer::new().format(peak.q_value).as_bytes());
                record.push_field(ryu::Buffer::new().format(peak.score).as_bytes());
                record.push_field(ryu::Buffer::new().format(peak.spectral_angle).as_bytes());
                for x in data {
                    record.push_field(ryu::Buffer::new().format(x).as_bytes());
                }
                Some(record)
            })
            .collect::<Vec<csv::ByteRecord>>();

        for record in records {
            wtr.write_record(&record)?;
        }
        wtr.flush()?;

        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    fn write_report(
        &self,
        features: &[Feature],
        areas: Option<HashMap<(PrecursorId, bool), (Peak, Vec<f64>), fnv::FnvBuildHasher>>,
        filenames: &[String],
    ) -> anyhow::Result<String> {
        let path = self.make_path("results.sage.report.html");

        let global_q_value_filter = 0.01;
        let predict_section_q_value_filter = 0.01;

        // Create a new report
        let mut report = Report::new(
            "Sage",
            &self.parameters.version,
            Some("https://github.com/lazear/sage/blob/master/figures/logo.png?raw=true"),
            "Sage Report",
        );

        /* Section 1: Introduction */
        {
            let mut intro_section = ReportSection::new("Results Overview");
            intro_section.add_content(html! {
                "The following files were processed:"
                ul {
                    @for filename in filenames {
                        li { (filename) }
                    }
                }
            });

            // Number of targets identified at global q-value filter at spectrum level per file
            let num_psm_targets_per_file: Vec<usize> = filenames
                .iter()
                .map(|filename| {
                    features
                        .iter()
                        .filter(|f| {
                            f.label == 1
                                && f.spectrum_q <= global_q_value_filter
                                && filenames[f.file_id] == filename.to_string()
                        })
                        .count()
                })
                .collect();

            // Number of peptides identified at global q-value filter at peptide level per file
            let mut num_peptide_targets_per_file: Vec<usize> = Vec::new();
            for filename in filenames {
                let mut peptides = HashSet::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.peptide_q <= global_q_value_filter
                        && filenames[f.file_id] == filename.to_string()
                }) {
                    peptides.insert(self.database[feature.peptide_idx].to_string());
                }
                num_peptide_targets_per_file.push(peptides.len());
            }

            // Number of proteins identified at global q-value filter at protein level per file
            let mut num_protein_targets_per_file: Vec<usize> = Vec::new();
            for filename in filenames {
                let mut proteins = HashSet::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.protein_q <= global_q_value_filter
                        && filenames[f.file_id] == filename.to_string()
                }) {
                    proteins.insert(
                        self.database[feature.peptide_idx]
                            .proteins(&self.database.decoy_tag, self.database.generate_decoys),
                    );
                }
                num_protein_targets_per_file.push(proteins.len());
            }

            // Total MS2 intensity at global q-value filter at each level per file
            let total_ms2_intensity_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    features
                        .iter()
                        .filter(|f| {
                            f.label == 1
                                && f.spectrum_q <= global_q_value_filter
                                && f.peptide_q <= global_q_value_filter
                                && f.protein_q <= global_q_value_filter
                                && filenames[f.file_id] == filename.to_string()
                        })
                        .map(|f| f.ms2_intensity)
                        .sum()
                })
                .collect();

            // Total LFQ (MS1) intensity at global q-value filter per file (if LFQ is enabled)
            let total_lfq_intensity_per_file: Vec<f32> = if let Some(areas) = &areas {
                let mut total_lfq_intensities = Vec::new();
                for i in 0..filenames.len() {
                    let mut intensities = Vec::new();
                    for ((id, decoy), (peak, data)) in areas {
                        if !decoy && peak.q_value <= global_q_value_filter {
                            intensities.push(data[i] as f32);
                        }
                    }
                    total_lfq_intensities.push(intensities.iter().sum());
                }
                total_lfq_intensities
            } else {
                vec![0.0; filenames.len()]
            };

            // Mmedian MS1 mass accuracy for each file, using feature.delta_mass
            let median_ms1_mass_accuracy_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut accuracies = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == filename.to_string()
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        accuracies.push(feature.delta_mass);
                    }
                    accuracies.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = accuracies.len() / 2;

                    if accuracies.is_empty() {
                        return std::f32::NAN;
                    }

                    if accuracies.len() % 2 == 0 {
                        if mid > 0 {
                            (accuracies[mid - 1] + accuracies[mid]) / 2.0
                        } else {
                            accuracies[mid]
                        }
                    } else {
                        accuracies[mid]
                    }
                })
                .collect();

            // Median MS2 mass accuracy for each file, using feature.average_ppm
            let median_ms2_mass_accuracy_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut accuracies = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == filename.to_string()
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        accuracies.push(feature.average_ppm);
                    }
                    accuracies.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = accuracies.len() / 2;

                    if accuracies.is_empty() {
                        return std::f32::NAN;
                    }

                    if accuracies.len() % 2 == 0 {
                        if mid > 0 {
                            (accuracies[mid - 1] + accuracies[mid]) / 2.0
                        } else {
                            accuracies[mid]
                        }
                    } else {
                        accuracies[mid]
                    }
                })
                .collect();

            // Median RT deviation for each file, using feature.delta_rt_model
            let median_rt_deviation_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut deviations = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == filename.to_string()
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        deviations.push(feature.delta_rt_model);
                    }
                    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = deviations.len() / 2;

                    if deviations.is_empty() {
                        return std::f32::NAN;
                    }

                    if deviations.len() % 2 == 0 {
                        if mid > 0 {
                            (deviations[mid - 1] + deviations[mid]) / 2.0
                        } else {
                            deviations[mid]
                        }
                    } else {
                        deviations[mid]
                    }
                })
                .collect();

            // Median IM deviation for each file, using feature.delta_ims_model
            let median_im_deviation_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut deviations = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == filename.to_string()
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        deviations.push(feature.delta_ims_model);
                    }
                    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let mid = deviations.len() / 2;

                    if deviations.is_empty() {
                        return std::f32::NAN;
                    }

                    if deviations.len() % 2 == 0 {
                        if mid > 0 {
                            (deviations[mid - 1] + deviations[mid]) / 2.0
                        } else {
                            deviations[mid]
                        }
                    } else {
                        deviations[mid]
                    }
                })
                .collect();

            // Average peptide length for each file
            let avg_peptide_length_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut lengths = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == filename.to_string()
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        lengths.push(feature.peptide_len as f32);
                    }
                    lengths.iter().sum::<f32>() / lengths.len() as f32
                })
                .collect();

            // Average peptide charge for each file
            let avg_peptide_charge_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut charges = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == filename.to_string()
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        charges.push(feature.charge as f32);
                    }
                    charges.iter().sum::<f32>() / charges.len() as f32
                })
                .collect();

            // Average number of matched peaks for each file
            let avg_matched_peaks_per_file: Vec<f32> = filenames
                .iter()
                .map(|filename| {
                    let mut peaks = Vec::new();
                    for feature in features.iter().filter(|f| {
                        filenames[f.file_id] == filename.to_string()
                            && f.label == 1
                            && f.spectrum_q <= global_q_value_filter
                    }) {
                        peaks.push(feature.matched_peaks as f32);
                    }
                    peaks.iter().sum::<f32>() / peaks.len() as f32
                })
                .collect();

            // Prepare html table to add to the report
            let table = html! {
                div class="table-container" {
                    table id="dataTable"  class="display" {
                        thead {
                            tr {
                                th { "File" }
                                th { "PSMs" }
                                th { "Peptides" }
                                th { "Proteins" }
                                th { "Total MS1 Intensity" }
                                th { "Total MS2 Intensity" }
                                th { "Median MS1 Delta Mass" }
                                th { "Median MS2 Delta Mass" }
                                th { "Median RT Deviation" }
                                th { "Median IM Deviation" }
                                th { "Average Peptide Length" }
                                th { "Average Peptide Charge" }
                                th { "Average Matched Peaks" }
                            }
                        }
                        tbody {
                            @for (i, filename) in filenames.iter().enumerate() {
                                tr {
                                    td { (filename) }
                                    td { (num_psm_targets_per_file[i]) }
                                    td { (num_peptide_targets_per_file[i]) }
                                    td { (num_protein_targets_per_file[i]) }
                                    td { (total_lfq_intensity_per_file[i]) }
                                    td { (total_ms2_intensity_per_file[i]) }
                                    td { (median_ms1_mass_accuracy_per_file[i]) }
                                    td { (median_ms2_mass_accuracy_per_file[i]) }
                                    td { (median_rt_deviation_per_file[i]) }
                                    td { (median_im_deviation_per_file[i]) }
                                    td { (avg_peptide_length_per_file[i]) }
                                    td { (avg_peptide_charge_per_file[i]) }
                                    td { (avg_matched_peaks_per_file[i]) }
                                }
                            }
                        }
                    }
                    button id="downloadCsv" { "Download as CSV" }
                }
            };

            intro_section.add_content(table);

            // Add boxplot of the LFQ intensities from areas if available
            if let Some(areas) = areas {
                let mut lfq_intensities: Vec<Vec<f64>> = Vec::new();
                for i in 0..filenames.len() {
                    let mut intensities = Vec::new();
                    for ((_id, decoy), (peak, data)) in &areas {
                        if !decoy && peak.q_value <= global_q_value_filter {
                            intensities.push(data[i].log2());
                        }
                    }
                    lfq_intensities.push(intensities);
                }

                let lfq_boxplot = plot_boxplot(
                    &lfq_intensities,
                    filenames.to_vec(),
                    &format!("LFQ Intensities ({:?}% Q-value)", global_q_value_filter),
                    "Run",
                    "Log2(Intensity)",
                )
                .unwrap();
                intro_section.add_plot(lfq_boxplot);
            }

            report.add_section(intro_section);
        }

        /* Section 2: Scoring QC */
        {
            let mut scoring_section = ReportSection::new("Scoring Quality Control");

            scoring_section.add_content(html! {
                "It is important to assess the quality of the scoring model to ensure that the model is performing as expected, and that we're not overfitting or violating any assumptions of the Target-Decoy approach. The plot below shows the distribution of discriminant scores for each PSM, colored by whether the PSM is a target or decoy. We would expect the target distributions to be bimodal, where the first mode represents false targets that should align with the decoy distribution, and the second mode represents true targets."
            });

            // Extract sage_discriminant_score and label from features
            let (scores, labels): (Vec<f64>, Vec<i32>) = features
                .iter()
                .map(|f| (f.discriminant_score as f64, f.label))
                .unzip();

            if !scores.is_empty() && scores.len() > 100 {
                let score_histogram =
                    plot_score_histogram(&scores, &labels, "LDA Score", "Score").unwrap();

                scoring_section.add_plot(score_histogram);

                let pp_plot = plot_pp(&scores, &labels, "PP Plot").unwrap();

                scoring_section.add_content(html! {
                    "The Probability-Probability (PP) plot is a diagnostic tool that can be used to assess the quality of the scoring model. It plots the empirical cumulative distribution function (ECDF) of the target distribution against the ECDF of the decoy distribution. See: Debrie, E. et. al. (2023) Journal of Proteome Research. for more information."
                });
                scoring_section.add_plot(pp_plot);

                let spectrum_q_histogram = plot_score_histogram(
                    &features
                        .iter()
                        .map(|f| f.spectrum_q as f64)
                        .collect::<Vec<f64>>(),
                    &labels,
                    "Spectrum Q-value",
                    "Q-value",
                )
                .unwrap();
                scoring_section.add_plot(spectrum_q_histogram);

                let peptide_q_histogram = plot_score_histogram(
                    &features
                        .iter()
                        .map(|f| f.peptide_q as f64)
                        .collect::<Vec<f64>>(),
                    &labels,
                    "Peptide Q-value",
                    "Q-value",
                )
                .unwrap();
                scoring_section.add_plot(peptide_q_histogram);

                let protein_q_histogram = plot_score_histogram(
                    &features
                        .iter()
                        .map(|f| f.protein_q as f64)
                        .collect::<Vec<f64>>(),
                    &labels,
                    "Protein Q-value",
                    "Q-value",
                )
                .unwrap();
                scoring_section.add_plot(protein_q_histogram);
            } else {
                scoring_section.add_content(html! {
                    div style="margin-top: 10px; margin-bottom: 10px; padding: 15px; background-color: #ffe6e6; border: 1px solid #ff9999; color: #cc0000; border-radius: 5px; white-space: pre-line;" {
                        p {
                            "There are not enough scores to plot the scoring quality control plots."
                        }
                    }
                });
            }

            report.add_section(scoring_section);
        }

        /* Section 3: Predicted Properties */
        {
            let mut predicted_properties_section = ReportSection::new("Predicted Properties");

            predicted_properties_section.add_content(html! {
                "The following plots show the predicted properties of target peptides. The plots show the predicted retention time and ion mobility if present. The predicted properties are used to assess the quality of the model and to identify potential outliers."
            });

            // Normalized experimental RT per file
            let mut rt_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut rts = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    rts.push(feature.rt as f64);
                }

                let min_rt = rts.iter().cloned().fold(f64::INFINITY, f64::min);
                let max_rt = rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                rts = rts
                    .iter()
                    .map(|rt| (rt - min_rt) / (max_rt - min_rt))
                    .collect();

                rt_per_file.push(rts);
            }

            // Predicted RT per file
            let mut predicted_rt_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut predicted_rts = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    predicted_rts.push(feature.predicted_rt as f64);
                }
                predicted_rt_per_file.push(predicted_rts);
            }

            let rt_scatter = plot_scatter(
                &rt_per_file,
                &predicted_rt_per_file,
                filenames.to_vec(),
                "Retention Time LR Model",
                "Retention Time",
                "Predicted Retention Time",
            )
            .unwrap();
            predicted_properties_section.add_plot(rt_scatter);

            // Experimental IMS per file
            let mut ims_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut imss = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    imss.push(feature.ims as f64);
                }

                ims_per_file.push(imss);
            }

            // Predicted IMS per file
            let mut predicted_ims_per_file: Vec<Vec<f64>> = Vec::new();
            for i in 0..filenames.len() {
                let mut predicted_imss = Vec::new();
                for feature in features.iter().filter(|f| {
                    f.label == 1
                        && f.spectrum_q <= predict_section_q_value_filter
                        && filenames[f.file_id] == filenames[i]
                }) {
                    predicted_imss.push(feature.predicted_ims as f64);
                }
                predicted_ims_per_file.push(predicted_imss);
            }

            if !ims_per_file.is_empty() && !predicted_ims_per_file.is_empty() {
                let ims_scatter = plot_scatter(
                    &ims_per_file,
                    &predicted_ims_per_file,
                    filenames.to_vec(),
                    "Ion Mobility LR Model",
                    "Ion Mobility",
                    "Predicted Ion Mobility",
                )
                .unwrap();
                predicted_properties_section.add_plot(ims_scatter);
            }

            report.add_section(predicted_properties_section);
        }

        /* Section 4: Configuration */
        {
            let mut config_section = ReportSection::new("Configuration");
            config_section.add_content(html! {
                style {
                    ".code-container {
                        background-color: #f5f5f5;
                        padding: 10px;
                        border-radius: 5px;
                        overflow-x: auto;
                        font-family: monospace;
                        white-space: pre-wrap;
                    }"
                }
                div class="code-container" {
                    pre {
                        code { (PreEscaped(serde_json::to_string_pretty(&self.parameters)?)) }
                    }
                }
            });
            report.add_section(config_section);
        }

        // Save the report to HTML file
        report.save_to_file(&path.to_string())?;

        Ok(path.to_string())
    }
}
