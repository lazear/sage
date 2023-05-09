use anyhow::Context;
use clap::{value_parser, Arg, Command, ValueHint};
use input::{Input, Search};
use log::info;
use rayon::prelude::*;
use sage_cloudpath::CloudPath;
use sage_core::database::IndexedDatabase;
use sage_core::mass::Tolerance;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectrum::{ProcessedSpectrum, SpectrumProcessor};
use sage_core::tmt::TmtQuant;
use std::time::Instant;

mod input;
mod output;

struct Runner {
    database: IndexedDatabase,
    parameters: input::Search,
    start: Instant,
}

#[derive(Debug)]
#[non_exhaustive]
enum RunnerError {
    ChunkProcessingError {
        chunk_idx: usize,
        source: Vec<RunnerError>,
    },
    FileProcessingError {
        source: anyhow::Error,
    },
}

impl std::fmt::Display for RunnerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RunnerError::ChunkProcessingError { chunk_idx, source } => {
                writeln!(f, "Processing chunk {chunk_idx} failed")?;
                for file_error in source {
                    writeln!(f, "\t{file_error}")?;
                }
                Ok(())
            }
            RunnerError::FileProcessingError { source } => write!(f, "{source:#}"),
        }
    }
}

impl std::error::Error for RunnerError {}

#[derive(Default)]
struct SageResults {
    ms1: Vec<ProcessedSpectrum>,
    features: Vec<Feature>,
    quant: Vec<TmtQuant>,
}

impl FromParallelIterator<SageResults> for SageResults {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = SageResults>,
    {
        par_iter
            .into_par_iter()
            .reduce(SageResults::default, |mut acc, x| {
                acc.features.extend(x.features);
                acc.quant.extend(x.quant);
                acc.ms1.extend(x.ms1);
                acc
            })
    }
}

impl FromIterator<SageResults> for SageResults {
    fn from_iter<I>(par_iter: I) -> Self
    where
        I: IntoIterator<Item = SageResults>,
    {
        par_iter
            .into_iter()
            .fold(SageResults::default(), |mut acc, x| {
                acc.features.extend(x.features);
                acc.quant.extend(x.quant);
                acc.ms1.extend(x.ms1);
                acc
            })
    }
}

impl Runner {
    pub fn new(parameters: Search) -> anyhow::Result<Self> {
        let start = Instant::now();
        let database = parameters.database.clone().build().with_context(|| {
            format!(
                "Failed to build database from `{}`",
                parameters.database.fasta
            )
        })?;
        info!(
            "generated {} fragments in {}ms",
            database.size(),
            (Instant::now() - start).as_millis()
        );
        Ok(Self {
            database,
            parameters,
            start,
        })
    }

    fn spectrum_fdr(&self, features: &mut [Feature]) -> usize {
        if sage_core::ml::linear_discriminant::score_psms(features).is_none() {
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
        spectra: Vec<ProcessedSpectrum>,
    ) -> SageResults {
        let features: Vec<_> = spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= self.parameters.min_peaks && spec.level == 2)
            .flat_map(|spec| scorer.score(spec, self.parameters.report_psms))
            .collect();

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
                sage_core::tmt::quantify(&spectra, isobaric, Tolerance::Ppm(-20.0, 20.0), level)
            })
            .unwrap_or_default();
        let ms1 = spectra.into_iter().filter(|s| s.level == 1).collect();

        SageResults {
            features,
            quant,
            ms1,
        }
    }

    fn process_chunk(
        &self,
        scorer: &Scorer,
        chunk: &[String],
        chunk_idx: usize,
        batch_size: usize,
    ) -> anyhow::Result<SageResults> {
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

        let (valids, errors): (Vec<_>, Vec<_>) = chunk
            .par_iter()
            .map(|path| {
                sage_core::read_mzml(path, sn)
                    .with_context(|| format!("Failed to process `{path}`"))
            })
            .partition(Result::is_ok);

        if !errors.is_empty() {
            let errors: Vec<_> = errors
                .into_iter()
                .map(Result::unwrap_err)
                .map(|err| RunnerError::FileProcessingError { source: err })
                .collect();

            anyhow::bail!(RunnerError::ChunkProcessingError {
                chunk_idx,
                source: errors
            });
        }

        let spectra: Vec<_> = valids.into_iter().map(Result::unwrap).collect();

        let io_time = Instant::now() - start;

        let count: usize = spectra.iter().map(|x| x.as_slice().len()).sum();

        let start = Instant::now();
        let results = spectra
            .into_par_iter()
            .enumerate()
            .map(|(idx, spectra)| {
                let sp = SpectrumProcessor::new(
                    self.parameters.max_peaks,
                    self.parameters.database.fragment_min_mz,
                    self.parameters.database.fragment_max_mz,
                    self.parameters.deisotope,
                    chunk_idx * batch_size + idx,
                );

                let spectra = spectra
                    .into_iter()
                    .map(move |spec| sp.process(spec))
                    .collect::<Vec<_>>();
                self.search_processed_spectra(scorer, spectra)
            })
            .collect::<SageResults>();
        let search_time = Instant::now() - start;
        info!(" - file IO: {:8} ms", io_time.as_millis());
        info!(
            " - search:  {:8} ms ({} spectra)",
            search_time.as_millis(),
            count
        );
        Ok(results)
    }

    pub fn batch_files(&self, scorer: &Scorer, batch_size: usize) -> anyhow::Result<SageResults> {
        let outputs: anyhow::Result<SageResults> = self
            .parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| self.process_chunk(scorer, chunk, chunk_idx, batch_size))
            .collect();
        outputs
    }

    pub fn run(mut self, parallel: usize) -> anyhow::Result<()> {
        let scorer = Scorer {
            db: &self.database,
            precursor_tol: self.parameters.precursor_tol,
            fragment_tol: self.parameters.fragment_tol,
            min_matched_peaks: self.parameters.min_matched_peaks,
            min_isotope_err: self.parameters.isotope_errors.0,
            max_isotope_err: self.parameters.isotope_errors.1,
            max_fragment_charge: self.parameters.max_fragment_charge,
            min_fragment_mass: self.parameters.database.fragment_min_mz,
            max_fragment_mass: self.parameters.database.fragment_max_mz,
            chimera: self.parameters.chimera,
            wide_window: self.parameters.wide_window,
        };

        // Collect all results into a single container
        let mut outputs = self
            .batch_files(&scorer, parallel)
            .context("Failed to process mzML files")?;

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

        if let Some(alignments) = alignments {
            if self.parameters.quant.lfq {
                let mut areas = sage_core::lfq::build_feature_map(
                    self.parameters.quant.lfq_settings,
                    &outputs.features,
                )
                .quantify(&self.database, &outputs.ms1, &alignments);

                let q_precursor = sage_core::fdr::picked_precursor(&mut areas);

                log::info!("discovered {} target MS1 peaks at 5% FDR", q_precursor);
                self.parameters
                    .output_paths
                    .push(self.write_lfq(areas, &filenames)?);
            }
        }

        log::info!(
            "discovered {} target peptide-spectrum matches at 1% FDR",
            q_spectrum
        );
        log::info!("discovered {} target peptides at 1% FDR", q_peptide);
        log::info!("discovered {} target proteins at 1% FDR", q_protein);
        log::trace!("writing outputs");

        if self.parameters.write_pin {
            self.parameters
                .output_paths
                .push(self.write_pin(&outputs.features, &filenames)?);
        }

        self.parameters
            .output_paths
            .push(self.write_features(outputs.features, &filenames)?);

        if !outputs.quant.is_empty() {
            self.parameters
                .output_paths
                .push(self.write_tmt(outputs.quant, &filenames)?);
        }

        let path = self.make_path("results.json");
        self.parameters.output_paths.push(path.to_string());
        println!("{}", serde_json::to_string_pretty(&self.parameters)?);

        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        path.write_bytes_sync(bytes)?;

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);

        Ok(())
    }
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("sage")
        .version(clap::crate_version!())
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .about("\u{1F52E} Sage \u{1F9D9} - Proteomics searching so fast it feels like magic!")
        .arg(
            Arg::new("parameters")
                .required(true)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Path to configuration parameters (JSON file)")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("mzml_paths")
                .num_args(1..)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Paths to mzML files to process. Overrides mzML files listed in the \
                     configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to FASTA database. Overrides the FASTA file \
                     specified in the configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("output_directory")
                .short('o')
                .long("output_directory")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path where search and quant results will be written. \
                     Overrides the directory specified in the configuration file.",
                )
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("batch-size")
                .long("batch-size")
                .value_parser(value_parser!(u16).range(1..)) // value_parse range not implemented for usize
                .help("Number of files to load and search in parallel (default = # of CPUs/2)")
                .value_hint(ValueHint::Other),
        )
        .arg(
            Arg::new("write-pin")
                .long("write-pin")
                .action(clap::ArgAction::SetTrue)
                .help("Write percolator-compatible `.pin` output files"),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    // get batch-size value first, or default value, as u16 then casted back to usize
    let parallel = matches
        .get_one::<u16>("batch-size")
        .copied()
        .unwrap_or_else(|| num_cpus::get() as u16 / 2) as usize;

    let input = Input::from_arguments(matches)?;

    let runner = input.build().and_then(Runner::new)?;

    runner.run(parallel)?;

    Ok(())
}
