use clap::{Arg, Command};
use log::info;
use rayon::prelude::*;
use sage_cloudpath::CloudPath;
use sage_core::database::{Builder, IndexedDatabase, Parameters};
use sage_core::mass::Tolerance;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectrum::{ProcessedSpectrum, SpectrumProcessor};
use sage_core::tmt::{Isobaric, TmtQuant};
use serde::{Deserialize, Serialize};
use std::time::{self, Instant};

#[derive(Serialize)]
/// Actual search parameters - may include overrides or default values not set by user
struct Search {
    database: Parameters,
    quant: Quant,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    isotope_errors: (i8, i8),
    deisotope: bool,
    chimera: bool,
    min_peaks: usize,
    max_peaks: usize,
    max_fragment_charge: Option<u8>,
    report_psms: usize,
    predict_rt: bool,
    parallel: bool,
    mzml_paths: Vec<String>,
    output_paths: Vec<String>,

    #[serde(skip_serializing)]
    output_directory: CloudPath,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
struct Input {
    database: Builder,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    report_psms: Option<usize>,
    chimera: Option<bool>,
    min_peaks: Option<usize>,
    max_peaks: Option<usize>,
    max_fragment_charge: Option<u8>,
    isotope_errors: Option<(i8, i8)>,
    deisotope: Option<bool>,
    quant: Option<Quant>,
    predict_rt: Option<bool>,
    output_directory: Option<String>,
    parallel: Option<bool>,
    mzml_paths: Option<Vec<String>>,
}

#[derive(Serialize, Deserialize, Default)]
struct Quant {
    tmt: Option<Isobaric>,
    lfq: Option<bool>,
}

impl Input {
    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        let mut file = std::fs::File::open(path.as_ref())?;
        serde_json::from_reader(&mut file).map_err(anyhow::Error::from)
    }

    pub fn build(self) -> anyhow::Result<Search> {
        let database = self.database.make_parameters();
        let isotope_errors = self.isotope_errors.unwrap_or((0, 0));
        if isotope_errors.0 > isotope_errors.1 {
            log::warn!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
        }

        let mzml_paths = self.mzml_paths.expect("'mzml_paths' must be provided!");

        let output_directory = match self.output_directory {
            Some(path) => {
                let path = path.parse::<CloudPath>()?;
                if let CloudPath::Local(p) = &path {
                    std::fs::create_dir_all(p)?;
                }
                path
            }
            None => CloudPath::Local(std::env::current_dir()?),
        };

        Ok(Search {
            database,
            quant: self.quant.unwrap_or_default(),
            parallel: self.parallel.unwrap_or(true),
            mzml_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            report_psms: self.report_psms.unwrap_or(1),
            max_peaks: self.max_peaks.unwrap_or(150),
            min_peaks: self.min_peaks.unwrap_or(15),
            max_fragment_charge: self.max_fragment_charge,
            isotope_errors: self.isotope_errors.unwrap_or((0, 0)),
            deisotope: self.deisotope.unwrap_or(true),
            chimera: self.chimera.unwrap_or(false),
            predict_rt: self.predict_rt.unwrap_or(true),
            output_paths: Vec::new(),
        })
    }
}

struct Runner {
    database: IndexedDatabase,
    parameters: Search,
    start: Instant,
}

#[derive(Default)]
struct SageResults {
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
                acc
            })
    }
}

impl Runner {
    pub fn new(start: Instant, parameters: Search) -> anyhow::Result<Self> {
        let database = parameters.database.clone().build()?;
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
        if sage_core::ml::linear_discriminant::score_psms(features).is_some() {
            features
                .par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));
        } else {
            log::warn!(
                "linear model fitting failed, falling back to poisson-based FDR calculation"
            );
            features.par_sort_unstable_by(|a, b| a.poisson.total_cmp(&b.poisson));
        }
        sage_core::ml::qvalue::spectrum_q_value(features)
    }

    // Create a path for `file_name` in the specified output directory, if it exists,
    // otherwise, write to current directory
    fn make_path<S: AsRef<str>>(&self, file_name: S) -> CloudPath {
        let mut path = self.parameters.output_directory.clone();
        path.push(file_name);
        path
    }

    fn serialize_feature(&self, feature: &Feature, filenames: &[String]) -> csv::ByteRecord {
        let mut record = csv::ByteRecord::new();
        record.push_field(feature.peptide.as_str().as_bytes());
        record.push_field(feature.proteins.as_str().as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.num_proteins).as_bytes());
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
        record.push_field(ryu::Buffer::new().format(feature.isotope_error).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_mass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.hyperscore).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_hyperscore)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_rt).as_bytes());
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
        record.push_field(ryu::Buffer::new().format(feature.q_value).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.peptide_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.protein_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ms1_intensity).as_bytes());
        record
    }

    fn write_features(
        &self,
        features: Vec<Feature>,
        filenames: &[String],
    ) -> anyhow::Result<String> {
        let path = self.make_path("results.sage.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        let headers = csv::ByteRecord::from(vec![
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
            "isotope_error",
            "precursor_ppm",
            "fragment_ppm",
            "hyperscore",
            "delta_hyperscore",
            "rt",
            "predicted_rt",
            "delta_rt",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "matched_intensity_pct",
            "scored_candidates",
            "poisson",
            "sage_discriminant_score",
            "posterior_error",
            "spectrum_fdr",
            "peptide_fdr",
            "protein_fdr",
            "ms1_intensity",
        ]);

        wtr.write_byte_record(&headers)?;
        for record in features
            .into_par_iter()
            .map(|feat| self.serialize_feature(&feat, filenames))
            .collect::<Vec<_>>()
        {
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    fn write_quant(&self, quant: &[TmtQuant], filenames: &[String]) -> anyhow::Result<String> {
        let path = self.make_path("quant.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = csv::ByteRecord::from(vec!["file", "scannr", "ion_injection_time"]);
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

    fn search_processed_spectra(
        &self,
        scorer: &Scorer,
        spectra: Vec<ProcessedSpectrum>,
    ) -> SageResults {
        let mut features: Vec<_> = spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= self.parameters.min_peaks && spec.level == 2)
            .flat_map(|spec| scorer.score(spec, self.parameters.report_psms))
            .collect();

        if self.parameters.predict_rt {
            let _ = sage_core::ml::retention_model::predict(&self.database, &mut features);
        }

        if self.parameters.quant.lfq.unwrap_or(false) {
            sage_core::lfq::quantify(&mut features, &spectra);
        }

        let quant = self
            .parameters
            .quant
            .tmt
            .as_ref()
            .map(|isobaric| {
                sage_core::tmt::quantify(&spectra, isobaric, Tolerance::Ppm(-20.0, 20.0), 3)
            })
            .unwrap_or_default();
        SageResults { features, quant }
    }

    fn process_file<P: AsRef<str>>(
        &self,
        scorer: &Scorer,
        path: P,
        file_id: usize,
    ) -> anyhow::Result<SageResults> {
        let sp = SpectrumProcessor::new(
            self.parameters.max_peaks,
            self.parameters.database.fragment_min_mz,
            self.parameters.database.fragment_max_mz,
            self.parameters.deisotope,
            file_id,
        );

        let spectra = sage_cloudpath::read_mzml(&path)?
            .into_par_iter()
            .map(|spec| sp.process(spec))
            .collect::<Vec<_>>();

        log::trace!("{}: read {} spectra", path.as_ref(), spectra.len());

        Ok(self.search_processed_spectra(scorer, spectra))
    }

    fn process_chunk(
        &self,
        scorer: &Scorer,
        chunk: &[String],
        chunk_idx: usize,
        batch_size: usize,
    ) -> SageResults {
        // Read all of the spectra at once - this can help prevent memory over-consumption issues
        let spectra = chunk
            .par_iter()
            .map(sage_cloudpath::read_mzml)
            .collect::<Vec<_>>();

        spectra
            .into_par_iter()
            .enumerate()
            .filter_map(|(idx, spectra)| match spectra {
                Ok(spectra) => {
                    log::info!("{}: read {} spectra", chunk[idx], spectra.len());
                    Some((idx, spectra))
                }
                Err(e) => {
                    log::error!("error while processing {}: {}", chunk[idx], e);
                    None
                }
            })
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
            .collect::<SageResults>()
    }

    pub fn batch_files(&self, scorer: &Scorer) -> SageResults {
        let batch_size = num_cpus::get() / 2;
        self.parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| self.process_chunk(scorer, chunk, chunk_idx, batch_size))
            .collect::<SageResults>()
    }

    pub fn run(mut self) -> anyhow::Result<()> {
        let scorer = Scorer::new(
            &self.database,
            self.parameters.precursor_tol,
            self.parameters.fragment_tol,
            self.parameters.isotope_errors.0,
            self.parameters.isotope_errors.1,
            self.parameters.max_fragment_charge,
            self.parameters.database.fragment_min_mz,
            self.parameters.database.fragment_max_mz,
            self.parameters.chimera,
        );

        //Collect all results into a single container
        let mut outputs = match self.parameters.parallel {
            true => self.batch_files(&scorer),
            false => self
                .parameters
                .mzml_paths
                .iter()
                .enumerate()
                .flat_map(|(file_id, path)| self.process_file(&scorer, path, file_id))
                .collect::<SageResults>(),
        };

        let q_spectrum = self.spectrum_fdr(&mut outputs.features);
        let q_peptide = sage_core::fdr::picked_peptide(&self.database, &mut outputs.features);
        let q_protein = sage_core::fdr::picked_protein(&self.database, &mut outputs.features);

        info!(
            "discovered {} peptide-spectrum matches at 1% FDR",
            q_spectrum
        );
        info!("discovered {} peptides at 1% FDR", q_peptide);
        info!("discovered {} proteins at 1% FDR", q_protein);

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

        self.parameters
            .output_paths
            .push(self.write_features(outputs.features, &filenames)?);
        if !outputs.quant.is_empty() {
            self.parameters
                .output_paths
                .push(self.write_quant(&outputs.quant, &filenames)?);
        }

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);

        let path = self.make_path("results.json");
        self.parameters.output_paths.push(path.to_string());
        println!("{}", serde_json::to_string_pretty(&self.parameters)?);

        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        path.write_bytes_sync(bytes)?;

        Ok(())
    }
}

fn main() -> anyhow::Result<()> {
    // let env = env_logger::Env::default().filter_or("SAGE_LOG", "info");
    env_logger::builder()
        .filter(None, log::LevelFilter::Error)
        .filter_module("sage", log::LevelFilter::Trace)
        .init();

    let start = time::Instant::now();

    let matches = Command::new("sage")
        .version(clap::crate_version!())
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .about("\u{1F52E} Sage \u{1F9D9} - Proteomics searching so fast it feels like magic!")
        .arg(
            Arg::new("parameters")
                .required(true)
                .help("The search parameters as a JSON file."),
        )
        .arg(Arg::new("mzml_paths").num_args(1..).help(
            "mzML files to analyze. Overrides mzML files listed in the \
                     parameter file.",
        ))
        .arg(Arg::new("fasta").short('f').long("fasta").help(
            "The FASTA protein database. Overrides the FASTA file \
                     specified in the parameter file.",
        ))
        .arg(
            Arg::new("output_directory")
                .short('o')
                .long("output_directory")
                .help(
                    "Where the search and quant results will be written. \
                     Overrides the directory specified in the parameter file.",
                ),
        )
        .arg(
            Arg::new("no-parallel")
                .long("no-parallel")
                .action(clap::ArgAction::SetFalse)
                .help(
                    "Turn off parallel file searching. \
                 Useful for memory constrained systems.",
                ),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    let path = matches
        .get_one::<String>("parameters")
        .expect("required parameters");
    let mut input = Input::load(path)?;

    // Handle JSON configuration overrides
    if let Some(output_directory) = matches.get_one::<String>("output_directory") {
        input.output_directory = Some(output_directory.into());
    }
    if let Some(fasta) = matches.get_one::<String>("fasta") {
        input.database.fasta = Some(fasta.into());
    }
    if let Some(mzml_paths) = matches.get_many::<String>("mzml_paths") {
        input.mzml_paths = Some(mzml_paths.into_iter().map(|p| p.into()).collect());
    }

    if let Some(no_parallel) = matches.get_one::<bool>("no-parallel").copied() {
        if !no_parallel {
            input.parallel = Some(false);
        }
    }

    let runner = input
        .build()
        .and_then(|search| Runner::new(start, search))?;

    runner.run()?;

    Ok(())
}
