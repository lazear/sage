use clap::{Arg, Command};
use log::info;
use rayon::prelude::*;
use sage::mass::Tolerance;
use sage::scoring::Scorer;
use sage::spectrum::SpectrumProcessor;
use sage::tmt::Isobaric;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::time::{self, Instant};

#[derive(Serialize)]
/// Actual search parameters - may include overrides or default values not set by user
struct Search {
    database: sage::database::Parameters,
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
    mzml_paths: Vec<PathBuf>,
    pin_paths: Vec<PathBuf>,

    #[serde(skip_serializing)]
    output_directory: Option<PathBuf>,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
struct Input {
    database: sage::database::Builder,
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
    output_directory: Option<PathBuf>,
    parallel: Option<bool>,
    mzml_paths: Option<Vec<PathBuf>>,
}

#[derive(Serialize, Deserialize, Default)]
struct Quant {
    tmt: Option<Isobaric>,
    lfq: Option<bool>,
}

impl Input {
    pub fn load<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let mut file = std::fs::File::open(path)?;
        serde_json::from_reader(&mut file).map_err(anyhow::Error::from)
    }

    pub fn build(self) -> anyhow::Result<Search> {
        let database = self.database.make_parameters();
        let isotope_errors = self.isotope_errors.unwrap_or((0, 0));
        if isotope_errors.0 > isotope_errors.1 {
            log::warn!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
        }

        let mzml_paths = self.mzml_paths.expect("'mzml_paths' must be provided!");

        if let Some(dir) = &self.output_directory {
            std::fs::create_dir_all(&dir)?;
        }

        Ok(Search {
            database,
            quant: self.quant.unwrap_or_default(),
            parallel: self.parallel.unwrap_or(true),
            mzml_paths,
            output_directory: self.output_directory,
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
            pin_paths: Vec::new(),
        })
    }
}

struct Runner {
    database: sage::database::IndexedDatabase,
    parameters: Search,
    start: Instant,
}

#[derive(Default)]
struct SageResults {
    features: Vec<sage::scoring::Percolator>,
    quant: Vec<sage::tmt::TmtQuant>,
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

    fn spectrum_fdr(&self, features: &mut [sage::scoring::Percolator]) -> usize {
        if sage::ml::linear_discriminant::score_psms(features).is_some() {
            features
                .par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));
        } else {
            log::warn!(
                "linear model fitting failed, falling back to poisson-based FDR calculation"
            );
            features.par_sort_unstable_by(|a, b| a.poisson.total_cmp(&b.poisson));
        }
        sage::ml::qvalue::spectrum_q_value(features)
    }

    // Create a path for `file_name` in the specified output directory, if it exists,
    // otherwise, write to current directory
    fn make_path<S: AsRef<str>>(&self, file_name: S) -> Option<PathBuf> {
        self.parameters
            .output_directory
            .clone()
            .or_else(|| std::env::current_dir().ok())
            .map(|mut directory| {
                directory.push(file_name.as_ref());
                directory
            })
    }

    fn write_features(&self, features: &[sage::scoring::Percolator]) -> anyhow::Result<PathBuf> {
        let path = self
            .make_path("search.pin")
            .expect("no output directory specified, and current directory is invalid!");
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)?;
        for feat in features {
            wtr.serialize(feat)?;
        }
        wtr.flush()?;
        Ok(path)
    }

    fn write_quant(&self, quant: &[sage::tmt::TmtQuant]) -> anyhow::Result<PathBuf> {
        let path = self
            .make_path("quant.csv")
            .expect("no output directory specified, and current directory is invalid!");

        let mut wtr = csv::WriterBuilder::new().from_path(&path)?;
        let mut headers = vec![
            "file_id".into(),
            "scannr".into(),
            "ion_injection_time".into(),
        ];
        headers.extend(
            self.parameters
                .quant
                .tmt
                .as_ref()
                .map(|tmt| tmt.headers())
                .expect("TMT quant cannot be performed without setting this parameter"),
        );

        wtr.write_record(&headers)?;

        for q in quant {
            let mut record = vec![
                q.file_id.to_string(),
                q.scannr.to_string(),
                q.ion_injection_time.to_string(),
            ];
            record.extend(q.peaks.iter().map(|x| x.to_string()));
            wtr.write_record(&record)?;
        }
        wtr.flush()?;
        Ok(path)
    }

    fn process_file<P: AsRef<Path>>(
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

        let spectra = sage::mzml::MzMlReader::read(&path)?
            .into_par_iter()
            .map(|spec| sp.process(spec))
            .collect::<Vec<_>>();

        log::trace!(
            "{}: read {} spectra",
            path.as_ref().display(),
            spectra.len()
        );

        let mut features: Vec<_> = spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= self.parameters.min_peaks && spec.level == 2)
            .flat_map(|spec| scorer.score(spec, self.parameters.report_psms))
            .collect();

        if self.parameters.predict_rt {
            let _ = sage::ml::retention_model::predict(&self.database, &mut features);
        }

        if self.parameters.quant.lfq.unwrap_or(false) {
            sage::lfq::quantify(&mut features, &spectra);
        }

        let quant = self
            .parameters
            .quant
            .tmt
            .as_ref()
            .map(|isobaric| sage::tmt::quantify(&spectra, isobaric, Tolerance::Ppm(-20.0, 20.0), 3))
            .unwrap_or_default();

        Ok(SageResults { features, quant })
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
            true => self
                .parameters
                .mzml_paths
                .par_iter()
                .enumerate()
                .flat_map(|(file_id, path)| self.process_file(&scorer, path, file_id))
                .collect::<SageResults>(),
            false => self
                .parameters
                .mzml_paths
                .iter()
                .enumerate()
                .flat_map(|(file_id, path)| self.process_file(&scorer, path, file_id))
                .collect::<SageResults>(),
        };

        let q_spectrum = self.spectrum_fdr(&mut outputs.features);
        let q_peptide = sage::fdr::picked_peptide(&self.database, &mut outputs.features);
        let q_protein = sage::fdr::picked_protein(&self.database, &mut outputs.features);

        info!(
            "discovered {} peptide-spectrum matches at 1% FDR",
            q_spectrum
        );
        info!("discovered {} peptides at 1% FDR", q_peptide);
        info!("discovered {} proteins at 1% FDR", q_protein);

        self.parameters
            .pin_paths
            .push(self.write_features(&outputs.features)?);
        if !outputs.quant.is_empty() {
            self.parameters
                .pin_paths
                .push(self.write_quant(&outputs.quant)?);
        }

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);

        let path = self
            .make_path("results.json")
            .expect("no output directory specified, and current directory is invalid!");
        let writer = std::fs::File::create(&path)?;
        self.parameters.pin_paths.push(path);
        println!("{}", serde_json::to_string_pretty(&self.parameters)?);
        serde_json::to_writer_pretty(writer, &self.parameters)?;

        Ok(())
    }
}

fn main() -> anyhow::Result<()> {
    let env = env_logger::Env::default().filter_or("SAGE_LOG", "info");
    env_logger::init_from_env(env);

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
