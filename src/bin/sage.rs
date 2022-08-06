use clap::{Arg, Command};
use log::info;
use rayon::prelude::*;
use sage::mass::Tolerance;
use sage::scoring::{assign_q_values, Scorer};
use sage::spectrum::SpectrumProcessor;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::time::{self, Instant};

#[derive(Serialize)]
/// Actual search parameters - may include overrides or default values not set by user
struct Search {
    database: sage::database::Parameters,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    isotope_errors: (i8, i8),
    min_peaks: usize,
    max_peaks: usize,
    report_psms: usize,
    process_files_parallel: bool,
    mzml_paths: Vec<PathBuf>,
    pin_paths: Vec<PathBuf>,
    search_time: f32,

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
    min_peaks: Option<usize>,
    max_peaks: Option<usize>,
    isotope_errors: Option<(i8, i8)>,
    process_files_parallel: Option<bool>,
    output_directory: Option<PathBuf>,
    mzml_paths: Vec<PathBuf>,
}

impl Search {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let mut file = std::fs::File::open(path)?;
        let request: Input = serde_json::from_reader(&mut file)?;
        let database = request.database.make_parameters();
        let isotope_errors = request.isotope_errors.unwrap_or((0, 0));
        if isotope_errors.0 > isotope_errors.1 {
            panic!("Minimum isotope_error value greater than maximum! Correct usage: `isotope_errors: [-1, 3]`");
        }
        Ok(Search {
            database,
            precursor_tol: request.precursor_tol,
            fragment_tol: request.fragment_tol,
            report_psms: request.report_psms.unwrap_or(1),
            max_peaks: request.max_peaks.unwrap_or(150),
            min_peaks: request.min_peaks.unwrap_or(15),
            isotope_errors: request.isotope_errors.unwrap_or((0, 0)),
            pin_paths: Vec::new(),
            mzml_paths: request.mzml_paths,
            process_files_parallel: request.process_files_parallel.unwrap_or(true),
            output_directory: request.output_directory,
            search_time: 0.0,
        })
    }
}

fn process_mzml_file<P: AsRef<Path>>(
    p: P,
    search: &Search,
    scorer: &Scorer,
) -> Result<PathBuf, Box<dyn std::error::Error + Send + Sync + 'static>> {
    let sp = SpectrumProcessor::new(search.max_peaks, search.database.fragment_max_mz);

    if p.as_ref()
        .extension()
        .expect("expecting .mzML files as input!")
        .to_ascii_lowercase()
        != "mzml"
    {
        panic!("expecting .mzML files as input!")
    }

    let mut scores = sage::mzml::MzMlReader::read_ms2(&p)?
        .into_par_iter()
        .filter(|spec| spec.mz.len() >= search.min_peaks)
        .filter_map(|spec| sp.process(spec))
        .flat_map(|spec| scorer.score(&spec, search.report_psms))
        .collect::<Vec<_>>();

    (&mut scores).par_sort_unstable_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));
    let passing_psms = assign_q_values(&mut scores);

    let mut path = p.as_ref().to_path_buf();
    path.set_extension("sage.pin");

    if let Some(mut directory) = search.output_directory.clone() {
        directory.push(path.file_name().expect("BUG: should be a filename!"));
        path = directory;
    }

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&path)?;

    let total_psms = scores.len();

    for (idx, mut score) in scores.into_iter().enumerate() {
        score.specid = idx;
        writer.serialize(score)?;
    }

    info!(
        "{:?}: assigned {} PSMs ({} with 1% FDR)",
        p.as_ref(),
        total_psms,
        passing_psms,
    );
    Ok(path)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let env = env_logger::Env::default().filter_or("sage_LOG", "info");
    env_logger::init_from_env(env);

    let start = time::Instant::now();

    let matches = Command::new("sage")
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .arg(Arg::new("parameters").required(true))
        .get_matches();

    let path = matches
        .get_one::<String>("parameters")
        .expect("required parameters");
    let mut search = Search::load(path)?;

    let db = search.database.clone().build()?;

    info!(
        "generated {} fragments in {}ms",
        db.size(),
        (Instant::now() - start).as_millis()
    );

    let scorer = Scorer::new(
        &db,
        search.precursor_tol,
        search.fragment_tol,
        search.isotope_errors.0,
        search.isotope_errors.1,
    );

    let output_paths = match search.process_files_parallel {
        true => search
            .mzml_paths
            .par_iter()
            .map(|ms2_path| process_mzml_file(ms2_path, &search, &scorer))
            .collect::<Vec<_>>(),
        false => search
            .mzml_paths
            .iter()
            .map(|ms2_path| process_mzml_file(ms2_path, &search, &scorer))
            .collect::<Vec<_>>(),
    };

    search.search_time = (Instant::now() - start).as_secs_f32();

    let mut failures = 0;
    search.pin_paths = search
        .mzml_paths
        .iter()
        .zip(output_paths.into_iter())
        .filter_map(|(input, output)| match output {
            Ok(path) => Some(path),
            Err(err) => {
                eprintln!(
                    "Encountered error while processing {}: {}",
                    input.as_path().to_string_lossy(),
                    err
                );
                failures += 1;
                None
            }
        })
        .collect::<Vec<_>>();

    let results = serde_json::to_string_pretty(&search)?;

    println!("{}", &results);
    std::fs::write("results.json", results)?;

    Ok(())
}
