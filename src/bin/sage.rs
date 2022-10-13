use clap::{Arg, Command};
use log::{info, warn};
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
    quant: Option<Isobaric>,
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
    chimera: Option<bool>,
    min_peaks: Option<usize>,
    max_peaks: Option<usize>,
    max_fragment_charge: Option<u8>,
    isotope_errors: Option<(i8, i8)>,
    deisotope: Option<bool>,
    quant: Option<Isobaric>,
    predict_rt: Option<bool>,
    output_directory: Option<PathBuf>,
    mzml_paths: Option<Vec<PathBuf>>,
}

impl Search {
    pub fn load<P: AsRef<Path>>(
        path: P,
        mzml_paths: Option<Vec<P>>,
        fasta: Option<P>,
        output_directory: Option<P>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let mut file = std::fs::File::open(path)?;
        let mut request: Input = serde_json::from_reader(&mut file)?;
        if let Some(f) = fasta {
            request.database.update_fasta(f.as_ref().to_path_buf())
        };
        let database = request.database.make_parameters();
        let isotope_errors = request.isotope_errors.unwrap_or((0, 0));
        if isotope_errors.0 > isotope_errors.1 {
            log::warn!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
        }
        let mzml_paths = match mzml_paths {
            Some(p) => p.into_iter().map(|f| f.as_ref().to_path_buf()).collect(),
            _ => request.mzml_paths.expect("'mzml_paths' must be provided!"),
        };

        let output_directory = match output_directory {
            Some(p) => Some(p.as_ref().to_path_buf()),
            _ => request.output_directory,
        };

        Ok(Search {
            database,
            mzml_paths,
            output_directory,
            quant: request.quant,
            precursor_tol: request.precursor_tol,
            fragment_tol: request.fragment_tol,
            report_psms: request.report_psms.unwrap_or(1),
            max_peaks: request.max_peaks.unwrap_or(150),
            min_peaks: request.min_peaks.unwrap_or(15),
            max_fragment_charge: request.max_fragment_charge,
            isotope_errors: request.isotope_errors.unwrap_or((0, 0)),
            deisotope: request.deisotope.unwrap_or(true),
            chimera: request.chimera.unwrap_or(false),
            predict_rt: request.predict_rt.unwrap_or(true),
            pin_paths: Vec::new(),
            search_time: 0.0,
        })
    }
}

fn process_mzml_file_sps<P: AsRef<Path>>(
    p: P,
    search: &Search,
    scorer: &Scorer,
) -> Result<PathBuf, Box<dyn std::error::Error + Send + Sync + 'static>> {
    let sp = SpectrumProcessor::new(
        search.max_peaks,
        search.database.fragment_min_mz,
        search.database.fragment_max_mz,
        search.deisotope,
    );

    // TODO:
    //  - more robust error checking
    //  - better error messages, there are several places we use try operator
    //      that could be replaced with a custom error type
    if p.as_ref()
        .extension()
        .expect("expecting .mzML files as input!")
        .to_ascii_lowercase()
        != "mzml"
    {
        panic!("expecting .mzML files as input!")
    }

    let spectra = sage::mzml::MzMlReader::read(&p)?
        .into_par_iter()
        .map(|spec| sp.process(spec))
        .collect::<Vec<_>>();
    log::trace!("{}: read {} spectra", p.as_ref().display(), spectra.len());

    let mut path = p.as_ref().to_path_buf();
    path.set_extension("quant.csv");

    if let Some(mut directory) = search.output_directory.clone() {
        // If directory doesn't exist, attempt to create it
        if !directory.exists() {
            std::fs::create_dir_all(&directory)?;
        }
        directory.push(path.file_name().expect("BUG: should be a filename!"));
        path = directory;
    }

    let mut scores = Vec::new();

    if let Some(quant) = search.quant.as_ref() {
        let mut wtr = csv::WriterBuilder::default().from_path(&path)?;
        let mut headers = vec![
            "scannr".into(),
            "ms3_injection".into(),
            "peptide".into(),
            "sps_purity".into(),
            "correct_precursors".into(),
        ];
        headers.extend(quant.headers());
        wtr.write_record(&headers)?;

        for spectrum in &spectra {
            if spectrum.level == 3 {
                if let Some(quant) = sage::tmt::quantify_sps(
                    scorer,
                    &spectra,
                    spectrum,
                    quant,
                    Tolerance::Ppm(-20.0, 20.0),
                ) {
                    let mut v = vec![
                        quant.hit.scannr.to_string(),
                        spectrum.ion_injection_time.to_string(),
                        quant.hit.peptide.clone(),
                        quant.hit_purity.ratio.to_string(),
                        quant.hit_purity.correct_precursors.to_string(),
                    ];
                    v.extend(
                        quant
                            .intensities
                            .iter()
                            .map(|peak| peak.map(|p| p.intensity.to_string()).unwrap_or_default()),
                    );
                    wtr.write_record(v)?;
                    scores.push(quant.hit);

                    if let Some(chimera) = quant.chimera {
                        scores.push(chimera);
                    }
                }
            }
        }
        wtr.flush()?;
    } else {
        scores = spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= search.min_peaks && spec.level == 2)
            .flat_map(|spec| scorer.score(spec, search.report_psms))
            .collect();
    }

    if sage::lda::score_psms(scorer.db, &mut scores, search.predict_rt).is_some() {
        (&mut scores)
            .par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));
    } else {
        log::warn!("linear model fitting failed, falling back to poisson-based FDR calculation");
        (&mut scores).par_sort_unstable_by(|a, b| a.poisson.total_cmp(&b.poisson));
    }

    let passing_psms = sage::lda::assign_q_values(&mut scores);

    let mut path = p.as_ref().to_path_buf();
    path.set_extension("sage.pin");

    if let Some(mut directory) = search.output_directory.clone() {
        // If directory doesn't exist, attempt to create it
        if !directory.exists() {
            std::fs::create_dir_all(&directory)?;
        }
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
    let output_directory = matches.get_one::<String>("output_directory");
    let fasta = matches.get_one::<String>("fasta");
    let mzml_paths = matches
        .get_many::<String>("mzml_paths")
        .map(|vals| vals.collect());

    let mut search = Search::load(path, mzml_paths, fasta, output_directory)?;

    let db = search.database.clone().build()?;

    info!(
        "generated {} fragments in {}ms",
        db.size(),
        (Instant::now() - start).as_millis()
    );

    if search.chimera && search.report_psms != 1 {
        warn!("chimeric search turned on, but report_psms is not 1 - overriding");
    }

    let scorer = Scorer::new(
        &db,
        search.precursor_tol,
        search.fragment_tol,
        search.isotope_errors.0,
        search.isotope_errors.1,
        search.max_fragment_charge,
        search.database.fragment_min_mz,
        search.database.fragment_max_mz,
        search.chimera,
    );

    let output_paths = search
        .mzml_paths
        .par_iter()
        .map(|ms2_path| process_mzml_file_sps(ms2_path, &search, &scorer))
        .collect::<Vec<_>>();

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
