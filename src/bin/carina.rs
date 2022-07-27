use carina::database::{IndexedDatabase, PeptideIx, Theoretical};
use carina::ion_series::Kind;
use carina::mass::{Tolerance, NEUTRON, PROTON};
use carina::spectrum::{read_ms2, ProcessedSpectrum, SpectrumProcessor};
use clap::{Arg, Command};
use log::info;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::time::{self, Instant};

#[derive(Copy, Clone)]
pub struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    hyperscore: f32,
}

#[derive(Serialize)]
pub struct Percolator<'db> {
    peptide: String,
    proteins: &'db str,
    specid: usize,
    scannr: u32,
    label: i32,
    expmass: f32,
    calcmass: f32,
    charge: u8,
    rt: f32,
    delta_mass: f32,
    hyperscore: f32,
    delta_hyperscore: f32,
    matched_peaks: u32,
    matched_intensity_pct: f32,
    delta_matched: f32,
    scored_candidates: usize,
    spectrum_z_score: f32,
    q_value: f32,
}

impl Score {
    /// Calculate the X!Tandem hyperscore
    /// * `fact_table` is a precomputed vector of factorials
    fn hyperscore(&self, fact_table: &[f32]) -> f32 {
        let i = (self.summed_b + 1.0) * (self.summed_y + 1.0);
        let m = fact_table[(self.matched_b as usize).min(fact_table.len() - 2)]
            * fact_table[(self.matched_y as usize).min(fact_table.len() - 2)];

        let score = i.ln() + m.ln();
        if score.is_finite() {
            score
        } else {
            f32::MAX
        }
    }

    pub fn new(peptide: &Theoretical) -> Self {
        Score {
            peptide: peptide.peptide_index,
            matched_b: 0,
            matched_y: 0,
            summed_b: 0.0,
            summed_y: 0.0,
            hyperscore: 0.0,
        }
    }
}

pub struct Scorer<'db> {
    db: &'db IndexedDatabase,
    search: &'db Search,
    factorial: [f32; 32],
}

impl<'db> Scorer<'db> {
    pub fn new(db: &'db IndexedDatabase, search: &'db Search) -> Self {
        let mut factorial = [1.0f32; 32];
        for i in 1..32 {
            factorial[i] = factorial[i - 1] * i as f32;
        }

        debug_assert!(factorial[3] == 6.0);
        debug_assert!(factorial[4] == 24.0);

        Scorer {
            db,
            search,
            factorial,
        }
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score<'s>(&self, query: &ProcessedSpectrum) -> Vec<Percolator<'db>> {
        // Create a new `IndexedQuery`

        let candidates = self
            .db
            .query(query, self.search.precursor_tol, self.search.fragment_tol);

        // Allocate space for all potential candidates - many potential candidates
        // will not have fragments matched, so we use `Option<Score>`
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut score_vector: Vec<Option<Score>> = vec![None; potential];

        let mut total_intensity = 0.0;
        let mut matches = 0;
        for (fragment_mz, intensity) in query.peaks.iter() {
            total_intensity += intensity;
            for frag in candidates.page_search(*fragment_mz) {
                let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                let mut sc = score_vector[idx].take().unwrap_or_else(|| Score::new(frag));

                match frag.kind {
                    Kind::B => {
                        sc.matched_b += 1;
                        sc.summed_b += intensity
                    }
                    Kind::Y => {
                        sc.matched_y += 1;
                        sc.summed_y += intensity
                    }
                }

                score_vector[idx] = Some(sc);
                matches += 1;
            }
        }

        if matches == 0 {
            return Vec::new();
        }

        // Now that we have processed all candidates, calculate the hyperscore
        let mut scores = score_vector
            .into_iter()
            .filter_map(|s| match s {
                Some(mut sc) => {
                    sc.hyperscore = sc.hyperscore(&self.factorial);
                    Some(sc)
                }
                None => None,
            })
            .collect::<Vec<_>>();

        scores.sort_unstable_by(|b, a| a.hyperscore.total_cmp(&b.hyperscore));

        // Calculate median & std deviation of hyperscores
        let median = scores[scores.len() / 2].hyperscore;
        let mut sd = scores
            .iter()
            .fold(0.0f32, |sum, x| sum + (x.hyperscore - median).powi(2));
        sd = (sd / scores.len().saturating_sub(1) as f32).sqrt();

        let mut reporting = Vec::new();

        for idx in 0..self.search.report_psms.min(scores.len()) {
            let better = scores[idx];
            let next = scores
                .get(idx + 1)
                .map(|score| score.hyperscore)
                .unwrap_or_default();

            let z_score = match sd.is_finite() && sd > 0.0 {
                true => (better.hyperscore - median) / sd,
                false => 0.1,
            };

            let peptide = self.db[better.peptide].peptide();
            reporting.push(Percolator {
                peptide: peptide.to_string(),
                proteins: &peptide.protein,
                specid: 0,
                scannr: query.scan,
                label: self.db[better.peptide].label(),
                expmass: query.monoisotopic_mass + PROTON,
                calcmass: peptide.monoisotopic + PROTON,
                charge: query.charge,
                rt: query.rt,
                delta_mass: (query.monoisotopic_mass - peptide.monoisotopic).abs(),
                hyperscore: better.hyperscore,
                delta_hyperscore: better.hyperscore - next,
                matched_peaks: (better.matched_b + better.matched_y) as u32,
                matched_intensity_pct: (better.summed_b + better.summed_y) / total_intensity,
                delta_matched: (better.matched_b + better.matched_y) as f32
                    / (matches as f32 / scores.len() as f32),
                scored_candidates: scores.len(),
                spectrum_z_score: z_score,
                q_value: 1.0,
            })
        }
        reporting
    }

    /// Assign q_values in place to a set of PSMs
    ///
    /// # Invariants
    /// * `scores` must be sorted in descending order (e.g. best PSM is first)
    pub fn assign_q_values(&self, scores: &mut [Percolator]) -> usize {
        // FDR Calculation:
        // * Sort by score, descending
        // * Estimate FDR
        // * Calculate q-value

        let mut decoy = 1;
        let mut target = 0;

        for score in scores.iter_mut() {
            match score.label == -1 {
                true => decoy += 1,
                false => target += 1,
            }
            score.q_value = decoy as f32 / target as f32;
        }

        // Reverse slice, and calculate the cumulative minimum
        let mut q_min = 1.0f32;
        let mut passing = 0;
        for score in scores.iter_mut().rev() {
            q_min = q_min.min(score.q_value);
            score.q_value = q_min;
            if q_min <= 0.01 {
                passing += 1;
            }
        }
        passing
    }
}

#[derive(Serialize)]
pub struct Search {
    database: carina::database::Parameters,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    max_fragment_charge: u8,
    min_peaks: usize,
    max_peaks: usize,
    report_psms: usize,
    process_files_parallel: bool,
    ms2_paths: Vec<PathBuf>,
    pin_paths: Vec<PathBuf>,
    search_time: f32,

    #[serde(skip_serializing)]
    output_directory: Option<PathBuf>,
}

#[derive(Deserialize)]
struct Input {
    database: carina::database::Builder,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    report_psms: Option<usize>,
    min_peaks: Option<usize>,
    max_peaks: Option<usize>,
    max_fragment_charge: Option<u8>,
    process_files_parallel: Option<bool>,
    output_directory: Option<PathBuf>,
    ms2_paths: Vec<PathBuf>,
}

impl Search {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let mut file = std::fs::File::open(path)?;
        let request: Input = serde_json::from_reader(&mut file)?;
        let database = request.database.make_parameters();
        Ok(Search {
            database,
            precursor_tol: request.precursor_tol,
            fragment_tol: request.fragment_tol,
            report_psms: request.report_psms.unwrap_or(1),
            max_peaks: request.max_peaks.unwrap_or(150),
            min_peaks: request.min_peaks.unwrap_or(15),
            max_fragment_charge: request.max_fragment_charge.unwrap_or(3),
            pin_paths: Vec::new(),
            ms2_paths: request.ms2_paths,
            process_files_parallel: request.process_files_parallel.unwrap_or(true),
            output_directory: request.output_directory,
            search_time: 0.0,
        })
    }
}

fn process_ms2_file<P: AsRef<Path>>(
    p: P,
    scorer: &Scorer,
) -> Result<PathBuf, Box<dyn std::error::Error + Send + Sync + 'static>> {
    let sp = SpectrumProcessor::new(
        scorer.search.max_peaks,
        scorer.search.max_fragment_charge,
        scorer.search.database.fragment_max_mz,
    );

    let mut scores = read_ms2(&p)?
        .into_par_iter()
        .filter(|spec| spec.peaks.len() >= scorer.search.min_peaks)
        .flat_map(|spec| scorer.score(&sp.process(spec)))
        .collect::<Vec<_>>();
    (&mut scores).par_sort_unstable_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));
    let passing_psms = scorer.assign_q_values(&mut scores);

    let mut path = p.as_ref().to_path_buf();
    path.set_extension("carina.pin");

    if let Some(mut directory) = scorer.search.output_directory.clone() {
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
    let env = env_logger::Env::default().filter_or("CARINA_LOG", "info");
    env_logger::init_from_env(env);

    let start = time::Instant::now();

    let matches = Command::new("carina")
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .arg(Arg::new("parameters").required(true))
        .get_matches();

    let path = matches
        .get_one::<String>("parameters")
        .expect("required parameters");
    let mut search = Search::load(path)?;

    let db = search.database.clone().build()?;

    let buckets = db.buckets();
    let mut avg_delta = 0.0;

    for i in 1..buckets.len() {
        let delta = buckets[i] - buckets[i - 1];
        avg_delta += delta;
    }
    dbg!(avg_delta / buckets.len() as f32);

    info!(
        "generated {} fragments in {}ms",
        db.size(),
        (Instant::now() - start).as_millis()
    );

    let scorer = Scorer::new(&db, &search);
    let start = Instant::now();

    let output_paths = match search.process_files_parallel {
        true => search
            .ms2_paths
            .par_iter()
            .map(|ms2_path| process_ms2_file(ms2_path, &scorer))
            .collect::<Vec<_>>(),
        false => search
            .ms2_paths
            .iter()
            .map(|ms2_path| process_ms2_file(ms2_path, &scorer))
            .collect::<Vec<_>>(),
    };

    search.search_time = (Instant::now() - start).as_secs_f32();

    let mut failures = 0;
    search.pin_paths = search
        .ms2_paths
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
