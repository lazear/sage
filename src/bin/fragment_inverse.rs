use carina::database::{IndexedDatabase, PeptideIx, Theoretical};
use carina::ion_series::Kind;
use carina::mass::{Residue, Tolerance, PROTON};
use carina::peptide::TargetDecoy;
use carina::spectrum::{read_spectrum, ProcessedSpectrum};
use clap::{App, Arg, ArgMatches, Command};
use indicatif::ParallelProgressIterator;
use log::info;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::{self, Instant};

#[derive(Copy, Clone)]
pub struct Score {
    peptide: PeptideIx,
    matched_b: u32,
    matched_y: u32,
    summed_b: f32,
    summed_y: f32,
    hyperlog: f32,
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
    deltascore: f32,
    matched_peaks: u32,
    summed_intensity: f32,
    total_candidates: usize,
    q_value: f32,
}

impl Score {
    fn hyperscore(&self, fact_table: &[f32]) -> f32 {
        let i = (self.summed_b + 1.0) * (self.summed_y + 1.0);
        let m = fact_table[(self.matched_b as usize).min(fact_table.len() - 1)]
            * fact_table[(self.matched_y as usize).min(fact_table.len() - 1)];
        let score = (i * m as f32).ln();
        if score.is_finite() {
            score
        } else {
            0.0
        }
    }

    pub fn new(peptide: &Theoretical) -> Self {
        Score {
            peptide: peptide.peptide_index,
            matched_b: 0,
            matched_y: 0,
            summed_b: 0.0,
            summed_y: 0.0,
            hyperlog: 0.0,
        }
    }
}

pub struct Scorer<'db> {
    db: &'db IndexedDatabase,
    report: usize,
    fragment_tol: Tolerance,
    precursor_tol: Tolerance,
    factorial: [f32; 64],
}

impl<'db> Scorer<'db> {
    pub fn new(
        db: &'db IndexedDatabase,
        report: usize,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
    ) -> Self {
        let mut factorial = [1.0f32; 64];
        for i in 1..64 {
            factorial[i] = factorial[i - 1] * i as f32;
        }

        debug_assert!(factorial[3] == 6.0);
        debug_assert!(factorial[4] == 24.0);

        Scorer {
            db,
            report,
            fragment_tol,
            precursor_tol,
            factorial,
        }
    }

    pub fn score<'s>(&self, query: &ProcessedSpectrum) -> Vec<Percolator<'db>> {
        let mut scores: HashMap<PeptideIx, Score> = HashMap::new();
        let candidates = self.db.query(query, self.precursor_tol, self.fragment_tol);

        for (idx, fragment_mz) in query.mz.iter().enumerate() {
            for frag in candidates.page_search(*fragment_mz) {
                let mut sc = scores
                    .entry(frag.peptide_index)
                    .or_insert_with(|| Score::new(frag));

                let intensity = query.int[idx];
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
            }
        }

        let mut scores = scores
            .into_values()
            .filter(|sc| sc.matched_b + sc.matched_y > 0)
            .map(|mut sc| {
                sc.hyperlog = sc.hyperscore(&self.factorial);
                sc
            })
            .collect::<Vec<_>>();

        // FDR Calculation:
        // * Sort by score, descending
        // * Estimate FDR
        // * Calculate q-value

        scores.sort_by(|b, a| a.hyperlog.total_cmp(&b.hyperlog));
        let mut q_values = vec![0.0; scores.len()];
        let mut decoy = 1;
        let mut target = 0;

        for (idx, score) in scores.iter().enumerate() {
            match self.db[score.peptide] {
                TargetDecoy::Target(_) => target += 1,
                TargetDecoy::Decoy(_) => decoy += 1,
            }
            q_values[idx] = decoy as f32 / target as f32;
        }

        // Reverse array, and calculate the cumulative minimum
        let mut q_min = 1.0f32;
        for idx in (0..q_values.len()).rev() {
            q_min = q_min.min(q_values[idx]);
            q_values[idx] = q_min;
        }

        let mut reporting = Vec::new();
        if scores.is_empty() {
            return reporting;
        }
        for idx in 0..self.report.min(scores.len()) {
            let better = scores[idx];
            let next = scores
                .get(idx + 1)
                .map(|score| score.hyperlog)
                .unwrap_or_default();

            let peptide = self.db[better.peptide].peptide();
            reporting.push(Percolator {
                peptide: peptide.to_string(),
                proteins: &peptide.protein,
                specid: 0,
                scannr: query.scan,
                label: self.db[better.peptide].label(),
                expmass: query.precursor_mz,
                calcmass: peptide.monoisotopic + PROTON,
                charge: query.charge,
                rt: query.rt,
                delta_mass: (query.precursor_mz - peptide.monoisotopic - PROTON),
                hyperscore: better.hyperlog,
                deltascore: better.hyperlog - next,
                matched_peaks: better.matched_b + better.matched_y,
                summed_intensity: better.summed_b + better.summed_y,
                total_candidates: scores.len(),
                q_value: q_values[idx],
            })
        }
        reporting
    }
}

#[derive(Serialize)]
struct SearchParameters {
    database: carina::database::Parameters,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    report_psms: usize,
    ms2_paths: Vec<String>,
}

impl SearchParameters {
    pub fn from_args<P: AsRef<Path>>(
        matches: ArgMatches,
        save_to: P,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let parameter_path = matches.get_one::<String>("parameters").expect("required");
        let data = std::fs::read_to_string(parameter_path)?;
        let builder: carina::database::Builder = serde_json::from_str(&data)?;
        let params = builder.make_parameters();

        let precursor_tol = Tolerance::Ppm(
            matches
                .get_one::<String>("precursor_ppm")
                .expect("required")
                .parse::<f32>()?,
        );
        let fragment_tol = Tolerance::Ppm(
            matches
                .get_one::<String>("fragment_ppm")
                .expect("required")
                .parse::<f32>()?,
        );
        let report_psms = matches
            .get_one::<String>("report")
            .expect("required")
            .parse::<usize>()?;
        let ms2_paths = matches
            .get_many::<String>("spectra")
            .expect("required")
            .cloned()
            .collect();

        let params = SearchParameters {
            database: params,
            precursor_tol,
            fragment_tol,
            report_psms,
            ms2_paths,
        };

        info!("logging search parameters to {:?}", save_to.as_ref());
        std::fs::write(save_to, serde_json::to_string_pretty(&params)?)?;
        Ok(params)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let env = env_logger::Env::default().filter_or("CARINA_LOG", "info");
    env_logger::init_from_env(env);

    let start = time::Instant::now();

    let matches = Command::new("carina")
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .arg(Arg::new("parameters").required(true))
        .arg(Arg::new("precursor_ppm").required(true))
        .arg(Arg::new("fragment_ppm").required(true))
        .arg(Arg::new("report").required(false).default_value("1"))
        .arg(Arg::new("spectra").multiple_values(true).required(true))
        .get_matches();

    let search = SearchParameters::from_args(matches, "search.json")?;

    let db = search.database.build()?;
    info!(
        "generated {} fragments in {}ms",
        db.size(),
        (Instant::now() - start).as_millis()
    );

    for ms2_path in &search.ms2_paths {
        let mut spectra = read_spectrum(ms2_path)?
            .into_iter()
            .map(ProcessedSpectrum::from)
            .collect::<Vec<_>>();

        info!("{:?}: read {} spectra", ms2_path, spectra.len());
        spectra.sort_by(|a, b| a.precursor_mz.total_cmp(&b.precursor_mz));

        let scorer = Scorer::new(&db, 1, Tolerance::Ppm(50.0), Tolerance::Ppm(10.0));
        let start = Instant::now();
        let scores: Vec<Percolator> = spectra
            .par_iter()
            .progress()
            .flat_map(|spectra| scorer.score(spectra))
            .collect();
        let duration = Instant::now() - start;

        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(format!("results.pin"))?;

        for (idx, mut score) in scores.into_iter().enumerate() {
            score.specid = idx;
            writer.serialize(score)?;
        }

        info!(
            "{:?}: processed {} scans in {} ms ({} scans/sec)",
            ms2_path,
            spectra.len(),
            duration.as_millis(),
            spectra.len() as f32 / duration.as_secs_f32()
        );
    }
    info!(
        "all searches completed in {:0.2} s",
        (Instant::now() - start).as_secs_f32()
    );

    Ok(())
}
