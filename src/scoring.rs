use crate::database::{IndexedDatabase, PeptideIx, Theoretical};
use crate::ion_series::Kind;
use crate::mass::{Tolerance, PROTON};
use crate::spectrum::ProcessedSpectrum;
use serde::Serialize;

/// Structure to hold temporary scores
#[derive(Copy, Clone)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    hyperscore: f32,
}

#[derive(Serialize)]
/// Features of a candidate peptide spectrum match
pub struct Percolator<'db> {
    /// Peptide sequence, including modifications e.g.: NC(+57.021)HK
    pub peptide: String,
    /// Name of *a* protein containing this peptide sequence
    pub proteins: &'db str,
    /// Arbitrary spectrum id
    pub specid: usize,
    /// MS2 scan number
    pub scannr: u32,
    /// Target/Decoy label, -1 is decoy, 1 is target
    pub label: i32,
    /// Experimental mass MH+
    pub expmass: f32,
    /// Calculated mass, MH+
    pub calcmass: f32,
    /// Reported precursor charge
    pub charge: u8,
    /// Retention time
    pub rt: f32,
    /// Difference between expmass and calcmass
    pub delta_mass: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f32,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_hyperscore: f32,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Difference between # of matched peaks for this candidate,
    /// and the # of matched peaks for the average candidate
    pub delta_matched: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: usize,
    /// Estimated spectrum z score using hyperscore
    pub spectrum_z_score: f32,
    /// Assigned q_value
    pub q_value: f32,
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
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    min_isotope_err: i8,
    max_isotope_err: i8,
    factorial: [f32; 32],
}

impl<'db> Scorer<'db> {
    pub fn new(
        db: &'db IndexedDatabase,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
        min_isotope_err: i8,
        max_isotope_err: i8,
    ) -> Self {
        let mut factorial = [1.0f32; 32];
        for i in 1..32 {
            factorial[i] = factorial[i - 1] * i as f32;
        }

        debug_assert!(factorial[3] == 6.0);
        debug_assert!(factorial[4] == 24.0);

        Scorer {
            db,
            precursor_tol,
            fragment_tol,
            min_isotope_err,
            max_isotope_err,
            factorial,
        }
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score(&self, query: &ProcessedSpectrum, report_psms: usize) -> Vec<Percolator<'db>> {
        // Create a new `IndexedQuery`

        let candidates = self.db.query(
            query,
            self.precursor_tol,
            self.fragment_tol,
            self.max_isotope_err,
            self.max_isotope_err,
        );

        // Allocate space for all potential candidates - many potential candidates
        // will not have fragments matched, so we use `Option<Score>`
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut score_vector: Vec<Option<Score>> = vec![None; potential];

        let mut total_intensity = 0.0;
        let mut matches = 0;
        for peak in query.peaks.iter() {
            total_intensity += peak.intensity;
            for frag in candidates.page_search(peak.mass) {
                let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                let mut sc = score_vector[idx].take().unwrap_or_else(|| Score::new(frag));

                match frag.kind {
                    Kind::B => {
                        sc.matched_b += 1;
                        sc.summed_b += peak.intensity
                    }
                    Kind::Y => {
                        sc.matched_y += 1;
                        sc.summed_y += peak.intensity
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

        for idx in 0..report_psms.min(scores.len()) {
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
}

/// Assign q_values in place to a set of PSMs, returning the number of PSMs
/// q <= 0.01
///
/// # Invariants
/// * `scores` must be sorted in descending order (e.g. best PSM is first)
pub fn assign_q_values(scores: &mut [Percolator]) -> usize {
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
