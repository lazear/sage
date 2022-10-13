use crate::database::{binary_search_slice, IndexedDatabase, PeptideIx, Theoretical};
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON};
use crate::peptide::Peptide;
use crate::spectrum::{Peak, Precursor, ProcessedSpectrum};
use serde::Serialize;

/// Structure to hold temporary scores
#[derive(Copy, Clone, Default, Debug)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    hyperscore: f64,
    ppm_difference: f32,
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug)]
struct PreScore {
    peptide: PeptideIx,
    matched: u32,
}

#[derive(Serialize, Clone, Debug)]
/// Features of a candidate peptide spectrum match
pub struct Percolator {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    /// Peptide sequence, including modifications e.g.: NC(+57.021)HK
    pub peptide: String,
    /// Internal peptide decoy index, used to match forward & reversed sequences
    /// See: Lin et al., DOI: 10.1021/acs.jproteome.2c00282
    pub peptide_decoy_idx: u16,
    /// Peptide length
    pub peptide_len: usize,
    /// Proteins containing this peptide sequence
    pub proteins: String,
    /// Number of proteins assigned to this peptide sequence
    pub num_proteins: usize,
    /// Spectrum id - empty until search is done
    pub specid: String,
    /// File identifier
    pub file_id: usize,
    /// MS2 scan number
    pub scannr: u32,
    /// Target/Decoy label, -1 is decoy, 1 is target
    pub label: i32,
    /// Experimental mass
    pub expmass: f32,
    /// Calculated mass
    pub calcmass: f32,
    /// Reported precursor charge
    pub charge: u8,
    /// Retention time
    pub rt: f32,
    /// Predicted RT, if enabled
    pub predicted_rt: f32,
    /// Difference between predicted & observed RT
    pub delta_rt: f32,
    /// Difference between expmass and calcmass
    pub delta_mass: f32,
    /// C13 isotope error
    pub isotope_error: f32,
    /// Average ppm delta mass for matched fragments
    pub average_ppm: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f64,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_hyperscore: f64,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    /// Longest b-ion series
    pub longest_b: usize,
    /// Longest y-ion series
    pub longest_y: usize,
    /// Longest y-ion series, divided by peptide length
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: usize,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f64,
    /// Combined score from linear discriminant analysis, used for FDR calc
    pub discriminant_score: f32,
    /// Posterior error probability for this PSM / local FDR
    pub posterior_error: f32,
    /// Assigned q_value
    pub q_value: f32,

    pub ms2_intensity: f32,
    pub ms1_intensity: f32,
    pub ms1_apex: f32,
    pub ms1_apex_rt: f32,
}

impl Score {
    /// Calculate the X!Tandem hyperscore
    /// * `fact_table` is a precomputed vector of factorials
    fn hyperscore(&self, fact_table: &[f64]) -> f64 {
        let i = (self.summed_b + 1.0) as f64 * (self.summed_y + 1.0) as f64;
        let m = fact_table[(self.matched_b as usize).min(fact_table.len() - 2)].ln()
            + fact_table[(self.matched_y as usize).min(fact_table.len() - 2)].ln();

        let score = i.ln() + m;
        if score.is_finite() {
            score
        } else {
            255.0
        }
    }
}

pub struct Scorer<'db> {
    pub db: &'db IndexedDatabase,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    /// Precursor isotope error lower bounds (e.g. -1)
    min_isotope_err: i8,
    /// Precursor isotope error upper bounds (e.g. 3)
    max_isotope_err: i8,
    max_fragment_charge: Option<u8>,
    min_fragment_mass: f32,
    max_fragment_mass: f32,
    factorial: [f64; 32],
    chimera: bool,
}

impl<'db> Scorer<'db> {
    pub fn new(
        db: &'db IndexedDatabase,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
        min_isotope_err: i8,
        max_isotope_err: i8,
        max_fragment_charge: Option<u8>,
        min_fragment_mass: f32,
        max_fragment_mass: f32,
        chimera: bool,
    ) -> Self {
        let mut factorial = [1.0f64; 32];
        for i in 1..32 {
            factorial[i] = factorial[i - 1] * i as f64;
        }

        debug_assert!(factorial[3] == 6.0);
        debug_assert!(factorial[4] == 24.0);

        Scorer {
            db,
            precursor_tol,
            fragment_tol,
            min_isotope_err,
            max_isotope_err,
            max_fragment_charge,
            min_fragment_mass,
            max_fragment_mass,
            factorial,
            chimera,
        }
    }

    pub fn score(&self, query: &ProcessedSpectrum, report_psms: usize) -> Vec<Percolator> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        match self.chimera {
            true => self.score_chimera(query),
            false => self.score_standard(query, report_psms),
        }
    }

    /// Preliminary Score, return # of matched peaks per candidate, sorted high to low
    fn matched_peaks(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        charge: u8,
    ) -> (usize, Vec<PreScore>) {
        let candidates = self.db.query(
            precursor_mass,
            self.precursor_tol,
            self.fragment_tol,
            self.min_isotope_err,
            self.max_isotope_err,
        );

        // Allocate space for all potential candidates - many potential candidates
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut score_vector = vec![PreScore::default(); potential];

        let mut matches = 0;
        let mut scored_candidates = 0;
        for peak in query.peaks.iter() {
            for charge in 1..charge {
                let mass = peak.mass * charge as f32;
                for frag in candidates.page_search(mass) {
                    let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                    let mut sc = &mut score_vector[idx];
                    if sc.matched == 0 {
                        scored_candidates += 1;
                    }
                    sc.peptide = frag.peptide_index;
                    sc.matched += 1;
                    matches += 1;
                }
            }
        }
        if matches == 0 {
            return (matches, Vec::new());
        }
        score_vector.sort_unstable_by(|a, b| b.matched.cmp(&a.matched));
        score_vector.truncate(scored_candidates);
        (matches, score_vector)
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum, report_psms: usize) -> Vec<Percolator> {
        let (precursor_mass, precursor_charge) = query
            .extract_ms1_precursor()
            .expect("missing MS1 precursor");

        // If user has configured max_fragment_charge, potentially override precursor
        // charge
        let charge = precursor_charge.min(
            self.max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        );

        let (matches, preliminary) = self.matched_peaks(query, precursor_mass, charge);

        let n_calculate = 50.max(report_psms * 2).min(preliminary.len());
        let mut score_vector = preliminary
            .iter()
            .take(n_calculate)
            .map(|pre| self.score_candidate(query, charge, pre.peptide))
            .collect::<Vec<_>>();

        score_vector.sort_unstable_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));

        let mut reporting = Vec::new();

        // Expected value for poisson distribution
        let lambda = matches as f64 / preliminary.len() as f64;

        for idx in 0..report_psms.min(score_vector.len()) {
            let better = score_vector[idx];

            let next = score_vector
                .get(idx + 1)
                .map(|score| score.hyperscore)
                .unwrap_or_default();

            let peptide = &self.db[better.peptide];

            // Poisson distribution probability mass function
            let k = (better.matched_b + better.matched_y) as usize;
            let mut poisson = lambda.powi(k as i32) * f64::exp(-lambda)
                / self.factorial[k.min(self.factorial.len() - 1)];

            if poisson.is_infinite() {
                poisson = 1E-325;
            }

            let (b, y) = self.rescore(query, charge, peptide);

            let mut isotope_error = 0.0;
            for i in self.min_isotope_err..=self.max_isotope_err {
                let c13 = i as f32 * NEUTRON;
                let (iso_tol_lo, iso_tol_hi) = self.precursor_tol.bounds(precursor_mass - c13);
                if peptide.monoisotopic >= iso_tol_lo && peptide.monoisotopic <= iso_tol_hi {
                    isotope_error = c13;
                }
            }

            let delta_mass = (precursor_mass - peptide.monoisotopic - isotope_error).abs() * 1E6
                / peptide.monoisotopic;

            let (num_proteins, proteins) = self.db.assign_proteins(peptide);

            reporting.push(Percolator {
                // Identifiers
                peptide_idx: better.peptide,
                peptide: peptide.to_string(),
                peptide_decoy_idx: peptide.idx,
                proteins,
                num_proteins,
                specid: String::default(),
                scannr: query.scan as u32,
                file_id: query.file_id,
                label: peptide.label(),
                expmass: precursor_mass,
                calcmass: peptide.monoisotopic,

                // Features
                charge: precursor_charge,
                rt: query.scan_start_time,
                delta_mass,
                isotope_error,
                average_ppm: better.ppm_difference / k as f32,
                hyperscore: better.hyperscore,
                delta_hyperscore: better.hyperscore - next,
                matched_peaks: k as u32,
                matched_intensity_pct: 100.0 * (better.summed_b + better.summed_y)
                    / query.total_intensity,
                poisson: poisson.log10(),
                longest_b: b,
                longest_y: y,
                longest_y_pct: y as f32 / (peptide.sequence.len() as f32),
                peptide_len: peptide.sequence.len(),
                scored_candidates: preliminary.len(),
                missed_cleavages: peptide.missed_cleavages,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                q_value: 1.0,
                predicted_rt: 0.0,
                delta_rt: 1.0,
                ms2_intensity: better.summed_b + better.summed_y,
                ms1_intensity: 0.0,
                ms1_apex: 0.0,
                ms1_apex_rt: 0.0,
            })
        }
        reporting
    }

    /// Return 2 PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed
    pub fn score_chimera(&self, query: &ProcessedSpectrum) -> Vec<Percolator> {
        let (_, charge) = query
            .extract_ms1_precursor()
            .expect("missing MS1 precursor");

        let mut scores = self.score_standard(query, 1);
        if scores.is_empty() {
            return scores;
        }

        let best = &scores[0];
        let mut subtracted = ProcessedSpectrum {
            peaks: Vec::new(),
            precursors: query
                .precursors
                .iter()
                .map(|prec| Precursor {
                    mz: prec.mz + 0.005,
                    ..*prec
                })
                .collect(),
            ..*query
        };

        let peptide = &self.db[best.peptide_idx];
        if !peptide.decoy {
            let mut theo = [Kind::B, Kind::Y]
                .iter()
                .flat_map(|kind| {
                    IonSeries::new(peptide, *kind).map(|ion| Theoretical {
                        peptide_index: PeptideIx(0),
                        fragment_mz: ion.monoisotopic_mass,
                    })
                })
                .collect::<Vec<_>>();
            theo.sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

            for peak in &query.peaks {
                let mut allow = true;
                for charge in 1..charge {
                    let (low, high) = self.fragment_tol.bounds(peak.mass * charge as f32);
                    let slice =
                        binary_search_slice(&theo, |a, x| a.fragment_mz.total_cmp(x), low, high);
                    for frag in &theo[slice.0..slice.1] {
                        if frag.fragment_mz >= low && frag.fragment_mz <= high {
                            allow = false;
                        }
                    }
                }
                if allow {
                    subtracted.peaks.push(*peak);
                }
            }

            scores.extend(self.score_standard(&subtracted, 1));
        }
        scores
    }

    /// Calculate full hyperscore for a given PSM
    fn score_candidate(
        &self,
        query: &ProcessedSpectrum,
        charge: u8,
        peptide_ix: PeptideIx,
    ) -> Score {
        let mut score = Score {
            peptide: peptide_ix,
            ..Default::default()
        };
        let peptide = &self.db[peptide_ix];

        // Regenerate theoretical ions
        let mut fragments = IonSeries::new(peptide, Kind::B)
            .chain(IonSeries::new(peptide, Kind::Y))
            .filter(|ion| {
                ion.monoisotopic_mass >= self.min_fragment_mass
                    && ion.monoisotopic_mass <= self.max_fragment_mass
            })
            .collect::<Vec<_>>();

        fragments.sort_unstable_by(|a, b| a.monoisotopic_mass.total_cmp(&b.monoisotopic_mass));

        for peak in query.peaks.iter() {
            for charge in 1..charge {
                let mass = peak.mass * charge as f32;
                let (lo, hi) = self.fragment_tol.bounds(mass);
                let window = binary_search_slice(
                    &fragments,
                    |frag, mz| frag.monoisotopic_mass.total_cmp(mz),
                    lo,
                    hi,
                );

                for frag in fragments[window.0..window.1]
                    .iter()
                    .filter(|frag| frag.monoisotopic_mass >= lo && frag.monoisotopic_mass <= hi)
                {
                    score.ppm_difference +=
                        (frag.monoisotopic_mass - mass).abs() * 1E6 / frag.monoisotopic_mass;

                    match frag.kind {
                        Kind::B => {
                            score.matched_b += 1;
                            score.summed_b += peak.intensity;
                        }
                        Kind::Y => {
                            score.matched_y += 1;
                            score.summed_y += peak.intensity;
                        }
                    }
                }
            }
        }

        score.hyperscore = score.hyperscore(&self.factorial);
        score
    }

    /// Calculate the longest continous chain of B or Y fragment ions
    fn longest_series(&self, peaks: &[Peak], charge: u8, kind: Kind, peptide: &Peptide) -> usize {
        let mut current_start = peaks.len();
        let mut run = 0;
        let mut longest_run = 0;
        'outer: for (idx, frag) in IonSeries::new(peptide, kind)
            .map(|ion| Theoretical {
                peptide_index: PeptideIx(0),
                fragment_mz: ion.monoisotopic_mass,
            })
            .enumerate()
        {
            for charge in 1..charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                if crate::spectrum::select_closest_peak(
                    peaks,
                    frag.fragment_mz / charge as f32,
                    self.fragment_tol,
                )
                .is_some()
                {
                    run += 1;
                    if current_start + run == idx {
                        longest_run = longest_run.max(run);
                        continue 'outer;
                    }
                }
            }
            current_start = idx;
            run = 0;
        }

        longest_run
    }

    /// Rescore a candidate peptide selected for final reporting:
    /// calculate the longest (b, y) continous fragment ion series
    fn rescore(
        &self,
        query: &ProcessedSpectrum,
        charge: u8,
        candidate: &Peptide,
    ) -> (usize, usize) {
        let b = self.longest_series(&query.peaks, charge, Kind::B, candidate);
        let y = self.longest_series(&query.peaks, charge, Kind::Y, candidate);
        (b, y)
    }
}

// /// Calculate percentage of MS1 signal in isolation window that is attributed
// /// to the selected precursor ion
// pub fn ion_interference(spectra: &[ProcessedSpectrum], ms2: &ProcessedSpectrum) -> Option<f32> {
//     // let (mz, charge) = ms2.extract_ms1_precursor();
//     let precursor = ms2.precursors.first()?;
//     let ms1_scan = crate::spectrum::find_spectrum_by_id(spectra, precursor.scan?)?;

//     let isolation_window = precursor
//         .isolation_window
//         .unwrap_or_else(|| Tolerance::Da(-0.5, 0.5));
//     let precursor_mass = (precursor.mz - PROTON) * precursor.charge? as f32;
//     let (lo, hi) = isolation_window.bounds(precursor_mass);

//     let (idx_lo, idx_hi) =
//         binary_search_slice(&ms1_scan.peaks, |peak, mz| peak.mass.total_cmp(mz), lo, hi);

//     let mut selected_intensity = precursor.intensity.unwrap_or(0.0);
//     let mut total_intensity = selected_intensity;
//     for peak in ms1_scan.peaks[idx_lo..idx_hi]
//         .iter()
//         .filter(|peak| peak.mass >= lo && peak.mass <= hi)
//     {
//         if peak.mass == precursor_mass {
//             selected_intensity = peak.intensity;
//         }
//         total_intensity += peak.intensity;
//     }
//     if total_intensity == 0.0 {
//         None
//     } else {
//         Some((total_intensity - selected_intensity) / total_intensity)
//     }
// }
