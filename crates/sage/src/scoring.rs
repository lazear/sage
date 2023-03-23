use std::ops::AddAssign;

use crate::database::{binary_search_slice, IndexedDatabase, PeptideIx, Theoretical};
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON, PROTON};
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
    precursor_charge: u8,
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug)]
struct PreScore {
    peptide: PeptideIx,
    matched: u16,
    precursor_charge: u8,
}

/// Store preliminary scores & stats for first pass search for a query spectrum
#[derive(Clone, Default)]
struct InitialHits {
    matched_peaks: usize,
    // Number of peptide candidates with > 0 matched peaks
    scored_candidates: usize,
    preliminary: Vec<PreScore>,
}

impl AddAssign<InitialHits> for InitialHits {
    fn add_assign(&mut self, rhs: InitialHits) {
        self.matched_peaks += rhs.matched_peaks;
        self.scored_candidates += rhs.scored_candidates;
        // Add the non-zero matches into the accumulator
        self.preliminary
            .extend(rhs.preliminary.into_iter().take(rhs.scored_candidates));
    }
}

#[derive(Serialize, Clone, Debug)]
/// Features of a candidate peptide spectrum match
pub struct Feature {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    /// Peptide sequence, including modifications e.g.: NC(+57.021)HK
    pub peptide: String,
    /// Peptide length
    pub peptide_len: usize,
    /// Proteins containing this peptide sequence
    pub proteins: String,
    /// Number of proteins assigned to this peptide sequence
    pub num_proteins: usize,
    /// Spectrum id
    pub spec_id: String,
    /// File identifier
    pub file_id: usize,
    /// PSM rank
    pub rank: u32,
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
    /// Globally aligned retention time
    pub aligned_rt: f32,
    /// Predicted RT, if enabled
    pub predicted_rt: f32,
    /// Difference between predicted & observed RT
    pub delta_rt_model: f32,
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
    pub longest_b: u32,
    /// Longest y-ion series
    pub longest_y: u32,
    /// Longest y-ion series, divided by peptide length
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: u32,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f64,
    /// Combined score from linear discriminant analysis, used for FDR calc
    pub discriminant_score: f32,
    /// Posterior error probability for this PSM / local FDR
    pub posterior_error: f32,
    /// Assigned q_value
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub protein_q: f32,

    pub ms2_intensity: f32,
    pub ms1_intensity: f32,
    // pub ms1_apex: f32,
    // pub ms1_apex_rt: f32,
}

/// Stirling's approximation for log factorial
fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

impl Score {
    /// Calculate the X!Tandem hyperscore
    /// * `fact_table` is a precomputed vector of factorials
    fn hyperscore(&self) -> f64 {
        let i = (self.summed_b + 1.0) as f64 * (self.summed_y + 1.0) as f64;
        // let m = fact_table[(self.matched_b as usize).min(fact_table.len() - 2)].ln()
        // + fact_table[(self.matched_y as usize).min(fact_table.len() - 2)].ln();

        let score = i.ln() + lnfact(self.matched_b) + lnfact(self.matched_y);
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
    /// What is the minimum number of matched b and y ion peaks to report PSMs for?
    pub min_matched_peaks: u16,
    /// Precursor isotope error lower bounds (e.g. -1)
    pub min_isotope_err: i8,
    /// Precursor isotope error upper bounds (e.g. 3)
    pub max_isotope_err: i8,
    pub max_fragment_charge: Option<u8>,
    pub min_fragment_mass: f32,
    pub max_fragment_mass: f32,
    // factorial: [f64; 32],
    pub chimera: bool,
}

impl<'db> Scorer<'db> {
    pub fn score(&self, query: &ProcessedSpectrum, report_psms: usize) -> Vec<Feature> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        match self.chimera {
            true => self.score_chimera(query),
            false => self.score_standard(query, report_psms),
        }
    }

    #[inline(always)]
    /// If user has configured max_fragment_charge, potentially override precursor
    /// charge
    fn max_fragment_charge(&self, precursor_charge: u8) -> u8 {
        precursor_charge.min(
            self.max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
    }

    /// Preliminary Score, return # of matched peaks per candidate, sorted high to low
    fn matched_peaks(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
    ) -> InitialHits {
        let candidates = self.db.query(
            precursor_mass,
            self.precursor_tol,
            self.fragment_tol,
            self.min_isotope_err,
            self.max_isotope_err,
        );

        let max_fragment_charge = self.max_fragment_charge(precursor_charge);

        // Allocate space for all potential candidates - many potential candidates
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut hits = InitialHits {
            matched_peaks: 0,
            scored_candidates: 0,
            preliminary: vec![PreScore::default(); potential],
        };

        for peak in query.peaks.iter() {
            for charge in 1..max_fragment_charge {
                let mass = peak.mass * charge as f32;
                for frag in candidates.page_search(mass) {
                    let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                    let mut sc = &mut hits.preliminary[idx];
                    if sc.matched == 0 {
                        hits.scored_candidates += 1;
                        sc.precursor_charge = precursor_charge;
                    }
                    sc.peptide = frag.peptide_index;
                    sc.matched += 1;
                    hits.matched_peaks += 1;
                }
            }
        }
        if hits.matched_peaks == 0 {
            return hits;
        }
        // Sort the preliminary scoring vector from high to low number of matched peaks
        hits.preliminary
            .sort_unstable_by(|a, b| b.matched.cmp(&a.matched));
        hits
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum, report_psms: usize) -> Vec<Feature> {
        let precursor = query.precursors.get(0).unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        let hits = if let Some(charge) = precursor.charge {
            // Charge state is already annotated for this precusor, only search once
            let precursor_mass = mz * charge as f32;
            self.matched_peaks(query, precursor_mass, charge)
        } else {
            // Not all selected ion precursors have charge states annotated -
            // assume it could be z=2, z=3, z=4 and search all three
            let mut hits = (2..5).fold(InitialHits::default(), |mut hits, precursor_charge| {
                let precursor_mass = mz * precursor_charge as f32;
                hits += self.matched_peaks(query, precursor_mass, precursor_charge);
                hits
            });

            // We have merged results from multiple initial searches,
            // which are now concatenated and need to be sorted again
            hits.preliminary.sort_by(|a, b| b.matched.cmp(&a.matched));
            hits
        };

        // Determine how many candidates to actually calculate hyperscore for.
        // Hyperscore is relatively computationally expensive, so we don't want
        // to calculate it for every possible candidate (100s - 10,000s depending on search)
        // when we are only going to report a few PSMs. But we also want to calculate
        // it for enough candidates that we don't accidentally miss the best hit!
        //
        // Given that hyperscore is dominated by the number of matched peaks, it seems
        // reasonable to assume that the highest hyperscore will belong to one of the
        // top 50 candidates sorted by # of matched peaks.
        let n_calculate = 50.clamp(
            (report_psms * 2).min(hits.preliminary.len()),
            hits.preliminary.len(),
        );
        let mut score_vector = hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
            .take(n_calculate)
            .map(|pre| self.score_candidate(query, pre.precursor_charge, pre.peptide))
            .filter(|s| (s.matched_b + s.matched_y) >= self.min_matched_peaks)
            .collect::<Vec<_>>();

        // Hyperscore is our primary score function for PSMs
        score_vector.sort_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));

        // Expected value for poisson distribution
        // (average # of matches peaks/peptide candidate)
        let lambda = hits.matched_peaks as f64 / hits.scored_candidates as f64;

        let mut reporting = Vec::new();
        for idx in 0..report_psms.min(score_vector.len()) {
            let better = score_vector[idx];
            let peptide = &self.db[better.peptide];
            let precursor_mass = mz * better.precursor_charge as f32;

            let next = score_vector
                .get(idx + 1)
                .map(|score| score.hyperscore)
                .unwrap_or_default();

            // Poisson distribution probability mass function
            let k = better.matched_b + better.matched_y;
            let mut poisson = lambda.powi(k as i32) * f64::exp(-lambda) / lnfact(k).exp();

            if poisson.is_infinite() {
                // Approximately the smallest positive non-zero value representable by f64
                poisson = 1E-325;
            }

            // Calculate the longest continuous b- and y-ion ladders
            let (b, y) = self.rescore(query, better.precursor_charge, peptide);

            let mut isotope_error = 0.0;
            for i in self.min_isotope_err..=self.max_isotope_err {
                let c13 = i as f32 * NEUTRON;
                let (iso_tol_lo, iso_tol_hi) = self.precursor_tol.bounds(precursor_mass - c13);
                if peptide.monoisotopic >= iso_tol_lo && peptide.monoisotopic <= iso_tol_hi {
                    isotope_error = c13;
                    break;
                }
            }

            let delta_mass = (precursor_mass - peptide.monoisotopic - isotope_error).abs() * 1E6
                / peptide.monoisotopic;

            let (num_proteins, proteins) = self.db.assign_proteins(peptide);

            reporting.push(Feature {
                // Identifiers
                peptide_idx: better.peptide,
                peptide: peptide.to_string(),
                proteins,
                num_proteins,
                spec_id: query.id.clone(),
                file_id: query.file_id,
                rank: idx as u32 + 1,
                label: peptide.label(),
                expmass: precursor_mass,
                calcmass: peptide.monoisotopic,
                // Features
                charge: better.precursor_charge,
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
                scored_candidates: hits.scored_candidates as u32,
                missed_cleavages: peptide.missed_cleavages,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                spectrum_q: 1.0,
                protein_q: 1.0,
                peptide_q: 1.0,
                predicted_rt: 0.0,
                aligned_rt: query.scan_start_time,
                delta_rt_model: 0.999,
                ms2_intensity: better.summed_b + better.summed_y,
                ms1_intensity: precursor.intensity.unwrap_or(0.0),
            })
        }
        reporting
    }

    /// Return 2 PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed
    pub fn score_chimera(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
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
                    spectrum_ref: prec.spectrum_ref.clone(),
                    ..*prec
                })
                .collect(),
            id: query.id.clone(),
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
                for charge in 1..best.charge {
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
        precursor_charge: u8,
        peptide_ix: PeptideIx,
    ) -> Score {
        let mut score = Score {
            peptide: peptide_ix,
            precursor_charge,
            ..Default::default()
        };
        let peptide = &self.db[peptide_ix];
        let max_fragment_charge = self.max_fragment_charge(precursor_charge);

        // Regenerate theoretical ions - initial database search might be
        // using only a subset of all possible ions (e.g. no b1/b2/y1/y2)
        // so we need to completely re-score this candidate
        let mut fragments = IonSeries::new(peptide, Kind::B)
            .chain(IonSeries::new(peptide, Kind::Y))
            .filter(|ion| {
                ion.monoisotopic_mass >= self.min_fragment_mass
                    && ion.monoisotopic_mass <= self.max_fragment_mass
            })
            .collect::<Vec<_>>();

        fragments.sort_unstable_by(|a, b| a.monoisotopic_mass.total_cmp(&b.monoisotopic_mass));

        for peak in query.peaks.iter() {
            for charge in 1..max_fragment_charge {
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

        score.hyperscore = score.hyperscore();
        score
    }

    /// Calculate the longest continuous chain of B or Y fragment ions
    fn longest_series(
        &self,
        peaks: &[Peak],
        max_fragment_charge: u8,
        kind: Kind,
        peptide: &Peptide,
    ) -> u32 {
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
            for charge in 1..max_fragment_charge {
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
                        longest_run = longest_run.max(run as u32);
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
        precursor_charge: u8,
        candidate: &Peptide,
    ) -> (u32, u32) {
        let max_fragment_charge = self.max_fragment_charge(precursor_charge);
        let b = self.longest_series(&query.peaks, max_fragment_charge, Kind::B, candidate);
        let y = self.longest_series(&query.peaks, max_fragment_charge, Kind::Y, candidate);
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
