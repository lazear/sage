use crate::database::{binary_search_slice, IndexedDatabase, PeptideIx, Theoretical};
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, PROTON};
use crate::peptide::TargetDecoy;
use crate::spectrum::{Deisotoped, ProcessedSpectrum, SpectrumProcessor};
use serde::Serialize;

/// Structure to hold temporary scores
#[derive(Copy, Clone, Default)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    matched_nl: u16,
    summed_b: f32,
    summed_y: f32,
    hyperscore: f32,
    ppm_difference: f32,
}

#[derive(Serialize, Debug)]
/// Features of a candidate peptide spectrum match
pub struct Percolator<'db> {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    /// Peptide sequence, including modifications e.g.: NC(+57.021)HK
    pub peptide: String,
    /// Peptide length
    pub peptide_len: usize,
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
    pub average_ppm: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f32,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_hyperscore: f32,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    pub matched_neutral_loss: u32,
    /// Longest b-ion series
    pub longest_b: usize,
    /// Longest y-ion series
    pub longest_y: usize,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: usize,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f32,

    pub discriminant_score: f32,

    pub posterior_error: f32,
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
            255.0
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
    chimera: bool,
}

impl<'db> Scorer<'db> {
    pub fn new(
        db: &'db IndexedDatabase,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
        min_isotope_err: i8,
        max_isotope_err: i8,
        chimera: bool,
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
            chimera,
        }
    }

    pub fn score(&self, query: &ProcessedSpectrum, report_psms: usize) -> Vec<Percolator> {
        match self.chimera {
            true => self.score_chimera(query),
            false => self.score_standard(query, report_psms),
        }
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(
        &self,
        query: &ProcessedSpectrum,
        report_psms: usize,
    ) -> Vec<Percolator<'db>> {
        // Create a new `IndexedQuery`

        let candidates = self.db.query(
            query,
            self.precursor_tol,
            self.fragment_tol,
            self.min_isotope_err,
            self.max_isotope_err,
        );

        // Allocate space for all potential candidates - many potential candidates
        // will not have fragments matched, so we use `Option<Score>`
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut score_vector: Vec<Score> = vec![Score::default(); potential];

        let mut total_intensity = 0.0;
        let mut matches = 0;
        for peak in query.peaks.iter() {
            total_intensity += peak.intensity;
            for charge in 1..query.charge {
                let mass = peak.mass * charge as f32;
                for frag in candidates.page_search(mass) {
                    let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                    let mut sc = &mut score_vector[idx];
                    sc.peptide = frag.peptide_index;
                    sc.ppm_difference += (frag.fragment_mz - mass).abs() * 1E6 / frag.fragment_mz;

                    match frag.kind {
                        Kind::B => {
                            sc.matched_b += 1;
                            sc.summed_b += peak.intensity;
                        }
                        Kind::Y => {
                            sc.matched_y += 1;
                            sc.summed_y += peak.intensity;
                        }
                    }

                    matches += 1;
                }
            }
            // for nl in [H2O, NH3] {
            //     let mass = peak.mass + nl;
            //     for frag in candidates.page_search(mass) {
            //         let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
            //         let mut sc = &mut score_vector[idx];
            //         match frag.kind {
            //             Kind::B => {
            //                 sc.summed_b += peak.intensity * 0.1;
            //             }
            //             Kind::Y => {
            //                 sc.summed_y += peak.intensity * 0.1;
            //             }
            //         }
            //         sc.matched_nl += 1;
            //     }
            // }
        }

        if matches == 0 {
            return Vec::new();
        }

        // Calculating hyperscore is expensive, and if we have 1500 candidates but only report 1 ...
        // we can probably get away with a fast approximation. Hyperscore & final poisson PMF are dominated by
        // the number of matched peaks - ideally the best spectrum will be a clear winner.
        // Sort our scores highest # of matched peaks to lowest, then eliminate all candidates with no matched peaks
        // Performance profiling showed that we spend about 10% of runtime just calculating hyperscore for all matches
        score_vector
            .sort_unstable_by(|a, b| (b.matched_b + b.matched_y).cmp(&(a.matched_b + a.matched_y)));
        let idx = score_vector
            .iter()
            .enumerate()
            .filter(|(_, s)| s.matched_b + s.matched_y == 0)
            .map(|(i, _)| i)
            .next()
            .unwrap_or(score_vector.len());

        score_vector.truncate(idx);
        let n = score_vector.len();

        // Calculate hyperscore for at least 50 hits (or 2x number of reported PSMs)
        // then sort that chunk of the vector - we will never actually look at the rest anyway!
        let actually_calculate = 50.max(report_psms * 2).min(n);

        for score in score_vector.iter_mut().take(actually_calculate) {
            score.hyperscore = score.hyperscore(&self.factorial);
        }

        score_vector[..actually_calculate]
            .sort_unstable_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));

        let mut reporting = Vec::new();

        let lambda = matches as f32 / score_vector.len() as f32;

        for idx in 0..report_psms.min(score_vector.len()) {
            let better = score_vector[idx];
            let next = score_vector
                .get(idx + 1)
                .map(|score| score.hyperscore)
                .unwrap_or_default();

            let peptide = self.db[better.peptide].peptide();
            let k = (better.matched_b + better.matched_y) as usize;

            // Poisson distribution probability mass function
            let poisson = lambda.powi(k as i32) * f32::exp(-lambda)
                / self.factorial[k.min(self.factorial.len() - 1)];

            let (b, y) = self.rescore(query, peptide);
            reporting.push(Percolator {
                // Identifiers
                peptide_idx: better.peptide,
                peptide: peptide.to_string(),
                proteins: &peptide.protein,
                specid: 0,
                scannr: query.scan,
                label: self.db[better.peptide].label(),
                expmass: query.monoisotopic_mass + PROTON,
                calcmass: peptide.monoisotopic + PROTON,

                // Features
                charge: query.charge,
                rt: query.rt,
                delta_mass: (query.monoisotopic_mass - peptide.monoisotopic).abs() * 1E6
                    / peptide.monoisotopic,
                average_ppm: better.ppm_difference / (better.matched_b + better.matched_y) as f32,
                hyperscore: better.hyperscore,
                delta_hyperscore: better.hyperscore - next,
                matched_peaks: k as u32,
                matched_neutral_loss: better.matched_nl as u32,
                matched_intensity_pct: (better.summed_b + better.summed_y) / total_intensity,
                poisson,
                longest_b: b,
                longest_y: y,
                peptide_len: peptide.sequence.len(),
                scored_candidates: n,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                q_value: 1.0,
            })
        }
        reporting
    }

    /// Return 2 PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed
    pub fn score_chimera(&self, query: &ProcessedSpectrum) -> Vec<Percolator<'db>> {
        let mut scores = self.score_standard(query, 1);
        if scores.is_empty() {
            return scores;
        }

        let best = &scores[0];
        let peptide = &self.db[best.peptide_idx];
        if let TargetDecoy::Target(peptide) = peptide {
            let mut subtracted = ProcessedSpectrum {
                monoisotopic_mass: query.monoisotopic_mass + 0.005,
                peaks: Vec::new(),
                ..*query
            };

            let mut theo = [Kind::B, Kind::Y]
                .iter()
                .flat_map(|kind| {
                    crate::ion_series::IonSeries::new(peptide, *kind).map(|ion| Theoretical {
                        peptide_index: PeptideIx(0),
                        fragment_mz: ion.monoisotopic_mass,
                        kind: ion.kind,
                    })
                })
                .collect::<Vec<_>>();
            theo.sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

            for peak in &query.peaks {
                let (low, high) = self.fragment_tol.bounds(peak.mass);
                let slice =
                    binary_search_slice(&theo, |a, x| a.fragment_mz.total_cmp(x), low, high);
                let mut allow = true;
                for frag in &theo[slice.0..slice.1] {
                    if frag.fragment_mz >= low && frag.fragment_mz <= high {
                        allow = false;
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

    /// Calculate the longest continous chain of B or Y fragment ions
    fn longest_series(&self, mz: &[f32], kind: Kind, peptide: &crate::peptide::Peptide) -> usize {
        let mut current_start = mz.len();
        let mut run = 0;
        let mut longest_run = 0;
        for (idx, frag) in crate::ion_series::IonSeries::new(peptide, kind)
            .map(|ion| Theoretical {
                peptide_index: PeptideIx(0),
                fragment_mz: ion.monoisotopic_mass,
                kind: ion.kind,
            })
            .enumerate()
        {
            let (lo, hi) = self.fragment_tol.bounds(frag.fragment_mz);
            let window = binary_search_slice(mz, |a, b| a.total_cmp(b), lo, hi);
            if mz[window.0..window.1]
                .iter()
                .filter(|&mz| *mz >= lo && *mz <= hi)
                .count()
                > 0
            {
                run += 1;
                if current_start + run == idx {
                    longest_run = longest_run.max(run);
                    continue;
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
        candidate: &crate::peptide::Peptide,
    ) -> (usize, usize) {
        let mut mz = query.peaks.iter().map(|peak| peak.mass).collect::<Vec<_>>();
        mz.sort_unstable_by(|a, b| a.total_cmp(b));

        let b = self.longest_series(&mz, Kind::B, candidate);
        let y = self.longest_series(&mz, Kind::Y, candidate);
        (b, y)
    }

    pub fn calculate_sps_purity(
        &self,
        processor: &crate::spectrum::SpectrumProcessor,
        spectra: &[crate::mzml::Spectrum],
        ms3_scan: &crate::mzml::Spectrum,
    ) -> Option<f32> {
        let first_precursor = ms3_scan
            .precursor
            .first()
            .expect("MS3 scan without at least one precursor!");

        let ms2 = crate::mzml::find_spectrum_by_id(
            spectra,
            first_precursor
                .scan
                .expect("MS3 scan without a MS2 precursor scan ID"),
        )
        .expect("Couldn't locate parent MS2 scan!");

        let query = processor.process(ms2).expect("Not an MS2 spectra!");

        let hit = self.score_chimera(&query);
        let chimera = hit.get(1); //.map(|p| &p.peptide);
        let hit = hit.first()?;

        let peptide = &self.db[hit.peptide_idx];

        let mut theoretical_peaks = IonSeries::new(peptide.peptide(), Kind::B).collect::<Vec<_>>();
        if let Some(c_terminal) = peptide.peptide().sequence.last() {
            match c_terminal {
                crate::mass::Residue::Mod('K', _) => {
                    theoretical_peaks.extend(IonSeries::new(peptide.peptide(), Kind::Y));
                }
                _ => {}
            }
        }

        theoretical_peaks
            .sort_unstable_by(|a, b| a.monoisotopic_mass.total_cmp(&b.monoisotopic_mass));
        // dbg!(&theoretical_peaks);

        eprintln!(
            "processing scan {}, assigned peptide {} {} [MS2 #{}] [chimera {:?} {}]",
            ms3_scan.scan_id,
            peptide.peptide(),
            hit.poisson,
            ms2.scan_id,
            chimera.map(|p| p.peptide.as_str()).unwrap_or_default(),
            chimera.map(|p| p.poisson).unwrap_or_default(),
        );

        let mut deisotoped = crate::spectrum::deisotope(&ms2.mz, &ms2.intensity, hit.charge, 5.0);
        crate::spectrum::path_compression(&mut deisotoped);

        let mut explained_intensity = 0.0;
        let mut interference = 0.0;

        for precursor in ms3_scan.precursor.iter() {
            // Spurious SPS precursor introduced by MSConvert
            if precursor.scan? != ms2.scan_id {
                // return None;
                continue;
            }

            let (tolerance_lo, tolerance_hi) = self.fragment_tol.bounds(precursor.mz);

            // This is really a m/z window, we haven't performed charge state deconvolution!
            // Select a window of MS2 peaks that have been sampled for MS3
            let isolation_window = Tolerance::Da(-1.0, 1.0).bounds(precursor.mz);

            let (idx_lo, idx_hi) = binary_search_slice(
                &deisotoped,
                |peak, key| peak.mz.total_cmp(key),
                isolation_window.0,
                isolation_window.1,
            );

            let window = &deisotoped[idx_lo..idx_hi];

            // Selected ion peak
            let mut best_peak = select_closest_peak(&window, |peak| peak.mz, precursor.mz);

            // We have already performed path compression on deisotoped peaks
            if let Some(parent) = window[best_peak].envelope {
                dbg!(window, best_peak);
                best_peak = parent - idx_lo;
            }

            // Sum of intensities in the MS3 isolation window that do not
            // correspond to either the selected ion peak, or one of it's isotopes
            interference += window
                .iter()
                .enumerate()
                .filter(|(idx, peak)| {
                    *idx != best_peak
                        && peak.envelope != Some(best_peak)
                        && peak.mz >= tolerance_lo
                        && peak.mz <= tolerance_hi
                })
                .map(|(_, peak)| peak.intensity)
                .sum::<f32>();

            // let purity = selected_intensity / summed_intensity;
            explained_intensity = window[best_peak].intensity;

            let matched = if let Some(charge) = window[best_peak].charge {
                let fragment_mass = (window[best_peak].mz - PROTON) * charge as f32;

                // Search for theoretical peaks matching this fragment
                let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(fragment_mass);

                let (idx_lo, idx_hi) = binary_search_slice(
                    &theoretical_peaks,
                    |ion, mz| ion.monoisotopic_mass.total_cmp(mz),
                    fragment_lo,
                    fragment_hi,
                );
                theoretical_peaks[idx_lo..idx_hi]
                    .iter()
                    .filter(|ion| {
                        ion.monoisotopic_mass >= isolation_window.0 - PROTON
                            && ion.monoisotopic_mass <= isolation_window.1 - PROTON
                    })
                    .count()
                    > 0
            } else {
                (1..hit.charge).any(|charge| {
                    let fragment_mass = (window[best_peak].mz - PROTON) * charge as f32;

                    // Search for theoretical peaks matching this fragment
                    let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(fragment_mass);

                    let (idx_lo, idx_hi) = binary_search_slice(
                        &theoretical_peaks,
                        |ion, mz| ion.monoisotopic_mass.total_cmp(mz),
                        fragment_lo,
                        fragment_hi,
                    );
                    let best_fragment = select_closest_peak(
                        &theoretical_peaks[idx_lo..idx_hi],
                        |theo| theo.monoisotopic_mass,
                        fragment_mass,
                    );
                    let best_fragment = theoretical_peaks[idx_lo + best_fragment];
                    eprintln!(
                        "\t\t{}/{} : {} ({:?})",
                        fragment_mass, charge, best_fragment.monoisotopic_mass, best_fragment.kind
                    );
                    best_fragment.monoisotopic_mass >= fragment_lo
                        && best_fragment.monoisotopic_mass <= fragment_hi
                })
            };

            eprintln!("\tprecursor {}; matched: {}", precursor.mz, matched);
        }

        Some(explained_intensity / (explained_intensity + interference))
    }
}

/// Linear search for closest peak - perform a binary search to narrow down
/// window first
fn select_closest_peak<T, F: Fn(&T) -> f32>(peaks: &[T], key: F, target: f32) -> usize {
    let mut min_idx = 0;
    let mut min_eps = f32::MAX;
    for (idx, peak) in peaks.iter().map(key).enumerate() {
        let eps = (target - peak).abs();
        if eps <= min_eps {
            min_eps = eps;
            min_idx = idx;
        }
    }
    min_idx
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
