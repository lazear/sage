//! TMT quantification

use crate::database::{binary_search_slice, Theoretical};
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, PROTON};
use crate::peptide::Peptide;
use crate::scoring::{Percolator, Scorer};
use crate::spectrum::{self, Peak, ProcessedSpectrum};

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Purity {
    pub ratio: f32,
    pub correct_precursors: usize,
    pub incorrect_precursors: usize,
}

fn purity_of_match<'a>(
    ms3: &ProcessedSpectrum,
    ms2: &ProcessedSpectrum,
    theoretical_peaks: &[Peak],
    max_charge: u8,
    isolation_tolerance: Tolerance,
    fragment_tolerance: Tolerance,
) -> Purity {
    let mut explained_intensity = 0.0;
    let mut interference = 0.0;
    let mut correct_precursors = 0;
    let mut incorrect_precursors = 0;

    for precursor in ms3.precursors.iter() {
        // Spurious SPS precursor introduced by MSConvert
        if precursor.scan.unwrap_or(0) != ms2.scan {
            // return None;
            continue;
        }

        // This is really a m/z window, we haven't performed charge state deconvolution!
        // Select a window of MS2 peaks that have been sampled for MS3
        let isolation_window = isolation_tolerance.bounds(precursor.mz - PROTON);

        let (idx_lo, idx_hi) = binary_search_slice(
            &ms2.peaks,
            |peak, key| peak.mass.total_cmp(key),
            isolation_window.0,
            isolation_window.1,
        );

        // Create slice surrounding SPS peaks selected
        let window = &ms2.peaks[idx_lo..idx_hi];
        interference += window
            .iter()
            .filter(|peak| peak.mass >= isolation_window.0 && peak.mass <= isolation_window.1)
            .map(|peak| peak.intensity)
            .sum::<f32>();

        // Pick the most intense peak within the isolation window, and decide whether it is
        // assignable to the best candidate peptide
        let assigned_to_candidate = if let Some(best_peak) =
            spectrum::select_closest_peak(&window, precursor.mz - PROTON, isolation_tolerance)
        {
            if ms3.scan == 3008 {
                dbg!(precursor.mz - PROTON, best_peak);
                dbg!(&theoretical_peaks);
            }

            let correct = (1..max_charge).any(|charge| {
                let mass = best_peak.mass * charge as f32;
                // Theoretical peaks don't have intensity assigned, so pick closest peak
                spectrum::select_closest_peak(&theoretical_peaks, mass, fragment_tolerance)
                    .is_some()
            });
            if correct {
                interference -= best_peak.intensity;
                explained_intensity += best_peak.intensity;
            };
            correct
        } else {
            eprintln!("couldn't select a peak within window - maybe deisotoped out of existence!!");
            dbg!(&ms2.peaks);
            false
        };

        if assigned_to_candidate {
            correct_precursors += 1;
        } else {
            incorrect_precursors += 1;
        }
    }

    Purity {
        ratio: explained_intensity / (explained_intensity + interference),
        correct_precursors,
        incorrect_precursors,
    }
}

fn mk_theoretical(peptide: &Peptide) -> Vec<Peak> {
    let mut theoretical_peaks = IonSeries::new(peptide, Kind::B)
        .map(|ion| Peak {
            mass: ion.monoisotopic_mass,
            intensity: 0.0,
        })
        .collect::<Vec<_>>();

    if let Some(c_terminal) = peptide.sequence.last() {
        if let crate::mass::Residue::Mod('K', _) = c_terminal {
            theoretical_peaks.extend(IonSeries::new(peptide, Kind::Y).map(|ion| Peak {
                mass: ion.monoisotopic_mass,
                intensity: 0.0,
            }));
        }
    }

    theoretical_peaks.sort_unstable_by(|a, b| a.mass.total_cmp(&b.mass));
    theoretical_peaks
}

#[derive(Debug)]
pub struct Quant<'a, 'b> {
    pub hit: Percolator<'a>,
    pub chimera: Option<Percolator<'a>>,
    pub hit_purity: Purity,
    pub chimera_purity: Option<Purity>,
    pub intensities: Vec<Option<&'b Peak>>,
}

fn quantify<'a>(
    ms3: &'a ProcessedSpectrum,
    labels: &[f32],
    label_tolerance: Tolerance,
) -> Vec<Option<&'a Peak>> {
    labels
        .iter()
        .map(|label| {
            spectrum::select_most_intense_peak(&ms3.peaks, label - PROTON, label_tolerance)
        })
        // .map(|label| spectrum::select_closest_peak(&ms3.peaks, label - PROTON, label_tolerance))
        .collect()
}

const TMT11PLEX: [f32; 11] = [
    126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
    130.141145, 131.138180, 131.144499,
];

pub fn quantify_sps<'a, 'b>(
    scorer: &'a Scorer<'a>,
    spectra: &[ProcessedSpectrum],
    ms3: &'b ProcessedSpectrum,
) -> Option<Quant<'a, 'b>> {
    let first_precursor = ms3
        .precursors
        .first()
        .expect("MS3 scan without at least one precursor!");

    let ms2 = spectrum::find_spectrum_by_id(
        spectra,
        first_precursor
            .scan
            .expect("MS3 scan without a MS2 precursor scan ID"),
    )
    .expect("Couldn't locate parent MS2 scan!");

    let ms2_charge = ms2
        .precursors
        .get(0)
        .and_then(|p| p.charge)
        .unwrap_or(2)
        .saturating_sub(1);
    let scores = scorer.score(&ms2, 1);
    let hit = scores.first()?.clone();
    let peptide = &scorer.db[hit.peptide_idx];
    let hit_purity = purity_of_match(
        ms3,
        ms2,
        &mk_theoretical(peptide.peptide()),
        ms2_charge,
        Tolerance::Da(-1.0, 1.0),
        Tolerance::Da(-1.0, 1.0),
    );

    // let chimera = scores.get(1).cloned();
    // let chimera_purity = match &chimera {
    //     None => None,
    //     Some(score) => Some(purity_of_match(
    //         ms3,
    //         ms2,
    //         &mk_theoretical(&scorer.db[score.peptide_idx].peptide()),
    //         2,
    //         Tolerance::Da(-1.0, 1.0),
    //         Tolerance::Da(-1.0, 1.0)
    //     )),
    // };

    Some(Quant {
        hit,
        hit_purity,
        chimera: None,
        chimera_purity: None,
        intensities: quantify(ms3, &TMT11PLEX, Tolerance::Ppm(-25.0, 25.0)),
    })
}
