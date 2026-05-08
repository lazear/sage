//! TMT quantification
#![allow(clippy::excessive_precision)]
use crate::mass::{Tolerance, PROTON};
use crate::scoring::Feature;
use crate::spectrum::{self, ProcessedSpectrum};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum Isobaric {
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
    Tmt18,
    User(Vec<f32>),
}

impl Isobaric {
    /// Return the monoisotopic mass of reporter ions
    pub fn reporter_masses(&self) -> &[f32] {
        match self {
            Isobaric::Tmt6 => &TMT6PLEX,
            Isobaric::Tmt10 => &TMT11PLEX[0..10],
            Isobaric::Tmt11 => &TMT11PLEX,
            Isobaric::Tmt16 => &TMT18PLEX[0..16],
            Isobaric::Tmt18 => &TMT18PLEX,
            Isobaric::User(labels) => labels,
        }
    }

    /// Return the monoisotopic mass of tag
    pub fn modification_mass(&self) -> Option<f32> {
        match self {
            Isobaric::Tmt6 | Isobaric::Tmt10 | Isobaric::Tmt11 => Some(229.162932),
            Isobaric::Tmt16 => Some(304.2071),
            Isobaric::Tmt18 => Some(304.2135),
            Isobaric::User(_) => None,
        }
    }

    /// Return a column name for each tag
    pub fn headers(&self) -> Vec<String> {
        match self {
            Isobaric::User(v) => v
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("user_{}", idx + 1))
                .collect(),
            _ => self
                .reporter_masses()
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("tmt_{}", idx + 1))
                .collect(),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Purity {
    pub ratio: f32,
    pub correct_precursors: usize,
    pub incorrect_precursors: usize,
}

#[derive(Debug)]
pub struct Quant<'ms3> {
    /// Top hit for this MS3 spectrum
    pub hit: Feature,
    /// Top chimeric/co-fragmenting hit for this spectrum
    pub chimera: Option<Feature>,
    /// SPS precursor purity for the top hit
    pub hit_purity: Purity,
    /// SPS precursor purity for the chimeric hit
    pub chimera_purity: Option<Purity>,
    /// Quanitified TMT reporter ion intensities
    pub intensities: Vec<Option<f32>>,
    /// MS3 spectrum
    pub spectrum: &'ms3 ProcessedSpectrum,
}

/// Return a vector containing the peaks closest to the m/zs defined in
/// `labels`, within a given tolerance window.
/// This function is MS-level agnostic, so it can be used for either MS2 or MS3
/// quant.
pub fn find_reporter_ions(
    masses: &[f32],
    intensities: &[f32],
    labels: &[f32],
    label_tolerance: Tolerance,
) -> Vec<Option<f32>> {
    labels
        .iter()
        .map(|&label| {
            spectrum::select_most_intense_peak(
                masses,
                intensities,
                label,
                label_tolerance,
                Some(-PROTON),
            )
            .map(|idx| intensities[idx])
        })
        .collect()
}

const TMT6PLEX: [f32; 6] = [
    126.127726, 127.124761, 128.134436, 129.131471, 130.141145, 131.138180,
];

const TMT11PLEX: [f32; 11] = [
    126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
    130.141145, 131.138180, 131.144499,
];

const TMT18PLEX: [f32; 18] = [
    126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
    130.141145, 131.138180, 131.144500, 132.141535, 132.147855, 133.144890, 133.151210, 134.148245,
    134.154565, 135.15160,
];

/// Search MS/MS and quantify isobaric tag intensities from an SPS-MS3 spectrum
///
/// * `scorer`: used for searching/scoring precursor MS2 spectrum
/// * `spectra`: a slice (generally entire mzML) of spectra, that can be searched for precursor spectra
/// * `ms3`: The MS3 spectrum to search and quantify
/// * `isobaric_labels`: specify label m/zs to be used
/// * `isobaric_tolerance`: specify label tolerance
// pub fn quantify_sps<'a, 'b>(
//     scorer: &'a Scorer<'a>,
//     spectra: &[ProcessedSpectrum],
//     ms3: &'b ProcessedSpectrum,
//     isobaric_labels: &Isobaric,
//     isobaric_tolerance: Tolerance,
// ) -> Option<Quant<'b>> {
//     let first_precursor = ms3
//         .precursors
//         .first()
//         .expect("MS3 scan without at least one precursor!");

//     let ms2 = spectrum::find_spectrum_by_id(
//         spectra,
//         first_precursor
//             .scan
//             .expect("MS3 scan without a MS2 precursor scan ID"),
//     )
//     .expect("Couldn't locate parent MS2 scan!");

//     let ms1_charge = ms2
//         .precursors
//         .get(0)
//         .and_then(|p| p.charge)
//         .unwrap_or(2)
//         .saturating_sub(1);

//     let scores = scorer.score_chimera(ms2);
//     let hit = scores.first()?.clone();
//     let peptide = &scorer.db[hit.peptide_idx];
//     let hit_purity = purity_of_match(
//         &ms3.precursors,
//         ms2,
//         &mk_theoretical(peptide),
//         ms1_charge,
//         scorer.fragment_tol,
//     );

//     let chimera = scores.get(1).cloned();
//     let chimera_purity = chimera.as_ref().map(|score| {
//         purity_of_match(
//             &ms3.precursors,
//             ms2,
//             &mk_theoretical(&scorer.db[score.peptide_idx]),
//             ms1_charge,
//             scorer.fragment_tol,
//         )
//     });

//     Some(Quant {
//         hit,
//         hit_purity,
//         chimera,
//         chimera_purity,
//         intensities: find_reporter_ions(
//             &ms3.peaks,
//             isobaric_labels.reporter_masses(),
//             isobaric_tolerance,
//         ),
//         spectrum: ms3,
//     })
// }

#[derive(Clone)]
pub struct TmtQuant {
    pub spec_id: String,
    pub file_id: usize,
    pub ion_injection_time: f32,
    pub peaks: Vec<f32>,
}

/// Quantify isobaric tags from an MS2 or MS3 spectrum
///
/// * `spectra`: a slice (generally entire mzML) of spectra, that can be searched
///   for precursor spectra
/// * `isobaric_labels`: specify label m/zs to be used
/// * `isobaric_tolerance`: specify label tolerance
/// * `level`: MSn level to extract isobaric peaks from
pub fn quantify(
    spectra: &[ProcessedSpectrum],
    isobaric_labels: &Isobaric,
    isobaric_tolerance: Tolerance,
    level: u8,
) -> Vec<TmtQuant> {
    spectra
        .par_iter()
        .filter(|spectrum| spectrum.level == level)
        .filter_map(|spectrum| {
            let spec_id = match level {
                1 => return None,
                2 => spectrum.id.clone(),
                _ => spectrum
                    .precursors
                    .first()
                    .and_then(|precursor| precursor.spectrum_ref.clone())
                    .unwrap_or_default(),
            };

            let peaks = find_reporter_ions(
                &spectrum.masses,
                &spectrum.intensities,
                isobaric_labels.reporter_masses(),
                isobaric_tolerance,
            )
            .into_iter()
            .map(|peak| peak.unwrap_or_default())
            .collect();

            Some(TmtQuant {
                spec_id,
                file_id: spectrum.file_id,
                ion_injection_time: spectrum.ion_injection_time,
                peaks,
            })
        })
        .collect()
}
