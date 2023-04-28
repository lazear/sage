use serde::Serialize;

use crate::database::binary_search_slice;
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::mzml::{Representation, Spectrum};

/// A charge-less peak at monoisotopic mass
#[derive(PartialEq, PartialOrd, Copy, Clone, Default, Debug, Serialize)]
pub struct Peak {
    pub mass: f32,
    pub intensity: f32,
}

/// A de-isotoped peak, that might have some charge state information
#[derive(PartialEq, PartialOrd, Debug, Copy, Clone)]
pub struct Deisotoped {
    pub mz: f32,
    // Cumulative intensity of all isotopic peaks in the envelope higher than this one
    pub intensity: f32,
    // Assigned charge
    pub charge: Option<u8>,
    // If `Some(idx)`, idx is the index of the parent isotopic envelope
    pub envelope: Option<usize>,
}

pub struct SpectrumProcessor {
    pub take_top_n: usize,
    pub max_fragment_mz: f32,
    pub min_fragment_mz: f32,
    pub deisotope: bool,
    pub file_id: usize,
}

#[derive(Default, Debug, Clone, Serialize)]
pub struct Precursor {
    pub mz: f32,
    pub intensity: Option<f32>,
    pub charge: Option<u8>,
    // pub scan: Option<usize>,
    pub spectrum_ref: Option<String>,
    pub isolation_window: Option<Tolerance>,
}

#[derive(Clone, Default, Debug, Serialize)]
pub struct ProcessedSpectrum {
    /// MSn level
    pub level: u8,
    /// Scan ID
    pub id: String,
    // pub scan: usize,
    /// File ID
    pub file_id: usize,
    /// Retention time
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Selected ions for precursors, if `level > 1`
    pub precursors: Vec<Precursor>,
    /// MS peaks, sorted by mass in ascending order
    pub peaks: Vec<Peak>,
    /// Total ion current
    pub total_intensity: f32,
}

/// Linear search for most intense peak
pub fn select_most_intense_peak(peaks: &[Peak], mz: f32, tolerance: Tolerance) -> Option<&Peak> {
    let (lo, hi) = tolerance.bounds(mz);
    let mut best_peak = None;
    let mut max_int = 0.0;
    for peak in peaks
        .iter()
        .filter(|peak| peak.mass >= lo && peak.mass <= hi)
    {
        if peak.intensity >= max_int {
            max_int = peak.intensity;
            best_peak = Some(peak);
        }
    }
    best_peak
}

/// Binary search followed by linear search to select the closest peak to `mz` within `tolerance` window
pub fn select_closest_peak(peaks: &[Peak], mz: f32, tolerance: Tolerance) -> Option<&Peak> {
    let (lo, hi) = tolerance.bounds(mz);
    let (i, j) = binary_search_slice(peaks, |peak, query| peak.mass.total_cmp(query), lo, hi);

    let mut best_peak = None;
    let mut min_eps = f32::MAX;
    for peak in peaks[i..j]
        .iter()
        .filter(|peak| peak.mass >= lo && peak.mass <= hi)
    {
        let eps = (peak.mass - mz).abs();
        if eps <= min_eps {
            min_eps = eps;
            best_peak = Some(peak);
        }
    }
    best_peak
}

// pub fn find_spectrum_by_id(
//     spectra: &[ProcessedSpectrum],
//     scan_id: usize,
// ) -> Option<&ProcessedSpectrum> {
//     // First try indexing by scan
//     if let Some(first) = spectra.get(scan_id.saturating_sub(1)) {
//         if first.scan == scan_id {
//             return Some(first);
//         }
//     }
//     // Fall back to binary search
//     let idx = spectra
//         .binary_search_by(|spec| spec.scan.cmp(&scan_id))
//         .ok()?;
//     spectra.get(idx)
// }

/// Deisotope a set of peaks by attempting to find C13 peaks under a given `ppm` tolerance
pub fn deisotope(mz: &[f32], int: &[f32], max_charge: u8, ppm: f32) -> Vec<Deisotoped> {
    let mut peaks = mz
        .iter()
        .zip(int.iter())
        .map(|(mz, int)| Deisotoped {
            mz: *mz,
            intensity: *int,
            envelope: None,
            charge: None,
        })
        .collect::<Vec<_>>();

    // Is the peak at index `i` an isotopic peak?
    for i in (0..mz.len()).rev() {
        // Two pointer approach, j is fast pointer
        let mut j = i.saturating_sub(1);
        while mz[i] - mz[j] <= NEUTRON + Tolerance::ppm_to_delta_mass(mz[i], ppm) {
            let delta = mz[i] - mz[j];
            let tol = Tolerance::ppm_to_delta_mass(mz[i], ppm);
            for charge in 1..=max_charge {
                let iso = NEUTRON / charge as f32;
                if (delta - iso).abs() <= tol && int[i] < int[j] {
                    // Make sure this peak isn't already part of an isotopic envelope
                    if let Some(existing) = peaks[i].charge {
                        if existing != charge {
                            continue;
                        }
                    }
                    peaks[j].intensity += peaks[i].intensity;
                    peaks[j].charge = Some(charge);
                    peaks[i].charge = Some(charge);
                    peaks[i].envelope = Some(j);
                }
            }
            j = j.saturating_sub(1);
            if j == 0 {
                break;
            }
        }
    }
    peaks
}

/// Path compression of isotopic envelope links
pub fn path_compression(peaks: &mut [Deisotoped]) {
    for idx in 0..peaks.len() {
        if let Some(parent) = peaks[idx].envelope {
            if let Some(upper) = peaks[parent].envelope {
                peaks[idx].envelope = Some(upper);
            }
            peaks[idx].intensity = 0.0;
        }
    }
}

impl ProcessedSpectrum {
    pub fn extract_ms1_precursor(&self) -> Option<(f32, u8)> {
        let precursor = self.precursors.get(0)?;
        let charge = precursor.charge?;
        let mass = (precursor.mz - PROTON) * charge as f32;
        Some((mass, charge))
    }

    pub fn in_isolation_window(&self, mz: f32) -> Option<bool> {
        let precursor = self.precursors.get(0)?;
        let (lo, hi) = precursor.isolation_window?.bounds(precursor.mz - PROTON);
        Some(mz >= lo && mz <= hi)
    }
}

impl SpectrumProcessor {
    /// Create a new [`SpectrumProcessor`]
    ///
    /// # Arguments
    /// * `take_top_n`: Keep only the top N most intense peaks from the spectrum
    /// * `min_fragment_mz`: Keep only fragments >= this m/z
    /// * `max_fragment_mz`: Keep only fragments <= this m/z
    /// * `deisotope`: Perform deisotoping & charge state deconvolution
    /// * `file_id`: Store this value in all [`ProcessedSpectrum`]
    pub fn new(
        take_top_n: usize,
        min_fragment_mz: f32,
        max_fragment_mz: f32,
        deisotope: bool,
        file_id: usize,
    ) -> Self {
        Self {
            take_top_n,
            min_fragment_mz,
            max_fragment_mz,
            deisotope,
            file_id,
        }
    }

    fn process_ms2(&self, should_deisotope: bool, spectrum: &Spectrum) -> Vec<Peak> {
        if spectrum.representation != Representation::Centroid {
            // Panic, because there's really nothing we can do with profile data
            panic!(
                "Scan {} contains profile data! Please convert to centroid",
                spectrum.id
            );
        }

        // If there is no precursor charge from the mzML file, then deisotope fragments up to z=3
        let charge = spectrum
            .precursors
            .get(0)
            .and_then(|p| p.charge)
            .unwrap_or(3);

        if should_deisotope {
            let mut peaks = deisotope(&spectrum.mz, &spectrum.intensity, charge, 10.0);
            peaks.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));

            peaks
                .into_iter()
                .filter(|peak| {
                    peak.envelope.is_none()
                        && peak.mz >= self.min_fragment_mz
                        && peak.mz <= self.max_fragment_mz
                })
                .map(|peak| {
                    // Convert from MH* to M
                    let mass = (peak.mz - PROTON) * peak.charge.unwrap_or(1) as f32;
                    Peak {
                        mass,
                        intensity: peak.intensity,
                    }
                })
                .take(self.take_top_n)
                .collect::<Vec<Peak>>()
        } else {
            let mut peaks = spectrum
                .mz
                .iter()
                .zip(spectrum.intensity.iter())
                .filter(|&(mz, _)| *mz >= self.min_fragment_mz && *mz <= self.max_fragment_mz)
                .map(|(mz, &intensity)| {
                    let mass = (mz - PROTON) * 1.0;
                    Peak { mass, intensity }
                })
                .collect::<Vec<_>>();
            peaks.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
            peaks.truncate(self.take_top_n);
            peaks
        }
    }

    pub fn process(&self, spectrum: Spectrum) -> ProcessedSpectrum {
        let mut peaks = match spectrum.ms_level {
            2 => self.process_ms2(self.deisotope, &spectrum),
            _ => spectrum
                .mz
                .iter()
                .zip(spectrum.intensity.iter())
                .map(|(mz, &intensity)| {
                    let mass = (mz - PROTON) * 1.0;
                    Peak { mass, intensity }
                })
                .collect::<Vec<_>>(),
        };

        peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
        let total_intensity = peaks.iter().map(|peak| peak.intensity).sum::<f32>();

        ProcessedSpectrum {
            level: spectrum.ms_level,
            id: spectrum.id,
            file_id: self.file_id,
            scan_start_time: spectrum.scan_start_time,
            ion_injection_time: spectrum.ion_injection_time,
            precursors: spectrum.precursors,
            peaks,
            total_intensity,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_deisotope() {
        let mut mz = [
            800.9,
            803.4080,
            804.4108,
            805.4106,
            806.4116,
            810.0,
            812.0,
            812.0 + NEUTRON / 2.0,
        ];
        let mut int = [1., 4., 3., 2., 1., 1., 9.0, 4.5];
        let mut peaks = deisotope(&mut mz, &mut int, 2, 5.0);

        assert_eq!(
            peaks,
            vec![
                Deisotoped {
                    mz: 800.9,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 803.4080,
                    intensity: 10.0,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 804.4108,
                    intensity: 6.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 805.4106,
                    intensity: 3.0,
                    charge: Some(1),
                    envelope: Some(2),
                },
                Deisotoped {
                    mz: 806.4116,
                    intensity: 1.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 810.0,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0,
                    intensity: 13.5,
                    charge: Some(2),
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0 + NEUTRON / 2.0,
                    intensity: 4.5,
                    charge: Some(2),
                    envelope: Some(6),
                }
            ]
        );

        path_compression(&mut peaks);
        assert_eq!(
            peaks,
            vec![
                Deisotoped {
                    mz: 800.9,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 803.4080,
                    intensity: 10.0,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 804.4108,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 805.4106,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 806.4116,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 810.0,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0,
                    intensity: 13.5,
                    charge: Some(2),
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0 + NEUTRON / 2.0,
                    intensity: 0.0,
                    charge: Some(2),
                    envelope: Some(6),
                }
            ]
        );
    }
}
