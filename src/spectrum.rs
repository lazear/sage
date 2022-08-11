use crate::mass::{Tolerance, NEUTRON, PROTON};

/// A charge-less peak at monoisotopic mass
#[derive(PartialEq, PartialOrd, Copy, Clone, Default, Debug)]
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
    take_top_n: usize,
    max_fragment_mz: f32,
    min_fragment_mz: f32,
    deisotope: bool,
}

#[derive(Clone, Default, Debug)]
pub struct ProcessedSpectrum {
    pub scan: u32,
    pub monoisotopic_mass: f32,
    pub charge: u8,
    pub rt: f32,
    pub peaks: Vec<Peak>,
}

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
                    peaks[j].intensity += peaks[i].intensity;
                    // Make sure this peak isn't already part of an isotopic envelope
                    if let Some(existing) = peaks[i].charge {
                        if existing != charge {
                            continue;
                        }
                    }
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

impl SpectrumProcessor {
    pub fn new(
        take_top_n: usize,
        min_fragment_mz: f32,
        max_fragment_mz: f32,
        deisotope: bool,
    ) -> Self {
        Self {
            take_top_n,
            min_fragment_mz,
            max_fragment_mz,
            deisotope,
        }
    }

    pub fn process(&self, s: &crate::mzml::Spectrum) -> Option<ProcessedSpectrum> {
        if s.ms_level != 2 {
            return None;
        }

        let precursor = s.precursor.get(0)?;

        // Calculate bounds for clearing precursor mz
        let precursor_mz = (precursor.mz - PROTON) * precursor.charge.unwrap_or(2) as f32;
        let precursor_charge = precursor.charge.unwrap_or(2);
        let (prec_lo, prec_hi) = Tolerance::Ppm(-1.5, 1.5).bounds(precursor_mz);

        let charge = precursor_charge.saturating_sub(1).max(1);

        let mut peaks = match self.deisotope {
            true => deisotope(&s.mz, &s.intensity, charge, 5.0),
            false => {
                s.mz.iter()
                    .zip(s.intensity.iter())
                    .map(|(mz, int)| Deisotoped {
                        mz: *mz,
                        intensity: *int,
                        charge: None,
                        envelope: None,
                    })
                    .collect()
            }
        };

        peaks.sort_unstable_by(|a, b| b.intensity.total_cmp(&a.intensity));

        let peaks = peaks
            .into_iter()
            .filter(|peak| peak.envelope.is_none())
            .take(self.take_top_n)
            .filter_map(|peak| {
                // Convert from MH* to M
                let fragment_mass = (peak.mz - PROTON) * peak.charge.unwrap_or(1) as f32;
                let mass_filter = fragment_mass <= self.max_fragment_mz
                    && fragment_mass >= self.min_fragment_mz
                    && (fragment_mass < prec_lo || fragment_mass > prec_hi);
                match mass_filter {
                    true => Some(Peak {
                        mass: fragment_mass,
                        intensity: peak.intensity,
                    }),
                    false => None,
                }
            })
            .collect::<Vec<Peak>>();

        Some(ProcessedSpectrum {
            scan: s.scan_id as u32,
            monoisotopic_mass: precursor_mz,
            charge: precursor_charge,
            rt: s.scan_start_time,
            peaks,
        })
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
