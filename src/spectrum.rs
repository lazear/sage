use crate::mass::{Tolerance, NEUTRON, PROTON};

/// A charge-less peak at monoisotopic mass
#[derive(PartialEq, PartialOrd, Copy, Clone, Default, Debug)]
pub struct Peak {
    pub mass: f32,
    pub intensity: f32,
}

/// A de-isotoped peak, that might have some charge state information
#[derive(PartialEq, PartialOrd, Debug, Copy, Clone)]
struct Deisotoped {
    mz: f32,
    intensity: f32,
    charge: Option<u8>,
    envelope: bool,
}

pub struct SpectrumProcessor {
    take_top_n: usize,
    max_fragment_mz: f32,
}

#[derive(Clone, Default, Debug)]
pub struct ProcessedSpectrum {
    pub scan: u32,
    pub monoisotopic_mass: f32,
    pub charge: u8,
    pub rt: f32,
    pub peaks: Vec<Peak>,
}

impl SpectrumProcessor {
    pub fn new(take_top_n: usize, max_fragment_mz: f32) -> Self {
        Self {
            take_top_n,
            max_fragment_mz,
        }
    }

    /// Deisotope a set of peaks by attempting to find C13 peaks under a given `ppm` tolerance
    fn deisotope(mz: &[f32], int: &[f32], max_charge: u8, ppm: f32) -> Vec<Deisotoped> {
        let mut peaks = mz
            .iter()
            .zip(int.iter())
            .map(|(mz, int)| Deisotoped {
                mz: *mz,
                intensity: *int,
                envelope: false,
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
                        peaks[i].envelope = true;
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

    pub fn process(&self, s: crate::mzml::Spectrum) -> Option<ProcessedSpectrum> {
        if s.ms_level != 2 {
            return None;
        }

        // Calculate bounds for clearing precursor mz
        let precursor = (s.precursor_mz? - PROTON) * s.precursor_charge? as f32;
        let (prec_lo, prec_hi) = Tolerance::Ppm(-1.5, 1.5).bounds(precursor);

        let charge = s.precursor_charge.unwrap_or(2).saturating_sub(1).max(1);

        let mut peaks = Self::deisotope(&s.mz, &s.intensity, charge, 5.0);

        peaks.sort_unstable_by(|a, b| b.intensity.total_cmp(&a.intensity));

        let peaks = peaks
            .into_iter()
            .filter(|peak| !peak.envelope)
            .take(self.take_top_n)
            .filter_map(|peak| {
                // Convert from MH* to M
                let fragment_mass = (peak.mz - PROTON) * peak.charge.unwrap_or(1) as f32;
                let mass_filter = fragment_mass <= self.max_fragment_mz
                    && (fragment_mass < prec_lo || fragment_mass > prec_hi);
                match mass_filter {
                    true => Some(Peak {
                        mass: fragment_mass,
                        intensity: peak.intensity.sqrt(),
                    }),
                    false => None,
                }
            })
            .collect::<Vec<Peak>>();

        Some(ProcessedSpectrum {
            scan: s.scan_id as u32,
            monoisotopic_mass: precursor,
            charge: s.precursor_charge?,
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
        let peaks = SpectrumProcessor::deisotope(&mut mz, &mut int, 2, 5.0);
        assert_eq!(
            peaks,
            vec![
                Deisotoped {
                    mz: 800.9,
                    intensity: 1.0,
                    charge: None,
                    envelope: false
                },
                Deisotoped {
                    mz: 803.4080,
                    intensity: 10.0,
                    charge: Some(1),
                    envelope: false
                },
                Deisotoped {
                    mz: 804.4108,
                    intensity: 6.0,
                    charge: Some(1),
                    envelope: true
                },
                Deisotoped {
                    mz: 805.4106,
                    intensity: 3.0,
                    charge: Some(1),
                    envelope: true
                },
                Deisotoped {
                    mz: 806.4116,
                    intensity: 1.0,
                    charge: Some(1),
                    envelope: true
                },
                Deisotoped {
                    mz: 810.0,
                    intensity: 1.0,
                    charge: None,
                    envelope: false
                },
                Deisotoped {
                    mz: 812.0,
                    intensity: 13.5,
                    charge: Some(2),
                    envelope: false
                },
                Deisotoped {
                    mz: 812.0 + NEUTRON / 2.0,
                    intensity: 4.5,
                    charge: Some(2),
                    envelope: true
                }
            ]
        );
    }
}
