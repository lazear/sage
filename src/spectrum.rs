use crate::mass::{Tolerance, PROTON};

pub struct SpectrumProcessor {
    take_top_n: usize,
    max_fragment_charge: u8,
    max_fragment_mz: f32,
}

#[derive(Clone, Default, Debug)]
pub struct ProcessedSpectrum {
    pub scan: u32,
    pub monoisotopic_mass: f32,
    pub charge: u8,
    pub rt: f32,
    pub peaks: Vec<(f32, f32)>,
}

impl SpectrumProcessor {
    pub fn new(take_top_n: usize, max_fragment_charge: u8, max_fragment_mz: f32) -> Self {
        Self {
            take_top_n,
            max_fragment_charge,
            max_fragment_mz,
        }
    }

    pub fn process(&self, s: crate::mzml::Spectrum) -> Option<ProcessedSpectrum> {
        if s.ms_level != 2 {
            return None;
        }

        let charge = self
            .max_fragment_charge
            .min(s.precursor_charge?.saturating_sub(1).max(1));

        let precursor = (s.precursor_mz? - PROTON) * s.precursor_charge? as f32;

        let (prec_lo, prec_hi) = Tolerance::Ppm(-1.5, 1.5).bounds(precursor);

        let mut peaks =
            s.mz.into_iter()
                .zip(s.intensity.into_iter())
                .collect::<Vec<_>>();

        peaks.sort_unstable_by(|(_, a), (_, b)| b.total_cmp(a));
        let n = peaks.len().min(self.take_top_n);
        let peaks = peaks[..n]
            .iter()
            .flat_map(|(mz, int)| {
                (1..=charge).filter_map(move |charge| {
                    // OK, this bit is kinda weird - to save memory, instead of calculating theoretical
                    // m/z's for different charge states, we instead resample the experimental spectra
                    //
                    // We assume that the m/z we are observing are potentially at charge state >1, and we
                    // want to convert them to charge state = 1 (well, really neutral monoisotopic mass)
                    // in order to simulate a calculated theoretical fragment with a higher charge
                    let fragment_mass = (mz - PROTON) * charge as f32;
                    match (fragment_mass <= self.max_fragment_mz)
                        && (fragment_mass > prec_hi || fragment_mass < prec_lo)
                    {
                        true => Some((fragment_mass, int.sqrt())),
                        false => None,
                    }
                })
            })
            .collect::<Vec<_>>();

        Some(ProcessedSpectrum {
            scan: s.scan_id as u32,
            monoisotopic_mass: precursor,
            charge: s.precursor_charge?,
            rt: s.scan_start_time,
            peaks,
        })
    }
}
