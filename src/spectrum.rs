use std::collections::BTreeMap;

use crate::mass::PROTON;

pub struct SpectrumProcessor {
    take_top_n: usize,
    max_fragment_charge: u8,
    max_fragment_mz: f32,
}

#[derive(Clone, Default, Debug)]
pub struct ProcessedSpectrum {
    pub scan: u32,
    pub precursor_mz: f32,
    pub charge: u8,
    pub rt: f32,
    pub mz: Vec<f32>,
    pub int: Vec<f32>,
}

impl SpectrumProcessor {
    pub fn new(take_top_n: usize, max_fragment_charge: u8, max_fragment_mz: f32) -> Self {
        Self {
            take_top_n,
            max_fragment_charge,
            max_fragment_mz,
        }
    }

    pub fn process(&self, s: Spectrum) -> ProcessedSpectrum {
        let charge = self.max_fragment_charge.min(s.charge);
        let mut mz = Vec::with_capacity(charge as usize * self.take_top_n);
        let mut int = Vec::with_capacity(charge as usize * self.take_top_n);
        for (Intensity(fragment_int), fragment_mz) in s.peaks.iter().rev().take(self.take_top_n) {
            for charge in 1..=charge {
                let fragment_mz = (fragment_mz * charge as f32) - (charge as f32 * PROTON);
                if fragment_mz < self.max_fragment_mz {
                    mz.push(fragment_mz);
                    int.push(fragment_int.sqrt());
                }
            }
        }

        // dbg!(mz.len());

        ProcessedSpectrum {
            scan: s.scan,
            precursor_mz: s.precursor_mass,
            charge: s.charge,
            rt: s.rt,
            mz,
            int,
        }
    }
}

/// An observed MS2 spectrum
#[derive(Clone, Default, Debug)]
pub struct Spectrum {
    scan: u32,
    rt: f32,
    precursor_mass: f32,
    charge: u8,
    peaks: BTreeMap<Intensity, f32>,
}

#[derive(Default, Clone, Debug, PartialEq, PartialOrd)]
struct Intensity(f32);

impl std::cmp::Eq for Intensity {
    fn assert_receiver_is_total_eq(&self) {}
}

impl std::cmp::Ord for Intensity {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// Read a `.ms2` file
pub fn read_spectrum<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<Vec<Spectrum>> {
    let buf = std::fs::read_to_string(path)?;
    let mut spectra = Vec::new();
    let mut current = Spectrum::default();

    for line in buf.lines() {
        if line.starts_with('H') {
            continue;
        } else if line.starts_with('S') {
            if !current.peaks.is_empty() {
                spectra.push(current);
                current = Spectrum::default();
            }

            current.scan = line
                .split_ascii_whitespace()
                .nth(1)
                .unwrap()
                .parse::<u32>()
                .unwrap();
        } else if line.starts_with('Z') {
            let mut ws = line.split_whitespace().skip(1);
            current.charge = ws.next().unwrap().parse::<u8>().unwrap();
            current.precursor_mass = ws.next().unwrap().parse::<f32>().unwrap();
        } else if line.starts_with('I') {
            if let Some(stripped) = line.strip_prefix("I	RetTime ") {
                current.rt = stripped.parse::<f32>().unwrap();
            }
            continue;
        } else {
            let mut ws = line.split_whitespace().map(|s| s.parse::<f32>());
            let mz = ws.next().unwrap().unwrap();
            let abundance = ws.next().unwrap().unwrap();
            current.peaks.insert(Intensity(abundance), mz);
        }
    }

    Ok(spectra)
}
