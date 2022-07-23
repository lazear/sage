use crate::mass::PROTON;

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

    pub fn process(&self, mut s: Spectrum) -> ProcessedSpectrum {
        let charge = self.max_fragment_charge.min(s.charge);

        s.peaks.sort_unstable_by(|(_, a), (_, b)| b.total_cmp(&a));
        let n = s.peaks.len().min(self.take_top_n);
        let mut peaks = s.peaks[..n]
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
                    match fragment_mass <= self.max_fragment_mz {
                        true => Some((fragment_mass, int.sqrt())),
                        false => None,
                    }
                })
            })
            .collect::<Vec<_>>();

        // Sort by m/z
        peaks.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));

        ProcessedSpectrum {
            scan: s.scan,
            monoisotopic_mass: s.precursor_mass - PROTON,
            charge: s.charge,
            rt: s.rt,
            peaks,
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
    pub peaks: Vec<(f32, f32)>,
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
pub fn read_ms2<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<Vec<Spectrum>> {
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
            current.peaks.push((mz, abundance));
        }
    }
    spectra.push(current);

    Ok(spectra)
}
