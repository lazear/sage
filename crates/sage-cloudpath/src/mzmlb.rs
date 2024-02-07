#![cfg(feature = "mzdata")]

use std::io;
use std::path::Path;

use mzdata::io::MzMLbReader as MzMLbReaderImpl;
use mzdata::prelude::*;
use mzdata::RawSpectrum as RawSpectrumImpl;
use sage_core::mass::Tolerance;
use sage_core::spectrum::{Precursor, RawSpectrum};

pub struct MzMLbReader {
    ms_level: Option<u8>,
    // If set to Some(level) and noise intensities are present in the MzML file,
    // divide intensities at this MS-level by noise to calculate S/N
    signal_to_noise: Option<u8>,

    file_id: usize,
}

impl MzMLbReader {
    /// Create a new [`MzMLbReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_file_id_and_level_filter(file_id: usize, ms_level: u8) -> Self {
        Self {
            ms_level: Some(ms_level),
            file_id,
            signal_to_noise: None,
        }
    }

    pub fn with_file_id(file_id: usize) -> Self {
        Self {
            ms_level: None,
            signal_to_noise: None,
            file_id,
        }
    }

    pub fn set_file_id(&mut self, file_id: usize) -> &mut Self {
        self.file_id = file_id;
        self
    }

    pub fn set_signal_to_noise(&mut self, sn: Option<u8>) -> &mut Self {
        self.signal_to_noise = sn;
        self
    }

    pub fn parse<B>(&self, b: B) -> Result<Vec<RawSpectrum>, io::Error>
    where
        B: AsRef<Path>,
    {
        let reader = MzMLbReaderImpl::new(&b)?;

        let spectra = reader
            .into_iter()
            .filter(|scan| {
                if let Some(ms_level) = self.ms_level {
                    return scan.ms_level() == ms_level;
                } else {
                    return true;
                }
            })
            .map(|scan| {
                let scan: RawSpectrumImpl = scan.into();
                let mut precusors = Vec::new();
                match scan.precursor() {
                    Some(p) => {
                        let p = Precursor {
                            mz: p.mz() as f32,
                            intensity: Some(p.ion.intensity),
                            charge: p.ion.charge.and_then(|v| Some(v as u8)),
                            spectrum_ref: p.precursor_id.clone(),
                            isolation_window: Some(Tolerance::Da(
                                p.isolation_window.lower_bound as f32,
                                p.isolation_window.upper_bound as f32,
                            )),
                        };
                        precusors.push(p)
                    }
                    None => {}
                }
                RawSpectrum {
                    file_id: self.file_id,
                    ms_level: scan.ms_level(),
                    id: scan.description.id.clone(),
                    precursors: precusors,
                    ion_injection_time: scan.acquisition().first_scan().unwrap().injection_time,
                    representation: match scan.description.signal_continuity {
                        mzdata::spectrum::SignalContinuity::Unknown => {
                            sage_core::spectrum::Representation::Profile
                        }
                        mzdata::spectrum::SignalContinuity::Centroid => {
                            sage_core::spectrum::Representation::Centroid
                        }
                        mzdata::spectrum::SignalContinuity::Profile => {
                            sage_core::spectrum::Representation::Profile
                        }
                    },
                    scan_start_time: scan.start_time() as f32,
                    total_ion_current: scan.peaks().tic(),
                    mz: scan.mzs().iter().map(|mz| (*mz) as f32).collect(),
                    intensity: scan.intensities().to_vec(),
                }
            })
            .collect();
        Ok(spectra)
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_read_mzmlb() -> io::Result<()> {
        let spectra = MzMLbReader::with_file_id(0).parse("../../tests/LQSRPAAPPAPGPGQLTLR.mzMLb")?;
        assert_eq!(spectra.len(), 1);
        let s = spectra.first().unwrap();
        assert_eq!(s.id, "controllerType=0 controllerNumber=1 scan=30069");
        assert_eq!(s.mz.len(), 299);
        assert_eq!(s.intensity.len(), 299);
        Ok(())
    }

}