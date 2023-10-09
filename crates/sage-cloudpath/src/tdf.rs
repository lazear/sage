use rayon::prelude::*;
use sage_core::spectrum::{Precursor, Representation, RawSpectrum};
use std::fmt::{Display, Formatter};
use timsrust;

#[derive(Default)]
pub struct TdfReader {}

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<str>,
        file_id: usize,
    ) -> Result<Vec<RawSpectrum>, TdfError> {
        let dda_spectra: Vec<timsrust::Spectrum> =
            timsrust::FileReader::new(path_name.as_ref().to_string()).read_all_spectra();
        let spectra: Vec<RawSpectrum> = (0..dda_spectra.len())
            .into_par_iter()
            .map(|index| {
                let dda_spectrum = &dda_spectra[index];
                let mut precursor: Precursor = Precursor::default();
                let dda_precursor: timsrust::Precursor = dda_spectrum.precursor.unwrap_as_precursor();
                precursor.mz = dda_precursor.mz as f32;
                precursor.charge = Option::from(dda_precursor.charge as u8);
                // precursor.ion_mobility = Option::from(dda_precursor.im as f32);
                precursor.intensity = Option::from(dda_precursor.intensity as f32);
                precursor.spectrum_ref = Option::from(dda_precursor.frame_index.to_string());
                let spectrum: RawSpectrum = RawSpectrum {
                    file_id: file_id,
                    precursors: vec![precursor],
                    representation: Representation::Centroid,
                    scan_start_time: dda_precursor.rt as f32,
                    ion_injection_time: dda_precursor.rt as f32,
                    total_ion_current: 0.0,
                    mz: dda_spectrum.mz_values.iter().map(|&x| x as f32).collect(),
                    ms_level: 2,
                    id: dda_precursor.index.to_string(),
                    // precursor_id: dda_precursor.index as u32,
                    // frame_id: dda_precursor.frame_index as u32,
                    intensity: dda_spectrum.intensities.iter().map(|&x| x as f32).collect(),
                };
                spectrum
            })
            .collect();
        Ok(spectra)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum TdfError {
    Unreadable,
}

impl Display for TdfError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            TdfError::Unreadable => f.write_str("Tdf Error : Malformed data."),
        }
    }
}
