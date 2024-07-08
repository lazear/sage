use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};

pub struct TdfReader;

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<str>,
        file_id: usize,
    ) -> Result<Vec<RawSpectrum>, timsrust::Error> {
        let spectrum_reader = timsrust::io::readers::SpectrumReader::new(path_name.as_ref());
        let spectra: Vec<RawSpectrum> = (0..spectrum_reader.len())
            .into_par_iter()
            .map(|index| {
                let dda_spectrum = spectrum_reader.get(index);
                let mut precursor: Precursor = Precursor::default();
                let dda_precursor: timsrust::ms_data::Precursor = dda_spectrum.precursor;
                precursor.mz = dda_precursor.mz as f32;
                precursor.charge = match dda_precursor.charge {
                    Some(x) => Some(x as u8),
                    None => None,
                };
                precursor.intensity = match dda_precursor.intensity {
                    Some(x) => Some(x as f32),
                    None => None,
                };
                precursor.spectrum_ref = Option::from(dda_precursor.frame_index.to_string());
                precursor.inverse_ion_mobility = Option::from(dda_precursor.im as f32);
                precursor.isolation_window = Option::from(Tolerance::Da(
                    -dda_spectrum.isolation_width as f32 / 2.0,
                    dda_spectrum.isolation_width as f32 / 2.0,
                ));
                let spectrum: RawSpectrum = RawSpectrum {
                    file_id,
                    precursors: vec![precursor],
                    representation: Representation::Centroid,
                    scan_start_time: dda_precursor.rt as f32 / 60.0,
                    ion_injection_time: dda_precursor.rt as f32,
                    total_ion_current: 0.0,
                    mz: dda_spectrum.mz_values.iter().map(|&x| x as f32).collect(),
                    ms_level: 2,
                    id: dda_spectrum.index.to_string(),
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
