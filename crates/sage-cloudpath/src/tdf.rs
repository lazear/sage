use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{
        BrukerSpectrumProcessor, FrameWindowSplittingStrategy, QuadWindowExpansionStrategy,
    },
    spectrum::{Precursor, RawSpectrum, Representation},
};

pub struct TdfReader;

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<str>,
        file_id: usize,
        bruker_spectrum_processor: BrukerSpectrumProcessor,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let spectrum_reader = timsrust::readers::SpectrumReader::build()
            .with_path(path_name.as_ref())
            .with_config(Self::parse_bruker_spectrum_config(
                bruker_spectrum_processor,
            ))
            .finalize()?;
        let spectra: Vec<RawSpectrum> = (0..spectrum_reader.len())
            .into_par_iter()
            .filter_map(|index| match spectrum_reader.get(index) {
                Ok(dda_spectrum) => match dda_spectrum.precursor {
                    Some(dda_precursor) => {
                        let mut precursor = Self::parse_precursor(dda_precursor);
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
                            intensity: dda_spectrum.intensities.iter().map(|&x| x as f32).collect(),
                        };
                        Some(spectrum)
                    }
                    None => None,
                },
                Err(_) => None,
            })
            .collect();
        Ok(spectra)
    }

    fn parse_precursor(dda_precursor: timsrust::Precursor) -> Precursor {
        let mut precursor: Precursor = Precursor::default();
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
        precursor
    }

    fn parse_bruker_spectrum_config(
        config: BrukerSpectrumProcessor,
    ) -> timsrust::readers::SpectrumReaderConfig {
        let mut spectrum_processing_params = timsrust::readers::SpectrumProcessingParams::default();
        if let Some(smoothing_window) = config.smoothing_window {
            spectrum_processing_params.smoothing_window = smoothing_window;
        }
        if let Some(centroiding_window) = config.centroiding_window {
            spectrum_processing_params.centroiding_window = centroiding_window;
        }
        if let Some(calibration_tolerance) = config.calibration_tolerance {
            spectrum_processing_params.calibration_tolerance = calibration_tolerance;
        }
        if let Some(calibrate) = config.calibrate {
            spectrum_processing_params.calibrate = calibrate;
        }
        let frame_splitting_params: timsrust::readers::FrameWindowSplittingConfiguration;
        if let Some(dia_strategy) = config.dia_strategy {
            frame_splitting_params = match dia_strategy {
                FrameWindowSplittingStrategy::Quadrupole(q) => {
                    timsrust::readers::FrameWindowSplittingConfiguration::Quadrupole(match q {
                        QuadWindowExpansionStrategy::None => {
                            timsrust::readers::QuadWindowExpansionStrategy::None
                        }
                        QuadWindowExpansionStrategy::Even(n) => {
                            timsrust::readers::QuadWindowExpansionStrategy::Even(n)
                        }
                        QuadWindowExpansionStrategy::UniformScan(x) => {
                            timsrust::readers::QuadWindowExpansionStrategy::UniformScan(x)
                        }
                        QuadWindowExpansionStrategy::UniformMobility(x) => {
                            timsrust::readers::QuadWindowExpansionStrategy::UniformMobility(x, None)
                        }
                    })
                }
                FrameWindowSplittingStrategy::Window(q) => {
                    timsrust::readers::FrameWindowSplittingConfiguration::Window(match q {
                        QuadWindowExpansionStrategy::None => {
                            timsrust::readers::QuadWindowExpansionStrategy::None
                        }

                        QuadWindowExpansionStrategy::Even(n) => {
                            timsrust::readers::QuadWindowExpansionStrategy::Even(n)
                        }

                        QuadWindowExpansionStrategy::UniformScan(x) => {
                            timsrust::readers::QuadWindowExpansionStrategy::UniformScan(x)
                        }
                        QuadWindowExpansionStrategy::UniformMobility(x) => {
                            timsrust::readers::QuadWindowExpansionStrategy::UniformMobility(x, None)
                        }
                    })
                }
            };
        } else {
            frame_splitting_params = timsrust::readers::FrameWindowSplittingConfiguration::default()
        }
        timsrust::readers::SpectrumReaderConfig {
            spectrum_processing_params,
            frame_splitting_params,
        }
    }
}
