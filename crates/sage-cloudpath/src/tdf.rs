use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use serde::{Deserialize, Serialize};
use std::{cmp::Ordering, path::Path};
use timsrust::core::{Converter, ScanIndex};
use timsrust::tdf::{self, Scan2ImConverter, Tof2MzConverter};
use timsrust::SpectrumReader;

// Serializable settings preserve the existing Sage parameter format.
#[derive(Default, Deserialize, Serialize, Debug, Clone, Copy)]
pub struct TimsrustSpectrumConfig {
    pub spectrum_processing_params: SpectrumProcessingParams,
    pub frame_splitting_params: FrameWindowSplittingConfiguration,
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
pub struct SpectrumProcessingParams {
    pub smoothing_window: u32,
    pub centroiding_window: u32,
    pub calibration_tolerance: f64,
    pub calibrate: bool,
}

impl Default for SpectrumProcessingParams {
    fn default() -> Self {
        let params = tdf::SpectrumProcessingParams::default();
        Self {
            smoothing_window: params.smoothing_window,
            centroiding_window: params.centroiding_window,
            calibration_tolerance: params.calibration_tolerance,
            calibrate: params.calibrate,
        }
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
pub enum FrameWindowSplittingConfiguration {
    Quadrupole(QuadWindowExpansionStrategy),
    Window(QuadWindowExpansionStrategy),
}

impl Default for FrameWindowSplittingConfiguration {
    fn default() -> Self {
        Self::Quadrupole(QuadWindowExpansionStrategy::Even(1))
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
pub enum QuadWindowExpansionStrategy {
    None,
    Even(usize),
    // The legacy converter is replaced by the dataset's calibration at runtime.
    UniformMobility(
        (f64, f64),
        #[serde(deserialize_with = "ignore_legacy_converter")] Option<()>,
    ),
    UniformScan((usize, usize)),
}

impl From<TimsrustSpectrumConfig> for tdf::SpectrumReaderConfig<timsrust::ImConverter> {
    fn from(config: TimsrustSpectrumConfig) -> Self {
        let params = config.spectrum_processing_params;
        let convert = |strategy| match strategy {
            QuadWindowExpansionStrategy::None => tdf::QuadWindowExpansionStrategy::None,
            QuadWindowExpansionStrategy::Even(n) => tdf::QuadWindowExpansionStrategy::Even(n),
            QuadWindowExpansionStrategy::UniformMobility(span_step, _) => {
                tdf::QuadWindowExpansionStrategy::UniformMobility(span_step, None)
            }
            QuadWindowExpansionStrategy::UniformScan(span_step) => {
                tdf::QuadWindowExpansionStrategy::UniformScan(span_step)
            }
        };
        Self {
            spectrum_processing_params: tdf::SpectrumProcessingParams {
                smoothing_window: params.smoothing_window,
                centroiding_window: params.centroiding_window,
                calibration_tolerance: params.calibration_tolerance,
                calibrate: params.calibrate,
            },
            frame_splitting_params: match config.frame_splitting_params {
                FrameWindowSplittingConfiguration::Quadrupole(strategy) => {
                    tdf::FrameWindowSplittingConfiguration::Quadrupole(convert(strategy))
                }
                FrameWindowSplittingConfiguration::Window(strategy) => {
                    tdf::FrameWindowSplittingConfiguration::Window(convert(strategy))
                }
            },
        }
    }
}

fn ignore_legacy_converter<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Option<()>, D::Error> {
    serde::de::IgnoredAny::deserialize(deserializer)?;
    Ok(None)
}

pub struct TdfReader;

#[derive(Deserialize, Serialize, Debug, Clone, Copy)]
pub struct BrukerMS1CentoidingConfig {
    pub mz_ppm: f32,
    pub ims_pct: f32,
}

impl Default for BrukerMS1CentoidingConfig {
    fn default() -> Self {
        BrukerMS1CentoidingConfig {
            mz_ppm: 5.0,
            ims_pct: 3.0,
        }
    }
}

#[derive(Default, Deserialize, Serialize, Debug, Clone, Copy)]
pub struct BrukerProcessingConfig {
    pub ms2: TimsrustSpectrumConfig,
    pub ms1: BrukerMS1CentoidingConfig,
}

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<Path>,
        file_id: usize,
        config: BrukerProcessingConfig,
        requires_ms1: bool,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let path = timsrust::TimsTofPath::new(path_name.as_ref().to_string_lossy().as_ref())?;
        let spectrum_reader = SpectrumReader::build()
            .with_path(path)
            .with_config(config.ms2.into())
            .finalize()?;
        let mut spectra = self.read_msn_spectra(file_id, &spectrum_reader)?;
        if requires_ms1 {
            let ms1s = self.read_ms1_spectra(&path_name, file_id, config.ms1)?;
            spectra.extend(ms1s);
        }

        Ok(spectra)
    }

    fn read_ms1_spectra(
        &self,
        path_name: impl AsRef<Path>,
        file_id: usize,
        config: BrukerMS1CentoidingConfig,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let start = std::time::Instant::now();
        let path = path_name.as_ref().to_string_lossy();
        let frame_reader = timsrust::tdf::TdfFrameReader::new(path.as_ref())?;
        let mz_converter = Tof2MzConverter::new(path.as_ref());
        let ims_converter = Scan2ImConverter::new(path.as_ref());
        let tol_ppm = config.mz_ppm;
        let im_tol_pct = config.ims_pct;

        let ms1_spectra: Vec<RawSpectrum> = frame_reader
            .parallel_filter(|f| f.info().ms_level() == timsrust::core::MSLevel::MS1)
            .map_init(
                || PeakBuffer::with_capacity(2 * MAX_PEAKS),
                |buffer, frame| match frame {
                    Ok(frame) => {
                        buffer.clear();
                        buffer.with_frame(&frame, &ims_converter, &mz_converter);

                        // Squash the mobility dimension
                        let (mz, (intensity, mobility)): (Vec<f32>, (Vec<f32>, Vec<f32>)) =
                            buffer.fastcentroid_frame(tol_ppm, im_tol_pct);

                        let scan_start_time = frame.info().rt_in_seconds() as f32 / 60.0;
                        let ion_injection_time = 100.0; // This is made up, in theory we can read
                                                        // if from the tdf file
                        let total_ion_current = intensity.iter().sum::<f32>();
                        let id = frame.index().to_string();

                        let spec = RawSpectrum {
                            file_id,
                            precursors: vec![],
                            representation: Representation::Centroid,
                            scan_start_time,
                            ion_injection_time,
                            mz,
                            ms_level: 1,
                            id,
                            intensity,
                            total_ion_current,
                            mobility: Some(mobility),
                        };
                        Some(spec)
                    }
                    Err(x) => {
                        log::error!("error parsing spectrum: {:?}", x);
                        None
                    }
                },
            )
            .flatten()
            .collect();
        log::info!(
            "read {} ms1 spectra in {:#?}",
            ms1_spectra.len(),
            start.elapsed()
        );
        Ok(ms1_spectra)
    }

    fn read_msn_spectra(
        &self,
        file_id: usize,
        spectrum_reader: &SpectrumReader,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        Ok(spectrum_reader
            .par_iter()
            .filter_map(|spectrum| spectrum.ok())
            .filter_map(|spectrum| Self::parse_spectrum(file_id, &spectrum))
            .collect())
    }

    fn parse_spectrum(
        file_id: usize,
        spectrum: &timsrust::core::Spectrum<timsrust::core::Mz>,
    ) -> Option<RawSpectrum> {
        let dda_precursor = spectrum.precursor().as_ref()?;
        let mut precursor = Self::parse_precursor(dda_precursor);
        let window = spectrum.isolation_window();
        let mz = f64::from(dda_precursor.mz());
        // timsrust 0.6 uses zero when MonoisotopicMz is NULL in the TDF.
        // Such precursors still have a measured quadrupole isolation center.
        let mz = if mz.is_finite() && mz > 0.0 {
            mz
        } else {
            f64::from(window.center())
        };
        precursor.mz = mz as f32;
        precursor.isolation_window = Some(Tolerance::Da(
            (f64::from(window.lower()) - mz) as f32,
            (f64::from(window.upper()) - mz) as f32,
        ));
        Some(RawSpectrum {
            file_id,
            precursors: vec![precursor],
            representation: Representation::Centroid,
            scan_start_time: f64::from(dda_precursor.rt()) as f32 / 60.0,
            ion_injection_time: f64::from(dda_precursor.rt()) as f32,
            total_ion_current: 0.0,
            mz: spectrum
                .mz_values()
                .iter()
                .map(|&x| f64::from(x) as f32)
                .collect(),
            ms_level: 2,
            id: spectrum.index().to_string(),
            intensity: spectrum.intensities().iter().map(|&x| x as f32).collect(),
            mobility: None,
        })
    }

    fn parse_precursor(dda_precursor: &timsrust::core::Precursor) -> Precursor {
        let mut precursor: Precursor = Precursor::default();
        precursor.mz = f64::from(dda_precursor.mz()) as f32;
        precursor.charge = dda_precursor
            .charge()
            .and_then(|x| u8::try_from(i8::from(x)).ok());
        precursor.intensity = dda_precursor.intensity().map(|x| x as f32);
        precursor.spectrum_ref = Option::from(dda_precursor.frame_index().to_string());
        precursor.inverse_ion_mobility = Option::from(f64::from(dda_precursor.im()) as f32);
        precursor
    }
}

#[derive(Clone, Copy)]
struct ImsPeak {
    mz: f32,
    intensity: f32,
    im: f32,
}
const MAX_PEAKS: usize = 10_000;

/// Buffer that gets re-used on each thread to store the intermediates
/// of the centroiding for a single frame.
#[derive(Clone)]
struct PeakBuffer {
    peaks: Vec<ImsPeak>,
    order: Vec<usize>,
    agg_buff: Vec<ImsPeak>,
}

impl PeakBuffer {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            peaks: Vec::with_capacity(capacity),
            order: Vec::with_capacity(capacity),
            agg_buff: Vec::with_capacity(MAX_PEAKS),
        }
    }

    fn with_frame(
        &mut self,
        frame: &timsrust::core::Frame,
        ims_converter: &Scan2ImConverter,
        mz_converter: &Tof2MzConverter,
    ) {
        let expect_len = frame.ions().tof_indices().len();
        self.expand_to_capacity(expect_len);

        let mz_iter = frame
            .ions()
            .tof_indices()
            .iter()
            .map(|&x| f64::from(mz_converter.convert(x)) as f32);
        let intensities_iter = frame
            .ions()
            .intensities()
            .iter()
            .map(|&x| u32::from(x) as f32);
        let imss_iter = Self::expand_mobility_iter(frame.ions().scan_offsets(), ims_converter);

        let peak_iter = mz_iter
            .zip(intensities_iter)
            .zip(imss_iter)
            .map(|((mz, intensity), im)| ImsPeak { mz, intensity, im });
        self.peaks.extend(peak_iter);
        assert_eq!(self.peaks.len(), expect_len);

        // sort by mz ... bc binary searching on the mz space
        // for neighbors is the fastest way to find neighbors that I have tried.
        self.peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

        // The "order" is sorted by intensity
        // This will be used later during the centroiding (for details check that implementation)
        self.order.extend(0..self.len());
        self.order.sort_unstable_by(|&a, &b| {
            self.peaks[b]
                .intensity
                .partial_cmp(&self.peaks[a].intensity)
                .unwrap_or(Ordering::Equal)
        });
    }

    fn clear(&mut self) {
        self.peaks.clear();
        self.order.clear();
        self.agg_buff.clear();
    }

    fn expand_to_capacity(&mut self, capacity: usize) {
        if capacity <= self.len() {
            return;
        }
        let diff = capacity - self.len();
        // Grow by whatever is the largest 20% of the current capacity
        // or the difference.
        let diff = diff.max(self.len() / 5);

        self.peaks.reserve(diff);
        self.order.reserve(diff);
        self.agg_buff.reserve(capacity);
    }

    fn len(&self) -> usize {
        self.peaks.len()
    }

    /// Expand the scan offset slice to mobilities.
    ///
    /// The scan offsets is in essence a run-length
    /// encoded vector of scan numbers that can be converter to the 1/k0
    /// values.
    ///
    /// Essentially ... the slice [0,4,5,5], would expand to
    /// [0,0,0,0,1]; 0 to 4 have index 0, 4 to 5 have index 1, 5 to 5 would
    /// have index 2 but its empty!
    ///
    /// Then this index can be converted using the Scan2ImConverter.convert
    ///
    /// ... This should problably be implemented and exposed in timsrust.
    fn expand_mobility_iter<'a>(
        scan_offsets: &'a [usize],
        ims_converter: &'a Scan2ImConverter,
    ) -> impl Iterator<Item = f32> + 'a {
        let ims_iter = scan_offsets
            .windows(2)
            .enumerate()
            .filter_map(|(i, w)| {
                let num = w[1] - w[0];
                if num == 0 {
                    return None;
                }
                let lo = w[0];
                let hi = w[1];

                let im = f64::from(
                    ims_converter.convert(ScanIndex::try_from(i).expect("scan index out of range")),
                ) as f32;
                Some((im, lo, hi))
            })
            .flat_map(|(im, lo, hi)| (lo..hi).map(move |_| im));
        ims_iter
    }

    /// Centroiding of the IM-containing spectra
    ///
    /// This is a very rudimentary centroiding algorithm but... it seems to work well.
    /// It iterativelty goes over the peaks in decreasing intensity order and
    /// accumulates the intensity of the peaks surrounding the peak. (sort of
    /// like the first pass in dbscan).
    ///
    /// The preserved mobility and mz are the ones from the apex peak.
    /// A more complex version where the weighted mean is preserved is possible
    /// but I have seen only marginal gains and a lot more complexity + time.
    ///
    /// This dramatically reduces the number of peaks in the spectra
    /// which saves a ton of memory and time when doing LFQ, since we
    /// iterate over each peak.
    fn fastcentroid_frame(
        &mut self,
        mz_tol_ppm: f32,
        im_tol_pct: f32,
    ) -> (Vec<f32>, (Vec<f32>, Vec<f32>)) {
        // Make sure the array is mz sorted ... I should delete
        // this assertions once I am confident of the implementation.
        // but tbh, its not that slow and its simple.
        assert!(
            self.peaks.windows(2).all(|x| x[0].mz <= x[1].mz),
            "mz_array is not sorted"
        );
        assert!(self.agg_buff.is_empty(), "agg_buff is not empty");

        let mut global_num_included = 0;

        let utol = mz_tol_ppm / 1e6;
        let im_tol = im_tol_pct / 100.0;

        for &idx in &self.order {
            if self.peaks[idx].intensity <= 0.0 {
                continue;
            }
            if self.agg_buff.len() > MAX_PEAKS {
                let curr_loc_int = self.peaks[idx].intensity;
                if curr_loc_int > 200.0 {
                    log::debug!(
                        "Reached limit of the agg buffer at index {}/{} curr int={}",
                        idx,
                        self.len(),
                        curr_loc_int
                    );
                }
                break;
            }

            let mz = self.peaks[idx].mz;
            let im = self.peaks[idx].im;
            let da_tol = mz * utol;
            let left_e = mz - da_tol;
            let right_e = mz + da_tol;

            let ss_start = self.peaks.partition_point(|&x| x.mz < left_e);
            let ss_end = self.peaks.partition_point(|&x| x.mz <= right_e);

            let abs_im_tol = im * im_tol;
            let left_im = im - abs_im_tol;
            let right_im = im + abs_im_tol;

            let mut curr_intensity = 0.0;

            let mut num_includable = 0;
            for i in ss_start..ss_end {
                let im_i = self.peaks[i].im;
                if (self.peaks[i].intensity > 0.0) && im_i >= left_im && im_i <= right_im {
                    curr_intensity += self.peaks[i].intensity;
                    self.peaks[i].intensity = -1.0;
                    num_includable += 1;
                }
            }

            assert!(num_includable > 0, "At least 'itself' should be included");

            self.agg_buff.push(ImsPeak {
                mz,
                intensity: curr_intensity,
                im,
            });
            global_num_included += num_includable;

            if global_num_included == self.len() {
                log::debug!("All peaks were included in the centroiding");
                break;
            }
        }

        self.agg_buff
            .sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
        // println!("Centroiding: Start len: {}; end len: {};", arr_len, result.len());
        // Ultra data is usually start: 40k end 10k,
        // HT2 data is usually start 400k end 40k, limiting to 10k
        // rarely leaves peaks with intensity > 200 ... ive never seen
        // it happen. -JSP 2025-Jan

        self.agg_buff
            .drain(..)
            .map(|x| (x.mz, (x.intensity, x.im)))
            .unzip()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use timsrust::core::{Charge, FrameIndex, Im, IsolationWindow, Mz, Rt, Spectrum};

    #[test]
    fn missing_monoisotopic_mz_uses_isolation_center() {
        let precursor = timsrust::core::Precursor::new(
            Mz::from(0.0),
            Im::from(1.2),
            Rt::from(120.0),
            ScanIndex::try_from(0usize).unwrap(),
            None,
            Some(1000.0),
            1,
            FrameIndex::try_from(7usize).unwrap(),
        );
        let spectrum = Spectrum::new(
            vec![10.0],
            1,
            Some(precursor),
            vec![Mz::from(100.0)],
            IsolationWindow::new_from_center(Mz::from(924.0), Mz::from(3.0), 30.0),
        );
        let raw = TdfReader::parse_spectrum(0, &spectrum).unwrap();
        assert_eq!(raw.precursors[0].mz, 924.0);
        assert_eq!(raw.precursors[0].charge, None);
        assert!(matches!(
            raw.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 1.5))
        ));
    }

    #[test]
    fn skips_spectra_without_precursors() {
        let spectrum = Spectrum::new(vec![], 0, None, vec![], IsolationWindow::default());
        assert!(TdfReader::parse_spectrum(0, &spectrum).is_none());
    }

    #[test]
    fn legacy_processing_config_roundtrips() {
        use timsrust::tdf::{
            FrameWindowSplittingConfiguration as F, QuadWindowExpansionStrategy as Q,
        };
        for container in ["Quadrupole", "Window"] {
            for strategy in [
                serde_json::json!("None"),
                serde_json::json!({"Even": 3}),
                serde_json::json!({"UniformScan": [20, 10]}),
                serde_json::json!({"UniformMobility": [[0.1, 0.05], null]}),
            ] {
                let mut json = serde_json::to_value(BrukerProcessingConfig::default()).unwrap();
                json["ms2"]["spectrum_processing_params"] = serde_json::json!({
                    "smoothing_window": 3, "centroiding_window": 5,
                    "calibration_tolerance": 0.2, "calibrate": true,
                });
                json["ms2"]["frame_splitting_params"] = serde_json::json!({container: strategy});
                let config: BrukerProcessingConfig = serde_json::from_value(json.clone()).unwrap();
                assert_eq!(serde_json::to_value(config).unwrap(), json);
                let mapped: timsrust::tdf::SpectrumReaderConfig<timsrust::ImConverter> =
                    config.ms2.into();
                assert_eq!(mapped.spectrum_processing_params.smoothing_window, 3);
                assert_eq!(mapped.spectrum_processing_params.centroiding_window, 5);
                assert_eq!(mapped.spectrum_processing_params.calibration_tolerance, 0.2);
                assert!(mapped.spectrum_processing_params.calibrate);
                let mapped_strategy = match (container, mapped.frame_splitting_params) {
                    ("Quadrupole", F::Quadrupole(s)) | ("Window", F::Window(s)) => s,
                    _ => panic!("changed splitting mode"),
                };
                match mapped_strategy {
                    Q::None => assert_eq!(strategy, serde_json::json!("None")),
                    Q::Even(n) => assert_eq!(n, 3),
                    Q::UniformScan(pair) => assert_eq!(pair, (20, 10)),
                    Q::UniformMobility(pair, converter) => {
                        assert_eq!(pair, (0.1, 0.05));
                        assert!(converter.is_none());
                    }
                }
            }
        }
    }
}
