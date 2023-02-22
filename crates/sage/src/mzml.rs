use crate::mass::Tolerance;
use crate::spectrum::Precursor;
use async_compression::tokio::bufread::ZlibDecoder;
use quick_xml::events::Event;
use quick_xml::Reader;
use tokio::io::{AsyncBufRead, AsyncReadExt};

#[derive(Default, Debug, Clone)]
pub struct Spectrum {
    pub ms_level: u8,
    pub id: String,
    // pub scan_id: Option<usize>,
    pub precursors: Vec<Precursor>,
    /// Profile or Centroided data
    pub representation: Representation,
    /// Scan start time
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Total ion current
    pub total_ion_current: f32,
    /// M/z array
    pub mz: Vec<f32>,
    /// Intensity array
    pub intensity: Vec<f32>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Representation {
    Profile,
    Centroid,
}

impl Default for Representation {
    fn default() -> Self {
        Self::Profile
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
/// Which tag are we inside?
enum State {
    Spectrum,
    Scan,
    BinaryDataArray,
    Binary,
    Precursor,
    SelectedIon,
}

#[derive(Copy, Clone, Debug)]
enum BinaryKind {
    Intensity,
    Mz,
}

#[derive(Copy, Clone, Debug)]
enum Dtype {
    F32,
    F64,
}

// MUST supply only one of the following
const ZLIB_COMPRESSION: &str = "MS:1000574";
const NO_COMPRESSION: &str = "MS:1000576";

// MUST supply only one of the following
const INTENSITY_ARRAY: &str = "MS:1000515";
const MZ_ARRAY: &str = "MS:1000514";

// MUST supply only one of the following
const FLOAT_64: &str = "MS:1000523";
const FLOAT_32: &str = "MS:1000521";

const MS_LEVEL: &str = "MS:1000511";
const PROFILE: &str = "MS:1000128";
const CENTROID: &str = "MS:1000127";
const TOTAL_ION_CURRENT: &str = "MS:1000285";

const SCAN_START_TIME: &str = "MS:1000016";
const ION_INJECTION_TIME: &str = "MS:1000927";

const SELECTED_ION_MZ: &str = "MS:1000744";
const SELECTED_ION_INT: &str = "MS:1000042";
const SELECTED_ION_CHARGE: &str = "MS:1000041";

const ISO_WINDOW_LOWER: &str = "MS:1000828";
const ISO_WINDOW_UPPER: &str = "MS:1000829";

#[derive(Default)]
pub struct MzMLReader {
    ms_level: Option<u8>,
}

impl MzMLReader {
    /// Create a new [`MzMlReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_level_filter(ms_level: u8) -> Self {
        Self {
            ms_level: Some(ms_level),
        }
    }

    /// Here be dragons -
    /// Seriously, this kinda sucks because it's a giant imperative, stateful loop.
    /// But I also don't want to spend any more time working on an mzML parser...
    pub async fn parse<B: AsyncBufRead + Unpin>(&self, b: B) -> Result<Vec<Spectrum>, MzMLError> {
        let mut reader = Reader::from_reader(b);
        let mut buf = Vec::new();

        let mut state = None;
        let mut compression = false;
        let mut output_buffer = Vec::with_capacity(4096);
        let mut binary_dtype = Dtype::F64;
        let mut binary_array = None;

        let mut spectrum = Spectrum::default();
        let mut precursor = Precursor::default();
        let mut iso_window_lo: Option<f32> = None;
        let mut iso_window_hi: Option<f32> = None;
        let mut spectra = Vec::new();

        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)?
                    .ok_or_else(|| MzMLError::Malformed)?
                    .value
            };
        }

        loop {
            match reader.read_event_into_async(&mut buf).await {
                Ok(Event::Start(ref ev)) => {
                    // State transition into child tag
                    state = match (ev.name().into_inner(), state) {
                        (b"spectrum", _) => Some(State::Spectrum),
                        (b"scan", Some(State::Spectrum)) => Some(State::Scan),
                        (b"binaryDataArray", Some(State::Spectrum)) => Some(State::BinaryDataArray),
                        (b"binary", Some(State::BinaryDataArray)) => Some(State::Binary),
                        (b"precursor", Some(State::Spectrum)) => Some(State::Precursor),
                        (b"selectedIon", Some(State::Precursor)) => Some(State::SelectedIon),
                        _ => state,
                    };
                    match ev.name().into_inner() {
                        b"spectrum" => {
                            let id = extract!(ev, b"id");
                            let id = std::str::from_utf8(&id)?;
                            spectrum.id = id.to_string();
                        }
                        b"precursor" => {
                            // Not all precursor fields have a spectrumRef
                            if let Some(scan) = ev.try_get_attribute(b"spectrumRef")? {
                                let scan = std::str::from_utf8(&scan.value)?;
                                precursor.spectrum_ref = Some(scan.to_string())
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::Empty(ref ev)) => match (state, ev.name().into_inner()) {
                    (Some(State::BinaryDataArray), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        match accession {
                            ZLIB_COMPRESSION => compression = true,
                            NO_COMPRESSION => compression = false,
                            FLOAT_64 => binary_dtype = Dtype::F64,
                            FLOAT_32 => binary_dtype = Dtype::F32,
                            INTENSITY_ARRAY => binary_array = Some(BinaryKind::Intensity),
                            MZ_ARRAY => binary_array = Some(BinaryKind::Mz),
                            _ => {
                                // Unknown CV - perhaps noise
                                binary_array = None;
                            }
                        }
                    }
                    (Some(State::Spectrum), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        match accession {
                            MS_LEVEL => {
                                let level = extract!(ev, b"value");
                                let level = std::str::from_utf8(&level)?.parse::<u8>()?;
                                if let Some(filter) = self.ms_level {
                                    if level != filter {
                                        spectrum = Spectrum::default();
                                        state = None;
                                    }
                                }
                                spectrum.ms_level = level;
                            }
                            PROFILE => spectrum.representation = Representation::Profile,
                            CENTROID => spectrum.representation = Representation::Centroid,
                            TOTAL_ION_CURRENT => {
                                let value = extract!(ev, b"value");
                                let value = std::str::from_utf8(&value)?.parse::<f32>()?;
                                if value == 0.0 {
                                    // No ion current, break out of current state
                                    spectrum = Spectrum::default();
                                    state = None;
                                } else {
                                    spectrum.total_ion_current = value;
                                }
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Precursor), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value)?;
                        match accession {
                            ISO_WINDOW_LOWER => {
                                iso_window_lo = Some(value.parse()?);
                            }
                            ISO_WINDOW_UPPER => {
                                iso_window_hi = Some(value.parse()?);
                            }
                            _ => {}
                        }
                    }
                    (Some(State::SelectedIon), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value)?;
                        match accession {
                            SELECTED_ION_CHARGE => {
                                precursor.charge = Some(value.parse()?);
                            }
                            SELECTED_ION_MZ => {
                                precursor.mz = value.parse()?;
                            }
                            SELECTED_ION_INT => {
                                precursor.intensity = Some(value.parse()?);
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Scan), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value)?;
                        match accession {
                            SCAN_START_TIME => {
                                spectrum.scan_start_time = value.parse()?;
                            }
                            ION_INJECTION_TIME => {
                                spectrum.ion_injection_time = value.parse()?;
                            }
                            _ => {}
                        }
                    }

                    _ => {}
                },
                Ok(Event::Text(text)) => {
                    if let Some(State::Binary) = state {
                        if let Some(filter) = self.ms_level {
                            if spectrum.ms_level != filter {
                                continue;
                            }
                        }
                        let raw = text.unescape()?;
                        // There are occasionally empty binary data arrays, or unknown CVs
                        if raw.is_empty() || binary_array.is_none() {
                            continue;
                        }
                        let decoded = base64::decode(raw.as_bytes())?;
                        let bytes = match compression {
                            false => &decoded,
                            true => {
                                let mut r = ZlibDecoder::new(decoded.as_slice());
                                let n = r.read_to_end(&mut output_buffer).await?;
                                &output_buffer[..n]
                            }
                        };

                        let array = match binary_dtype {
                            Dtype::F32 => {
                                let mut buf: [u8; 4] = [0; 4];
                                bytes
                                    .chunks(4)
                                    .filter(|chunk| chunk.len() == 4)
                                    .map(|chunk| {
                                        buf.copy_from_slice(chunk);
                                        f32::from_le_bytes(buf)
                                    })
                                    .collect::<Vec<f32>>()
                            }
                            Dtype::F64 => {
                                let mut buf: [u8; 8] = [0; 8];
                                bytes
                                    .chunks(8)
                                    .map(|chunk| {
                                        buf.copy_from_slice(chunk);
                                        f64::from_le_bytes(buf) as f32
                                    })
                                    .collect::<Vec<f32>>()
                            }
                        };
                        output_buffer.clear();

                        match binary_array {
                            Some(BinaryKind::Intensity) => {
                                spectrum.intensity = array;
                            }
                            Some(BinaryKind::Mz) => {
                                spectrum.mz = array;
                            }
                            None => {}
                        }

                        binary_array = None;
                    }
                }
                Ok(Event::End(ev)) => {
                    state = match (state, ev.name().into_inner()) {
                        (Some(State::Binary), b"binary") => Some(State::BinaryDataArray),
                        (Some(State::BinaryDataArray), b"binaryDataArray") => Some(State::Spectrum),
                        (Some(State::SelectedIon), b"selectedIon") => Some(State::Precursor),
                        (Some(State::Precursor), b"precursor") => {
                            if precursor.mz != 0.0 {
                                precursor.isolation_window = match (iso_window_lo, iso_window_hi) {
                                    (Some(lo), Some(hi)) => Some(Tolerance::Da(-lo, hi)),
                                    _ => None,
                                };
                                spectrum.precursors.push(precursor);
                                precursor = Precursor::default();
                            }
                            Some(State::Spectrum)
                        }
                        (Some(State::Scan), b"scan") => Some(State::Spectrum),
                        (_, b"spectrum") => {
                            let allow = self
                                .ms_level
                                .as_ref()
                                .map(|&level| level == spectrum.ms_level)
                                .unwrap_or(true);
                            if allow {
                                spectra.push(spectrum);
                            }
                            spectrum = Spectrum::default();
                            None
                        }
                        _ => state,
                    };
                }
                Ok(Event::Eof) => break,
                Ok(_) => {}
                Err(err) => {
                    log::error!("unhandled XML error while parsing mzML: {}", err)
                }
            }
            buf.clear();
        }
        Ok(spectra)
    }
}

#[derive(Debug)]
pub enum MzMLError {
    Malformed,
    UnsupportedCV(String),
    XMLError(quick_xml::Error),
    IOError(std::io::Error),
}

impl std::fmt::Display for MzMLError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MzMLError::Malformed => f.write_str("MzMLError: malformed cvParam"),
            MzMLError::UnsupportedCV(s) => write!(f, "MzMLError: unsupported cvParam {}", s),
            MzMLError::IOError(s) => write!(f, "MzMLError: IO error {}", s),
            MzMLError::XMLError(s) => write!(f, "MzMLError: XML error {}", s),
        }
    }
}

impl std::error::Error for MzMLError {}

impl From<std::io::Error> for MzMLError {
    fn from(residual: std::io::Error) -> Self {
        Self::IOError(residual)
    }
}

impl From<quick_xml::Error> for MzMLError {
    fn from(residual: quick_xml::Error) -> Self {
        Self::XMLError(residual)
    }
}

impl From<std::str::Utf8Error> for MzMLError {
    fn from(_: std::str::Utf8Error) -> Self {
        Self::Malformed
    }
}

impl From<std::num::ParseFloatError> for MzMLError {
    fn from(_: std::num::ParseFloatError) -> Self {
        Self::Malformed
    }
}

impl From<std::num::ParseIntError> for MzMLError {
    fn from(_: std::num::ParseIntError) -> Self {
        Self::Malformed
    }
}

impl From<base64::DecodeError> for MzMLError {
    fn from(_: base64::DecodeError) -> Self {
        Self::Malformed
    }
}
