use std::io::{BufRead, Read};

use crate::mass::Tolerance;
use crate::spectrum::Precursor;
use flate2::read::ZlibDecoder;
use quick_xml::events::Event;
use quick_xml::Reader;
// use tokio::io::{AsyncBufRead, AsyncRead, AsyncReadExt, AsyncWriteExt, BufReader};

#[derive(Debug)]
pub enum IOError {
    BadCompression,
    Malformed,
    InvalidUri,
    UnsupportedCV(String),
    XMLError(quick_xml::Error),
    IOError(std::io::Error),
    OtherError(Box<dyn std::error::Error + Send + Sync + 'static>),
}

impl std::fmt::Display for IOError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IOError::BadCompression => f.write_str("IOError: bad zlib compression!"),
            IOError::Malformed => f.write_str("IOError: malformed cvParam"),
            IOError::InvalidUri => f.write_str("IOError: invalid URI"),
            IOError::UnsupportedCV(s) => write!(f, "IOError: unsupported cvParam {}", s),
            IOError::IOError(s) => write!(f, "IOError: IO error {}", s),
            IOError::XMLError(s) => write!(f, "IOError: XML error {}", s),
            IOError::OtherError(s) => write!(f, "IOError: XML error {}", s),
        }
    }
}

impl std::error::Error for IOError {}

#[derive(Default, Debug, Clone)]
pub struct Spectrum {
    pub ms_level: u8,
    pub scan_id: usize,
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
pub struct MzMlReader {
    ms_level: Option<u8>,
}

impl MzMlReader {
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
    pub fn parse<B: BufRead>(&self, b: B) -> Result<Vec<Spectrum>, IOError> {
        let mut reader = Reader::from_reader(b);
        let mut buf = Vec::new();

        let mut state = None;
        let mut compression = false;
        let mut output_buffer = Vec::with_capacity(4096);
        let mut binary_dtype = Dtype::F64;
        let mut binary_array = BinaryKind::Intensity;

        let mut spectrum = Spectrum::default();
        let mut precursor = Precursor::default();
        let mut iso_window_lo: Option<f32> = None;
        let mut iso_window_hi: Option<f32> = None;
        let mut spectra = Vec::new();

        let scan_id_regex = regex::Regex::new(r#"scan=(\d+)"#).expect("this is a valid regex");

        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)
                    .map_err(IOError::XMLError)?
                    .ok_or_else(|| IOError::Malformed)?
                    .value
            };
        }

        loop {
            match reader.read_event_into(&mut buf) {
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
                            let id = std::str::from_utf8(&id).map_err(|_| IOError::Malformed)?;
                            match scan_id_regex.captures(id).and_then(|c| c.get(1)) {
                                Some(m) => {
                                    spectrum.scan_id =
                                        m.as_str().parse().map_err(|_| IOError::Malformed)?;
                                }
                                None => {
                                    // fallback, try and extract from index
                                    let index = extract!(ev, b"index");
                                    let index = std::str::from_utf8(&index)
                                        .map_err(|_| IOError::Malformed)?
                                        .parse::<usize>()
                                        .map_err(|_| IOError::Malformed)?;
                                    spectrum.scan_id = index + 1;
                                }
                            }
                        }
                        b"precursor" => {
                            // Not all precursor fields have a spectrumRef
                            if let Some(scan) = ev
                                .try_get_attribute(b"spectrumRef")
                                .map_err(IOError::XMLError)?
                            {
                                let scan = std::str::from_utf8(&scan.value)
                                    .map_err(|_| IOError::Malformed)?;
                                precursor.scan = scan_id_regex
                                    .captures(scan)
                                    .and_then(|c| c.get(1))
                                    .map(|m| m.as_str().parse::<usize>())
                                    .transpose()
                                    .map_err(|_| IOError::Malformed)?;
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::Empty(ref ev)) => match (state, ev.name().into_inner()) {
                    (Some(State::BinaryDataArray), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession =
                            std::str::from_utf8(&accession).map_err(|_| IOError::Malformed)?;
                        match accession {
                            ZLIB_COMPRESSION => compression = true,
                            NO_COMPRESSION => compression = false,
                            FLOAT_64 => binary_dtype = Dtype::F64,
                            FLOAT_32 => binary_dtype = Dtype::F32,
                            INTENSITY_ARRAY => binary_array = BinaryKind::Intensity,
                            MZ_ARRAY => binary_array = BinaryKind::Mz,
                            _ => {
                                return Err(IOError::UnsupportedCV(accession.to_string()));
                            }
                        }
                    }
                    (Some(State::Spectrum), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession =
                            std::str::from_utf8(&accession).map_err(|_| IOError::Malformed)?;
                        match accession {
                            MS_LEVEL => {
                                let level = extract!(ev, b"value");
                                let level = std::str::from_utf8(&level)
                                    .map_err(|_| IOError::Malformed)?
                                    .parse::<u8>()
                                    .map_err(|_| IOError::Malformed)?;
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
                                let value = std::str::from_utf8(&value)
                                    .map_err(|_| IOError::Malformed)?
                                    .parse::<f32>()
                                    .map_err(|_| IOError::Malformed)?;
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
                        let accession =
                            std::str::from_utf8(&accession).map_err(|_| IOError::Malformed)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value).map_err(|_| IOError::Malformed)?;
                        match accession {
                            ISO_WINDOW_LOWER => {
                                iso_window_lo = Some(value.parse().map_err(|_| IOError::Malformed)?)
                            }
                            ISO_WINDOW_UPPER => {
                                iso_window_hi = Some(value.parse().map_err(|_| IOError::Malformed)?)
                            }
                            _ => {}
                        }
                    }
                    (Some(State::SelectedIon), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession =
                            std::str::from_utf8(&accession).map_err(|_| IOError::Malformed)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value).map_err(|_| IOError::Malformed)?;
                        match accession {
                            SELECTED_ION_CHARGE => {
                                precursor.charge =
                                    Some(value.parse().map_err(|_| IOError::Malformed)?)
                            }
                            SELECTED_ION_MZ => {
                                precursor.mz = value.parse().map_err(|_| IOError::Malformed)?
                            }
                            SELECTED_ION_INT => {
                                precursor.intensity =
                                    Some(value.parse().map_err(|_| IOError::Malformed)?)
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Scan), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession =
                            std::str::from_utf8(&accession).map_err(|_| IOError::Malformed)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value).map_err(|_| IOError::Malformed)?;
                        match accession {
                            SCAN_START_TIME => {
                                spectrum.scan_start_time =
                                    value.parse().map_err(|_| IOError::Malformed)?
                            }
                            ION_INJECTION_TIME => {
                                spectrum.ion_injection_time =
                                    value.parse().map_err(|_| IOError::Malformed)?
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
                        let raw = text.unescape().map_err(IOError::XMLError)?;
                        // There are occasionally empty binary data arrays...
                        if raw.is_empty() {
                            continue;
                        }
                        let decoded =
                            base64::decode(raw.as_bytes()).map_err(|_| IOError::Malformed)?;
                        let bytes = match compression {
                            false => &decoded,
                            true => {
                                let mut r = ZlibDecoder::new(decoded.as_slice());
                                let n = r
                                    .read_to_end(&mut output_buffer)
                                    .map_err(IOError::IOError)?;
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
                            BinaryKind::Intensity => {
                                spectrum.intensity = array;
                            }
                            BinaryKind::Mz => {
                                spectrum.mz = array;
                            }
                        }
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
                    dbg!(err);
                }
            }
            buf.clear();
        }
        Ok(spectra)
    }
}
