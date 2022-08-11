use quick_xml::events::Event;
use quick_xml::Reader;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum MzMlError {
    BadCompression,
    Malformed,
    UnsupportedCV(String),
}

impl std::fmt::Display for MzMlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MzMlError::BadCompression => f.write_str("MzMlError: bad zlib compression!"),
            MzMlError::Malformed => f.write_str("MzMlError: malformed cvParam"),
            MzMlError::UnsupportedCV(s) => write!(f, "MzMlError: unsupported cvParam {}", s),
        }
    }
}

impl std::error::Error for MzMlError {}

#[derive(Default, Debug, Copy, Clone)]
pub struct Precursor {
    pub mz: f32,
    pub intensity: Option<f32>,
    pub charge: Option<u8>,
    pub scan: Option<usize>,
}

#[derive(Default, Debug, Clone)]
pub struct Spectrum {
    pub ms_level: usize,
    pub scan_id: usize,
    pub precursor: Vec<Precursor>,
    pub representation: Representation,

    // Scan start time
    pub scan_start_time: f32,
    // Ion injection time
    pub ion_injection_time: f32,
    // M/z array
    pub mz: Vec<f32>,
    // Intensity array
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

const SELECTED_ION_MZ: &str = "MS:1000744";
const SELECTED_ION_INT: &str = "MS:1000042";
const SELECTED_ION_CHARGE: &str = "MS:1000041";

#[derive(Default)]
pub struct MzMlReader {
    ms_level: Option<usize>,
}

impl MzMlReader {
    /// Create a new [`MzMlReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_level_filter(ms_level: usize) -> Self {
        Self {
            ms_level: Some(ms_level),
        }
    }

    /// Convenience method
    pub fn read_ms2<P: AsRef<Path>>(
        p: P,
    ) -> Result<Vec<Spectrum>, Box<dyn std::error::Error + Send + Sync + 'static>> {
        let file = std::fs::File::open(p)?;
        let reader = BufReader::new(file);
        Self::with_level_filter(2).parse(reader)
    }

    pub fn parse<B: BufRead>(
        &self,
        b: B,
    ) -> Result<Vec<Spectrum>, Box<dyn std::error::Error + Send + Sync + 'static>> {
        let mut reader = Reader::from_reader(b);
        let mut buf = Vec::new();

        let mut state = None;
        let mut compression = false;
        let mut binary_dtype = Dtype::F64;
        let mut binary_array = BinaryKind::Intensity;

        let mut spectrum = Spectrum::default();
        let mut precursor = Precursor::default();
        let mut spectra = Vec::new();

        let scan_id_regex = regex::Regex::new(r#"scan=(\d+)"#)?;

        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)?
                    .ok_or_else(|| MzMlError::Malformed)?
                    .value
            };
        }

        loop {
            match reader.read_event(&mut buf) {
                Ok(Event::Start(ref ev)) => {
                    // State transition into child tag
                    state = match (ev.name(), state) {
                        (b"spectrum", _) => Some(State::Spectrum),
                        (b"scan", Some(State::Spectrum)) => Some(State::Scan),
                        (b"binaryDataArray", Some(State::Spectrum)) => Some(State::BinaryDataArray),
                        (b"binary", Some(State::BinaryDataArray)) => Some(State::Binary),
                        (b"selectedIon", Some(State::Spectrum)) => Some(State::SelectedIon),
                        _ => state,
                    };
                    match ev.name() {
                        b"spectrum" => {
                            let id = extract!(ev, b"id");
                            let id = std::str::from_utf8(&id)?;
                            match scan_id_regex.captures(id).and_then(|c| c.get(1)) {
                                Some(m) => {
                                    spectrum.scan_id = m.as_str().parse()?;
                                }
                                None => {
                                    // fallback, try and extract from index
                                    let index = extract!(ev, b"index");
                                    let index = std::str::from_utf8(&index)?.parse::<usize>()?;
                                    spectrum.scan_id = index + 1;
                                }
                            }
                        }
                        b"precursor" => {
                            // Not all precursor fields have a spectrumRef
                            if let Some(scan) = ev.try_get_attribute(b"spectrumRef")? {
                                let scan = std::str::from_utf8(&scan.value)?;
                                precursor.scan = scan_id_regex
                                    .captures(scan)
                                    .and_then(|c| c.get(1))
                                    .map(|m| m.as_str().parse::<usize>())
                                    .transpose()?;
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::Empty(ref ev)) => match (state, ev.name()) {
                    (Some(State::BinaryDataArray), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        match accession {
                            ZLIB_COMPRESSION => compression = true,
                            NO_COMPRESSION => compression = false,
                            FLOAT_64 => binary_dtype = Dtype::F64,
                            FLOAT_32 => binary_dtype = Dtype::F32,
                            INTENSITY_ARRAY => binary_array = BinaryKind::Intensity,
                            MZ_ARRAY => binary_array = BinaryKind::Mz,
                            _ => {
                                return Err(Box::new(MzMlError::UnsupportedCV(
                                    accession.to_string(),
                                )))
                            }
                        }
                    }
                    (Some(State::Spectrum), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        match accession {
                            MS_LEVEL => {
                                let level = extract!(ev, b"value");
                                let level = std::str::from_utf8(&level)?.parse::<usize>()?;
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
                            _ => {}
                        }
                    }
                    (Some(State::SelectedIon), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value)?;
                        match accession {
                            SELECTED_ION_CHARGE => precursor.charge = Some(value.parse()?),
                            SELECTED_ION_MZ => precursor.mz = value.parse()?,
                            SELECTED_ION_INT => precursor.intensity = Some(value.parse()?),
                            _ => {}
                        }
                    }
                    (Some(State::Scan), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        let accession = std::str::from_utf8(&accession)?;
                        let value = extract!(ev, b"value");
                        let value = std::str::from_utf8(&value)?;
                        match accession {
                            "MS:1000016" => spectrum.scan_start_time = value.parse()?,
                            "MS:1000927" => spectrum.ion_injection_time = value.parse()?,
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
                        let raw = text.unescaped()?;
                        let decoded = base64::decode(raw)?;
                        let bytes = match compression {
                            false => decoded,
                            true => miniz_oxide::inflate::decompress_to_vec_zlib(&decoded)
                                .map_err(|_| MzMlError::BadCompression)?,
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
                    state = match (state, ev.name()) {
                        (Some(State::Binary), b"binary") => Some(State::BinaryDataArray),
                        (Some(State::BinaryDataArray), b"binaryDataArray") => Some(State::Spectrum),
                        (Some(State::SelectedIon), b"selectedIon") => {
                            if precursor.mz != 0.0 {
                                spectrum.precursor.push(precursor);
                                precursor = Precursor::default();
                            }
                            Some(State::Spectrum)
                        }
                        (Some(State::Scan), b"scan") => Some(State::Spectrum),
                        (_, b"spectrum") => {
                            match self.ms_level {
                                Some(filter) => {
                                    if spectrum.ms_level == filter {
                                        spectra.push(spectrum);
                                    }
                                }
                                _ => {
                                    spectra.push(spectrum);
                                }
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
