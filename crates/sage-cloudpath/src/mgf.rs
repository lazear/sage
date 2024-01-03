use regex::Regex;
use sage_core::mass::Tolerance;
use sage_core::spectrum::RawSpectrum;
use sage_core::spectrum::{Precursor, Representation};
use std::panic::Location;

#[derive(Clone)]
pub struct DefaultParams {
    is_query_start: bool,
    regex_for_charge: Regex,
    file_id: usize,
    tol: Option<f32>,
    tol_unit: Option<String>,
    charge_array: Option<Vec<u8>>,
}

impl Default for DefaultParams {
    fn default() -> Self {
        Self {
            is_query_start: false,
            file_id: 0,
            regex_for_charge: Regex::new(r"(\d)\+?").unwrap(),
            tol: None,
            tol_unit: None,
            charge_array: None,
        }
    }
}

impl DefaultParams {
    pub fn default_with_file_id(file_id: usize) -> Self {
        Self {
            file_id,
            ..Default::default()
        }
    }
}

#[derive(Default, Clone)]
pub struct QueryData {
    default_params: DefaultParams,
    spectra: Vec<RawSpectrum>,

    id: String,
    precursors: Vec<Precursor>,
    precursor_tol: Option<f32>,
    precursor_tol_unit: Option<String>,
    precursor_charge_array: Option<Vec<u8>>,
    rt_in_minutes: Option<f32>,
    ion_mz_array: Vec<f32>,
    ion_intensity_array: Vec<f32>,
}

impl QueryData {
    pub fn default_with_params(default_params: DefaultParams) -> Self {
        Self {
            default_params,
            ..Default::default()
        }
    }
    pub fn init(&mut self) {
        self.id = String::default();
        self.precursors = Vec::new();
        self.precursor_tol = self.default_params.tol;
        self.precursor_tol_unit = self.default_params.tol_unit.clone();
        self.precursor_charge_array = self.default_params.charge_array.clone();
        self.rt_in_minutes = None;
        self.ion_mz_array = Vec::new();
        self.ion_intensity_array = Vec::new();
    }

    pub fn get_isolation_window(&mut self) -> Option<Tolerance> {
        if let Some(tol_value) = self.precursor_tol {
            if let Some(unit_str) = self.precursor_tol_unit.as_deref() {
                match unit_str {
                    "Da" => return Some(Tolerance::Da(-tol_value.abs(), tol_value.abs())),
                    "ppm" => return Some(Tolerance::Ppm(-tol_value.abs(), tol_value.abs())),
                    _ => return None,
                }
            }
        }
        None
    }

    //precursors with isolation window and charge
    pub fn get_precursors_with_charge(&mut self) -> Vec<Precursor> {
        let mut new_precursors = Vec::new();
        let isolation_window = self.get_isolation_window();

        for precursor in &mut self.precursors {
            precursor.isolation_window = isolation_window;

            if let Some(charge_array) = &self.precursor_charge_array {
                for &charge in charge_array.iter() {
                    let mut precursor_with_charge = precursor.clone();
                    precursor_with_charge.charge = Some(charge);
                    new_precursors.push(precursor_with_charge);
                }
            } else {
                new_precursors.push(precursor.clone());
            }
        }
        new_precursors
    }

    pub fn default_spectrum(&self, file_id: usize) -> RawSpectrum {
        RawSpectrum {
            file_id,
            ms_level: 2,
            representation: Representation::Centroid,
            ..Default::default()
        }
    }

    pub fn check_spectrum(&self, spectrum: &RawSpectrum) -> Result<bool, MgfError> {
        if spectrum.id.is_empty()
            || spectrum.precursors.is_empty()
            || spectrum.mz.is_empty()
            || spectrum.mz.len() != spectrum.intensity.len()
        {
            return Err(MgfError::Malformed {
                location: *Location::caller(),
            });
        }
        Ok(true)
    }
}

pub struct DefaultParser;

impl DefaultParser {
    pub fn get_parsers(&self) -> Vec<fn(&str, &mut DefaultParams) -> Result<bool, MgfError>> {
        vec![
            Self::parse_begin,
            Self::parse_tol,
            Self::parse_tol_unit,
            Self::parse_charge,
        ]
    }
    pub fn parse_begin(line: &str, default_params: &mut DefaultParams) -> Result<bool, MgfError> {
        if line.starts_with("BEGIN IONS") {
            default_params.is_query_start = true;
            return Ok(true);
        }
        Ok(false)
    }
    pub fn parse_tol(line: &str, default_params: &mut DefaultParams) -> Result<bool, MgfError> {
        if let Some(tol_str) = line.strip_prefix("TOL=") {
            if let Ok(tol) = tol_str.parse::<f32>() {
                default_params.tol = Some(tol);
            }
            return Ok(true);
        }
        Ok(false)
    }
    pub fn parse_tol_unit(
        line: &str,
        default_params: &mut DefaultParams,
    ) -> Result<bool, MgfError> {
        if let Some(tol_unit_str) = line.strip_prefix("TOLU=") {
            default_params.tol_unit = Some(tol_unit_str.to_string());
            return Ok(true);
        }
        Ok(false)
    }
    pub fn parse_charge(line: &str, default_params: &mut DefaultParams) -> Result<bool, MgfError> {
        let regex_for_charge = &default_params.regex_for_charge;

        if let Some(charge_str) = line.strip_prefix("CHARGE=") {
            let mut charge_array: Vec<u8> = Vec::new();
            for cap in regex_for_charge.captures_iter(charge_str) {
                if let Some(charge) = cap[0].chars().next().unwrap().to_digit(10) {
                    charge_array.push(charge as u8);
                }
            }
            default_params.charge_array = Some(charge_array);
            return Ok(true);
        }
        Ok(false)
    }
}

pub struct QueryParser;
impl QueryParser {
    pub fn get_parsers(&self) -> Vec<fn(&str, &mut QueryData) -> Result<bool, MgfError>> {
        vec![
            Self::parse_mz,
            Self::parse_end,
            Self::parse_pepmass,
            Self::parse_title,
            Self::parse_charge,
            Self::parse_tol,
            Self::parse_tol_unit,
            Self::parse_rt,
        ]
    }

    pub fn parse_pepmass(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(pepmass_str) = line.strip_prefix("PEPMASS=") {
            let mut precursor = Precursor::default();
            let mut pepmass = pepmass_str.split_ascii_whitespace();
            if let Some(mz_str) = pepmass.next() {
                match mz_str.parse::<f32>() {
                    Ok(mz) => precursor.mz = mz,
                    Err(_) => {
                        return Err(MgfError::Malformed {
                            location: *Location::caller(),
                        })
                    }
                }
            }
            if let Some(intensity_str) = pepmass.next() {
                if let Ok(intensity) = intensity_str.parse::<f32>() {
                    precursor.intensity = Some(intensity);
                }
            }
            query_data.precursors.push(precursor);
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_charge(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        let regex_for_charge = &query_data.default_params.regex_for_charge;

        if let Some(charge_str) = line.strip_prefix("CHARGE=") {
            let mut charge_array = Vec::new();
            for cap in regex_for_charge.captures_iter(charge_str) {
                if let Some(charge) = cap[0].chars().next().unwrap().to_digit(10) {
                    charge_array.push(charge as u8);
                }
            }
            query_data.precursor_charge_array = Some(charge_array);
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_rt(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(rt_str) = line.strip_prefix("RTINSECONDS=") {
            if let Ok(rt_in_seconds) = rt_str.parse::<f32>() {
                let rt_in_minutes = rt_in_seconds / 60.0;
                query_data.rt_in_minutes = Some(rt_in_minutes);
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub fn parse_title(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(id_str) = line.strip_prefix("TITLE=") {
            query_data.id = id_str.to_string();
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_tol(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(tol_str) = line.strip_prefix("TOL=") {
            if let Ok(tol) = tol_str.parse::<f32>() {
                query_data.precursor_tol = Some(tol);
            }
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_tol_unit(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if let Some(tol_unit_str) = line.strip_prefix("TOLU=") {
            query_data.precursor_tol_unit = Some(tol_unit_str.to_string());
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_mz(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if line.chars().nth(0).unwrap_or_default().is_numeric() {
            let mut mz_intensity = line.split_ascii_whitespace();
            if let Some(mz_str) = mz_intensity.next() {
                match mz_str.parse::<f32>() {
                    Ok(mz) => query_data.ion_mz_array.push(mz),
                    Err(_) => {
                        return Err(MgfError::Malformed {
                            location: *Location::caller(),
                        })
                    }
                }
            }
            if let Some(intensity_str) = mz_intensity.next() {
                if let Ok(intensity) = intensity_str.parse::<f32>() {
                    query_data.ion_intensity_array.push(intensity);
                }
            } else {
                query_data.ion_intensity_array.push(1.0)
            }
            return Ok(true);
        }
        Ok(false)
    }

    pub fn parse_end(line: &str, query_data: &mut QueryData) -> Result<bool, MgfError> {
        if line.starts_with("END IONS") {
            let mut spectrum = query_data.default_spectrum(query_data.default_params.file_id);

            spectrum.id = query_data.id.to_string();
            spectrum.precursors = query_data.get_precursors_with_charge();
            spectrum.scan_start_time = query_data.rt_in_minutes.unwrap_or_default();
            spectrum.total_ion_current = query_data.ion_intensity_array.iter().sum();
            spectrum.mz = std::mem::take(&mut query_data.ion_mz_array);
            spectrum.intensity = std::mem::take(&mut query_data.ion_intensity_array);

            match query_data.check_spectrum(&spectrum) {
                Ok(_) => query_data.spectra.push(spectrum),
                Err(err) => eprintln!("{}", err),
            }
            query_data.init();

            return Ok(true);
        }
        Ok(false)
    }
}

pub struct MgfReader {
    file_id: usize,
}

impl MgfReader {
    pub fn with_file_id(file_id: usize) -> Self {
        Self { file_id }
    }

    pub fn parse(&self, contents: String) -> Result<Vec<RawSpectrum>, MgfError> {
        let default_parsers = DefaultParser.get_parsers();
        let query_parsers = QueryParser.get_parsers();

        let mut default_params = DefaultParams::default_with_file_id(self.file_id);
        let mut lines = contents.as_str().lines();

        // embedded parameters
        while !default_params.is_query_start {
            let line = lines.next().unwrap().trim();
            for parser in &default_parsers {
                match parser(line, &mut default_params) {
                    Ok(true) => break,
                    Ok(false) => continue,
                    Err(err) => eprintln!("{}", err),
                }
            }
        }

        let mut query_data = QueryData::default_with_params(default_params);

        // query
        for line in lines {
            if line.is_empty() {
                continue;
            }
            let line = line.trim();
            for parser in &query_parsers {
                match parser(line, &mut query_data) {
                    Ok(true) => break,
                    Ok(false) => {}
                    Err(err) => eprintln!("{}", err),
                }
            }
        }
        Ok(query_data.spectra)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum MgfError {
    #[error("malformed MGF: {location}")]
    Malformed { location: Location<'static> },
    #[error("unsupported cvParam {0}")]
    UnsupportedCV(String),
    #[error("io error: {0}")]
    IOError(#[from] std::io::Error),
    #[error("utf8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("error parsing float: {0}")]
    FloatError(#[from] std::num::ParseFloatError),
    #[error("error parsing int: {0}")]
    IntError(#[from] std::num::ParseIntError),
    #[error("error decoding base64: {0}")]
    Base64Error(#[from] base64::DecodeError),
}

#[cfg(test)]
mod test {
    use sage_core::{
        mass::Tolerance,
        spectrum::{RawSpectrum, Representation},
    };

    use super::{MgfError, MgfReader};

    fn make_ions_section_spectrum_0() -> String {
        let s = r#"
        BEGIN IONS
        TITLE=spectrum 0
        RTINSECONDS=0.8963232289
        PEPMASS=367.069682741984 56700.5185546875
        CHARGE=2+ and 3+
        TOL=10
        TOLU=ppm
        148.2041016 
        169.5001831 4608.2421875
        226.0483246 5335.4907226563
        228.3407898 30918.244140625
        322.5945435 5311.5737304688
        1144.66272 6260.8315429688
        END IONS
        "#;
        return String::from(s);
    }

    fn run_asserts_for_spectrum_0(s: &RawSpectrum) {
        assert_eq!(s.id, "spectrum 0");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.representation, Representation::Centroid);
        assert_eq!(s.precursors.len(), 2);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert_eq!(s.precursors[1].charge, Some(3));
        assert!((s.precursors[0].mz - 367.069682741984).abs() < 0.0001);
        assert_eq!(s.precursors[0].intensity, Some(56700.5185546875));
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Ppm(-10.0, 10.0))
        );
        assert!((s.precursors[1].mz - 367.069682741984).abs() < 0.0001);
        assert_eq!(s.precursors[1].intensity, Some(56700.5185546875));
        assert_eq!(
            s.precursors[1].isolation_window,
            Some(Tolerance::Ppm(-10.0, 10.0))
        );
        assert!((s.scan_start_time - 0.8963232289 / 60.0).abs() < 0.0001);
        assert_eq!(s.ion_injection_time, 0.0);
        assert_eq!(s.intensity.len(), s.mz.len());
        assert!((s.mz[3] - 228.3407898).abs() < 0.0001);
        assert!((s.intensity[0] - 1.0).abs() < 0.0001);
    }

    #[tokio::test]
    async fn parse_spectrum() -> Result<(), MgfError> {
        let s = make_ions_section_spectrum_0();
        let mut spectra = MgfReader::with_file_id(0).parse(s)?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();

        run_asserts_for_spectrum_0(&s);
        Ok(())
    }

    #[tokio::test]
    async fn parse_two_spectra() -> Result<(), MgfError> {
        let mut content = "# a comment at the beginning of the file".to_string();
        content.push_str(&make_ions_section_spectrum_0());
        content.push_str("\n\n");
        content.push_str(&make_ions_section_spectrum_0());

        let spectra = MgfReader::with_file_id(0).parse(content)?;
        assert_eq!(spectra.len(), 2);
        spectra
            .iter()
            .for_each(|spec: &RawSpectrum| run_asserts_for_spectrum_0(spec));
        Ok(())
    }

    #[tokio::test]
    /// Example taken from https://www.matrixscience.com/help/data_file_help.html
    async fn parse_mgf_matrixscience_example_1() -> Result<(), MgfError> {
        let s = r#"
        COM=10 pmol digest of Sample X15
        ITOL=1
        ITOLU=Da
        MODS=Carbamidomethyl (C)
        IT_MODS=Oxidation (M)
        MASS=Monoisotopic
        USERNAME=Lou Scene
        USEREMAIL=leu@altered-state.edu
        CHARGE=2+ and 3+
        BEGIN IONS
        TITLE=Spectrum 1
        PEPMASS=983.6
        846.60 73
        846.80 44
        847.60 67
        1640.10 291
        1640.60 54
        1895.50 49
        END IONS

        BEGIN IONS
        TITLE=Spectrum 2
        PEPMASS=1084.9
        SCANS=3
        RTINSECONDS=25
        345.10 237
        370.20 128
        460.20 108
        1673.30 1007
        1674.00 974
        1675.30 79
        END IONS
        "#;
        let mut spectra = MgfReader::with_file_id(0).parse(s.to_string())?;
        assert_eq!(spectra.len(), 2);

        let s = spectra.pop().unwrap();
        assert_eq!(s.precursors.len(), 2);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert_eq!(s.precursors[1].charge, Some(3));
        assert_eq!(s.precursors[0].isolation_window, None);
        Ok(())
    }

    #[tokio::test]
    /// Example taken from https://www.matrixscience.com/help/data_file_help.html
    async fn parse_mgf_matrixscience_example_2() -> Result<(), MgfError> {
        let s = r#"
        # following lines define parameters.
        # NB no spaces allowed on either side of the = symbol
        COM=My favourite protein has been eaten by an enzyme
        CLE=Trypsin
        CHARGE=2+
        # following line will be treated as a peptide mass
        1024.6
        # following line is a sequence query, which must
        # conform precisely to sequence query syntax rules
        2321 seq(n-ACTL) comp(2[C])
        # so is this
        1896 ions(345.6:24.7,347.8:45.4, ... ,1024.7:18.7)
        # An MS/MS ions query is delimited by the tags
        # BEGIN IONS and END IONS. Space(s)
        # are used to separate mass and intensity values
        BEGIN IONS
        TITLE=The first peptide - dodgy peak detection, so extra wide tolerance
        PEPMASS=896.05 25674.3
        CHARGE=3+
        TOL=3
        TOLU=Da
        SEQ=n-AC[DHK]
        COMP=2[H]0[M]3[DE]*[K]
        240.1 3
        242.1 12
        245.2 32
        1623.7 55
        1624.7 23
        END IONS
        "#;
        let mut spectra = MgfReader::with_file_id(0).parse(s.to_string())?;
        assert_eq!(spectra.len(), 1);

        let s = spectra.pop().unwrap();
        assert_eq!(s.precursors.len(), 1);
        assert_eq!(s.precursors[0].charge, Some(3));
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-3.0, 3.0))
        );
        Ok(())
    }
}
