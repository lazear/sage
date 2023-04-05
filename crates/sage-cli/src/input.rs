use clap::ArgMatches;
use sage_cloudpath::CloudPath;
use sage_core::{
    database::{Builder, Parameters},
    lfq::LfqSettings,
    mass::Tolerance,
    tmt::Isobaric,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize)]
/// Actual search parameters - may include overrides or default values not set by user
pub struct Search {
    pub database: Parameters,
    pub quant: QuantSettings,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub min_peaks: usize,
    pub max_peaks: usize,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u16,
    pub report_psms: usize,
    pub predict_rt: bool,
    pub parallel: bool,
    pub mzml_paths: Vec<String>,
    pub output_paths: Vec<String>,

    #[serde(skip_serializing)]
    pub output_directory: CloudPath,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
pub struct Input {
    database: Builder,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    report_psms: Option<usize>,
    chimera: Option<bool>,
    min_peaks: Option<usize>,
    max_peaks: Option<usize>,
    max_fragment_charge: Option<u8>,
    min_matched_peaks: Option<u16>,
    isotope_errors: Option<(i8, i8)>,
    deisotope: Option<bool>,
    quant: Option<QuantOptions>,
    predict_rt: Option<bool>,
    output_directory: Option<String>,
    parallel: Option<bool>,
    mzml_paths: Option<Vec<String>>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct LfqOptions {
    peak_scoring: Option<sage_core::lfq::PeakScoringStrategy>,
    integration: Option<sage_core::lfq::IntegrationStrategy>,
    spectral_angle: Option<f64>,
    ppm_tolerance: Option<f32>,
}

impl From<LfqOptions> for LfqSettings {
    fn from(value: LfqOptions) -> LfqSettings {
        let default = LfqSettings::default();
        let settings = LfqSettings {
            peak_scoring: value.peak_scoring.unwrap_or(default.peak_scoring),
            integration: value.integration.unwrap_or(default.integration),
            spectral_angle: value.spectral_angle.unwrap_or(default.spectral_angle).abs(),
            ppm_tolerance: value.ppm_tolerance.unwrap_or(default.ppm_tolerance).abs(),
        };
        if settings.ppm_tolerance > 20.0 {
            log::warn!("lfq_settings.ppm_tolerance is higher than expected");
        }
        if settings.spectral_angle < 0.50 {
            log::warn!("lfq_settings.spectral_angle is lower than expected");
        }

        settings
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TmtOptions {
    level: Option<u8>,
    sn: Option<bool>,
}

#[derive(Copy, Clone, Serialize)]
pub struct TmtSettings {
    pub level: u8,
    pub sn: bool,
}

impl From<TmtOptions> for TmtSettings {
    fn from(value: TmtOptions) -> Self {
        let default = Self::default();
        Self {
            level: value.level.unwrap_or(default.level),
            sn: value.sn.unwrap_or(default.sn),
        }
    }
}

impl Default for TmtSettings {
    fn default() -> Self {
        Self {
            level: 3,
            sn: false,
        }
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub struct QuantOptions {
    pub tmt: Option<Isobaric>,
    #[serde(rename = "tmt_settings")]
    pub tmt_options: Option<TmtOptions>,

    pub lfq: Option<bool>,
    #[serde(rename = "lfq_settings")]
    pub lfq_options: Option<LfqOptions>,
}

#[derive(Serialize, Default)]
pub struct QuantSettings {
    pub tmt: Option<Isobaric>,
    pub tmt_settings: TmtSettings,
    pub lfq: bool,
    pub lfq_settings: LfqSettings,
}

impl From<QuantOptions> for QuantSettings {
    fn from(value: QuantOptions) -> Self {
        Self {
            tmt: value.tmt,
            tmt_settings: value.tmt_options.map(Into::into).unwrap_or_default(),

            lfq: value.lfq.unwrap_or(false),
            lfq_settings: value.lfq_options.map(Into::into).unwrap_or_default(),
        }
    }
}

impl Input {
    pub fn from_arguments(matches: ArgMatches) -> anyhow::Result<Self> {
        let path = matches
            .get_one::<String>("parameters")
            .expect("required parameters");
        let mut input = Input::load(path)?;

        // Handle JSON configuration overrides
        if let Some(output_directory) = matches.get_one::<String>("output_directory") {
            input.output_directory = Some(output_directory.into());
        }
        if let Some(fasta) = matches.get_one::<String>("fasta") {
            input.database.fasta = Some(fasta.into());
        }
        if let Some(mzml_paths) = matches.get_many::<String>("mzml_paths") {
            input.mzml_paths = Some(mzml_paths.into_iter().map(|p| p.into()).collect());
        }

        if let Some(no_parallel) = matches.get_one::<bool>("no-parallel").copied() {
            if !no_parallel {
                input.parallel = Some(false);
            }
        }

        Ok(input)
    }

    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        sage_core::read_json(path).map_err(anyhow::Error::from)
    }

    fn check_tolerances(tolerance: &Tolerance) {
        match tolerance {
            Tolerance::Ppm(lo, hi) => {
                if hi.abs() > lo.abs() {
                    log::warn!(
                        "Tolerances are applied to experimental masses, not theoretical: [{} - {}]",
                        lo,
                        hi
                    );
                }
            }
            Tolerance::Da(lo, hi) => {
                if hi.abs() > lo.abs() {
                    log::warn!(
                        "Tolerances are applied to experimental masses, not theoretical: [{} - {}]",
                        lo,
                        hi
                    );
                }
            }
        };
    }

    pub fn build(mut self) -> anyhow::Result<Search> {
        let database = self.database.make_parameters();

        Self::check_tolerances(&self.fragment_tol);
        Self::check_tolerances(&self.precursor_tol);

        let isotope_errors = self.isotope_errors.unwrap_or((0, 0));
        if isotope_errors.0 > isotope_errors.1 {
            log::warn!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
        }

        if !self.predict_rt.unwrap_or(true)
            && self.quant.as_ref().and_then(|q| q.lfq).unwrap_or(false)
        {
            log::warn!(
                "`predict_rt: false` and `lfq: true` are incompatible. Setting `predict_rt: true`"
            );
            self.predict_rt = Some(true);
        }

        let mzml_paths = self.mzml_paths.expect("'mzml_paths' must be provided!");

        let output_directory = match self.output_directory {
            Some(path) => {
                let path = path.parse::<CloudPath>()?;
                if let CloudPath::Local(p) = &path {
                    std::fs::create_dir_all(p)?;
                }
                path
            }
            None => CloudPath::Local(std::env::current_dir()?),
        };

        Ok(Search {
            database,
            quant: self.quant.map(Into::into).unwrap_or_default(),
            parallel: self.parallel.unwrap_or(true),
            mzml_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            report_psms: self.report_psms.unwrap_or(1),
            max_peaks: self.max_peaks.unwrap_or(150),
            min_peaks: self.min_peaks.unwrap_or(15),
            min_matched_peaks: self.min_matched_peaks.unwrap_or(4),
            max_fragment_charge: self.max_fragment_charge,
            isotope_errors: self.isotope_errors.unwrap_or((0, 0)),
            deisotope: self.deisotope.unwrap_or(true),
            chimera: self.chimera.unwrap_or(false),
            predict_rt: self.predict_rt.unwrap_or(true),
            output_paths: Vec::new(),
        })
    }
}
