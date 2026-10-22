use anyhow::{ensure, Context};
use clap::ArgMatches;
use sage_cloudpath::tdf::BrukerProcessingConfig;
use sage_cloudpath::Url;
use sage_core::scoring::ScoreType;
use sage_core::{
    database::{Builder, Parameters},
    lfq::LfqSettings,
    mass::Tolerance,
    tmt::Isobaric,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Clone)]
/// Actual search parameters - may include overrides or default values not set by user
pub struct Search {
    pub version: String,
    pub database: Parameters,
    pub quant: QuantSettings,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub wide_window: bool,
    pub min_peaks: usize,
    pub max_peaks: usize,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u16,
    pub report_psms: usize,
    pub predict_rt: bool,
    pub mzml_paths: Vec<Url>,
    pub output_paths: Vec<Url>,
    pub bruker_config: BrukerProcessingConfig,
    pub protein_grouping: bool,
    pub protein_grouping_peptide_fdr: f32,

    #[serde(skip_serializing)]
    pub output_directory: Url,

    #[serde(skip_serializing)]
    pub write_pin: bool,

    #[serde(skip_serializing)]
    pub write_report: bool,

    #[serde(skip_serializing)]
    pub annotate_matches: bool,

    pub score_type: ScoreType,
}

/// Like [`serde_json::to_value`], but `f32`s are written using their shortest
/// representation (`0.01`) instead of being widened (`0.009999999776482582`)
fn to_json<T: Serialize>(value: &T) -> serde_json::Result<serde_json::Value> {
    serde_json::from_str(&serde_json::to_string(value)?)
}

#[derive(Default)]
/// Presets for the configuration template written by `sage config`
pub struct Template {
    /// Placeholder paths are used if these are not set
    pub fasta: Option<String>,
    pub mzml_paths: Vec<String>,
    pub tmt: Option<Isobaric>,
    pub lfq: bool,
    pub open: bool,
    pub wide_window: bool,
}

impl Template {
    /// Minimal configuration, where everything not specified takes Sage's
    /// default value
    pub fn build(&self) -> anyhow::Result<serde_json::Value> {
        let mut config = serde_json::json!({
            "database": {
                "fasta": "path/to/database.fasta",
                "enzyme": {
                    "missed_cleavages": 2,
                    "min_len": 7,
                    "max_len": 35,
                },
                "static_mods": { "C": {"mass": 57.0215, "name": "Unimod:4"} },
                "variable_mods": { "M": [{"mass": 15.9949, "name": "Unimod:35"}] },
            },
            "precursor_tol": { "ppm": [-20.0, 20.0] },
            "fragment_tol": { "ppm": [-20.0, 20.0] },
            "isotope_errors": [-1, 3],
            "mzml_paths": ["path/to/file.mzML"],
        });

        if let Some(fasta) = &self.fasta {
            config["database"]["fasta"] = fasta.as_str().into();
        }
        if !self.mzml_paths.is_empty() {
            config["mzml_paths"] = self.mzml_paths.clone().into();
        }

        if let Some(tmt) = &self.tmt {
            // Labels on lysine and peptide N-termini: TMTpro for 16/18-plex,
            // otherwise TMT 6/10/11-plex
            let label = match tmt {
                Isobaric::Tmt16 | Isobaric::Tmt18 => {
                    serde_json::json!({"mass": 304.207146, "name": "Unimod:2016"})
                }
                _ => serde_json::json!({"mass": 229.162932, "name": "Unimod:737"}),
            };
            config["database"]["static_mods"]["K"] = label.clone();
            config["database"]["static_mods"]["^"] = label;
            config["quant"]["tmt"] = to_json(tmt)?;
            config["quant"]["tmt_settings"] = to_json(&TmtSettings::default())?;
        }
        if self.lfq {
            config["quant"]["lfq"] = true.into();
            config["quant"]["lfq_settings"] = to_json(&LfqSettings::default())?;
        }
        if self.open {
            config["precursor_tol"] = serde_json::json!({ "da": [-500.0, 100.0] });
            config["isotope_errors"] = serde_json::json!([0, 0]);
        }
        if self.wide_window {
            config["wide_window"] = true.into();
            config["chimera"] = true.into();
            config["report_psms"] = 5.into();
        }
        Ok(config)
    }
}

impl Search {
    /// Fully resolved search parameters for the default [`Template`]
    pub fn defaults() -> anyhow::Result<Self> {
        let mut input: Input = serde_json::from_value(Template::default().build()?)?;
        // Placeholder paths can't be resolved
        input.mzml_paths = Some(Vec::new());
        input.build()
    }
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
pub struct Input {
    pub database: Builder,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub report_psms: Option<usize>,
    pub chimera: Option<bool>,
    pub wide_window: Option<bool>,
    pub min_peaks: Option<usize>,
    pub max_peaks: Option<usize>,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: Option<u16>,
    pub precursor_charge: Option<(u8, u8)>,
    pub override_precursor_charge: Option<bool>,
    pub isotope_errors: Option<(i8, i8)>,
    pub deisotope: Option<bool>,
    pub quant: Option<QuantOptions>,
    pub predict_rt: Option<bool>,
    pub output_directory: Option<String>,
    pub mzml_paths: Option<Vec<String>>,
    pub bruker_config: Option<BrukerProcessingConfig>,
    pub protein_grouping: Option<bool>,
    pub protein_grouping_peptide_fdr: Option<f32>,

    pub annotate_matches: Option<bool>,
    pub write_pin: Option<bool>,
    pub write_report: Option<bool>,
    pub score_type: Option<ScoreType>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct LfqOptions {
    pub peak_scoring: Option<sage_core::lfq::PeakScoringStrategy>,
    pub integration: Option<sage_core::lfq::IntegrationStrategy>,
    pub spectral_angle: Option<f64>,
    pub ppm_tolerance: Option<f32>,
    pub mobility_pct_tolerance: Option<f32>,
    pub combine_charge_states: Option<bool>,
    pub peptide_q_value: Option<f32>,
}

impl From<LfqOptions> for LfqSettings {
    fn from(value: LfqOptions) -> LfqSettings {
        let default = LfqSettings::default();
        let settings = LfqSettings {
            peak_scoring: value.peak_scoring.unwrap_or(default.peak_scoring),
            integration: value.integration.unwrap_or(default.integration),
            spectral_angle: value.spectral_angle.unwrap_or(default.spectral_angle).abs(),
            ppm_tolerance: value.ppm_tolerance.unwrap_or(default.ppm_tolerance).abs(),
            peptide_q_value: value.peptide_q_value.unwrap_or(default.peptide_q_value),
            mobility_pct_tolerance: value
                .mobility_pct_tolerance
                .unwrap_or(default.mobility_pct_tolerance),
            combine_charge_states: value
                .combine_charge_states
                .unwrap_or(default.combine_charge_states),
        };
        if settings.ppm_tolerance > 20.0 {
            log::warn!("lfq_settings.ppm_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance > 4.0 {
            log::warn!("lfq_settings.mobility_pct_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance < 0.05 {
            log::warn!("lfq_settings.mobility_pct_tolerance is smaller than expected");
        }
        if settings.spectral_angle < 0.50 {
            log::warn!("lfq_settings.spectral_angle is lower than expected");
        }
        if settings.peptide_q_value > 0.01 {
            log::info!("lfq_settings.peptide_q_value is higher than expected, expect increased runtime and memory usage");
        }
        if settings.peptide_q_value < 0.01 {
            log::warn!("lfq_settings.peptide_q_value is lower than expected, not all identified peptides will have MS1 intensities extracted");
        }

        settings
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TmtOptions {
    pub level: Option<u8>,
    pub sn: Option<bool>,
}

#[derive(Copy, Clone, Serialize, Debug)]
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

#[derive(Serialize, Default, Clone)]
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
        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

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

        if let Some(write_pin) = matches.get_one::<bool>("write-pin").copied() {
            input.write_pin = Some(write_pin);
        }

        if let Some(write_report) = matches.get_one::<bool>("write-report").copied() {
            input.write_report = Some(write_report);
        }

        if let Some(annotate_matches) = matches.get_one::<bool>("annotate-matches").copied() {
            input.annotate_matches = Some(annotate_matches);
        }

        // avoid to later panic if these parameters are not set (but doesn't check if files exist)

        ensure!(
            input.database.fasta.is_some(),
            "`database.fasta` must be set. For more information try '--help'"
        );
        ensure!(
            input
                .mzml_paths
                .as_ref()
                .map(|p| p.len())
                .unwrap_or_default()
                > 0,
            "`mzml_paths` must be set. For more information try '--help'"
        );

        Ok(input)
    }

    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        sage_cloudpath::util::read_json(path).map_err(anyhow::Error::from)
    }

    /// Like [`Input::load`], but also returns the paths of any keys that are
    /// not recognized (and would otherwise be silently ignored). Reads from
    /// stdin if `path` is `-`
    pub fn load_checked<S: AsRef<str>>(path: S) -> anyhow::Result<(Self, Vec<String>)> {
        let value = match path.as_ref() {
            "-" => serde_json::from_reader(std::io::stdin().lock())?,
            path => sage_cloudpath::util::read_json(path)?,
        };
        Self::from_value_checked(value)
    }

    fn from_value_checked(value: serde_json::Value) -> anyhow::Result<(Self, Vec<String>)> {
        let mut unknown = Vec::new();
        let input = serde_ignored::deserialize(value, |path| {
            // Options are rendered as `?` path segments
            let path = path.to_string().replace(".?", "");
            // Output-only fields, present when re-using a `results.json`
            if !matches!(path.as_str(), "version" | "output_paths") {
                unknown.push(path);
            }
        })?;
        Ok((input, unknown))
    }

    /// Validate a configuration file without running a search, returning the
    /// number of problems found and the fully resolved configuration. Paths
    /// that may be supplied on the command line (`database.fasta`,
    /// `mzml_paths`) are allowed to be missing.
    pub fn lint<S: AsRef<str>>(path: S) -> anyhow::Result<(usize, serde_json::Value)> {
        let (mut input, unknown) = Self::load_checked(path)?;
        let fasta = input.database.fasta.clone();
        let mzml_paths = input.mzml_paths.clone();
        let output_directory = input.output_directory.take();
        let mut problems = unknown.len();
        for key in &unknown {
            log::error!("unrecognized parameter `{key}`");
        }

        let exists = |path: &str| match sage_cloudpath::try_parse_url(path) {
            Some(url) if url.scheme() != "file" => true,
            Some(url) => url.to_file_path().map_or(false, |p| p.exists()),
            None => std::path::Path::new(path).exists(),
        };

        let mut scan_fasta = false;
        match &input.database.fasta {
            Some(fasta) if !exists(fasta) => {
                log::error!("`database.fasta` does not exist: {fasta}");
                problems += 1;
            }
            Some(_) => scan_fasta = true,
            None => {
                log::warn!("`database.fasta` is not set, it must be passed with --fasta");
                input.database.fasta = Some(String::new());
            }
        }
        match &mut input.mzml_paths {
            Some(paths) if !paths.is_empty() => {
                // Drop missing files, since resolving them in `build` would fail
                paths.retain(|path| {
                    let ok = exists(path);
                    if !ok {
                        log::error!("`mzml_paths` entry does not exist: {path}");
                        problems += 1;
                    }
                    ok
                });
            }
            _ => {
                log::warn!(
                    "`mzml_paths` is not set, input files must be passed on the command line"
                );
                input.mzml_paths = Some(Vec::new());
            }
        }

        // Don't create the output directory as a side effect of linting
        let search = input.build()?;
        if scan_fasta {
            problems += Self::check_decoys(&search.database)?;
        }

        let mut resolved = to_json(&search)?;
        let config = resolved
            .as_object_mut()
            .expect("struct serializes to object");
        // Output-only fields
        config.shift_remove("version");
        config.shift_remove("output_paths");
        // Report paths as written, rather than placeholders or resolved URLs
        match fasta {
            Some(fasta) => resolved["database"]["fasta"] = fasta.into(),
            None => {
                resolved["database"]
                    .as_object_mut()
                    .map(|db| db.shift_remove("fasta"));
            }
        }
        match mzml_paths {
            Some(paths) => resolved["mzml_paths"] = paths.into(),
            None => {
                resolved
                    .as_object_mut()
                    .map(|c| c.shift_remove("mzml_paths"));
            }
        }
        // Not serialized as part of `Search`
        if let Some(output_directory) = output_directory {
            resolved["output_directory"] = output_directory.into();
        }
        resolved["write_pin"] = search.write_pin.into();
        resolved["write_report"] = search.write_report.into();
        resolved["annotate_matches"] = search.annotate_matches.into();

        Ok((problems, resolved))
    }

    /// Make sure the FASTA contains target sequences, and that decoys are
    /// either present in the database or will be generated by Sage
    fn check_decoys(database: &Parameters) -> anyhow::Result<usize> {
        let url = sage_cloudpath::to_url(&database.fasta)?;
        let tag = &database.decoy_tag;
        // Don't let `read_fasta` filter out decoys, so we can count them
        let fasta = sage_cloudpath::util::read_fasta(&url, tag, false)
            .with_context(|| format!("Failed to read `{}`", database.fasta))?;
        let total = fasta.targets.len();
        let decoys = fasta
            .targets
            .iter()
            .filter(|(acc, _)| acc.contains(tag.as_str()))
            .count();

        if total == 0 {
            log::error!("`database.fasta` contains no protein sequences");
            return Ok(1);
        }
        if decoys == total {
            log::error!(
                "all {total} entries in `database.fasta` contain the decoy tag `{tag}`, \
                 leaving no target sequences - check `database.decoy_tag`"
            );
            return Ok(1);
        }
        match (database.generate_decoys, decoys) {
            (true, 0) => log::info!("{total} target proteins, decoys will be generated"),
            (true, _) => log::info!(
                "{decoys} entries containing the decoy tag `{tag}` will be ignored, \
                 and decoys generated from the remaining {} target proteins",
                total - decoys
            ),
            (false, 0) => {
                log::error!(
                    "`database.generate_decoys` is false, but no entries in `database.fasta` \
                     contain the decoy tag `{tag}` - FDR cannot be estimated"
                );
                return Ok(1);
            }
            (false, _) => log::info!(
                "{} target and {decoys} decoy proteins (tag `{tag}`) in `database.fasta`",
                total - decoys
            ),
        }
        Ok(0)
    }

    fn check_mass_tolerances(tolerance: &Tolerance) {
        let (lo, hi) = match tolerance {
            Tolerance::Ppm(lo, hi) => (*lo, *hi),
            Tolerance::Pct(lo, hi) => {
                log::warn!(
                    "Pct tolerances are very rarely used for mass tolerances, did you mean ppm?"
                );
                (*lo, *hi)
            }
            Tolerance::Da(lo, hi) => (*lo, *hi),
        };
        if hi.abs() > lo.abs() {
            log::warn!(
                "Tolerances are applied to experimental masses, not theoretical: [{}, {}]. Example open search config: [-800, 150]",
                lo,
                hi
            );
        }
        if lo > 0.0 {
            log::warn!(
                "The `left` tolerance should probably be negative, for example: [{}, {}]",
                -lo,
                hi.abs()
            )
        }
        if hi < 0.0 {
            log::warn!(
                "The `right` tolerance should probably be positive, for example: [{}, {}]",
                -lo.abs(),
                hi
            )
        }
    }

    pub fn build(mut self) -> anyhow::Result<Search> {
        let database = self
            .database
            .make_parameters()
            .context("invalid modification configuration")?;

        Self::check_mass_tolerances(&self.fragment_tol);
        Self::check_mass_tolerances(&self.precursor_tol);

        if let Some(isotope_errors) = self.isotope_errors {
            if isotope_errors.0 > isotope_errors.1 {
                log::error!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
                std::process::exit(1);
            }
        }
        if let Some(charges) = self.precursor_charge {
            if charges.0 > charges.1 {
                log::error!(
                    "Precursor charges should be specified [low, high], user provided: [{}, {}]",
                    charges.0,
                    charges.1
                );
                std::process::exit(1);
            }
        }

        if !self.predict_rt.unwrap_or(true)
            && self.quant.as_ref().and_then(|q| q.lfq).unwrap_or(false)
        {
            log::warn!(
                "`predict_rt: false` and `lfq: true` are incompatible. Setting `predict_rt: true`"
            );
            self.predict_rt = Some(true);
        }

        let mzml_paths = self
            .mzml_paths
            .expect("'mzml_paths' must be provided!")
            .iter()
            .map(|s| sage_cloudpath::to_url(s))
            .collect::<Result<Vec<_>, _>>()?;

        let output_directory = match self.output_directory {
            Some(path) => {
                match sage_cloudpath::try_parse_url(&path) {
                    Some(mut url) => {
                        // Valid URL, might still be a local directory that doesn't exist
                        if url.scheme() == "file" {
                            let path = url.to_file_path().expect("url scheme is file");
                            std::fs::create_dir_all(path)?;
                        }

                        if !url.path().ends_with("/") {
                            url.set_path(&format!("{}/", url.path()));
                        }
                        url
                    }
                    None => {
                        // Treat as a local path (covers Windows `C:\...` which
                        // otherwise parses as a URL with scheme `c`).
                        let path = std::path::Path::new(&path);
                        std::fs::create_dir_all(path)?;
                        Url::from_directory_path(path.canonicalize()?).expect("valid path")
                    }
                }
            }
            None => {
                let dir = std::env::current_dir()?;
                Url::from_directory_path(dir).expect("valid path")
            }
        };

        let score_type = self.score_type.unwrap_or(ScoreType::SageHyperScore);

        Ok(Search {
            version: clap::crate_version!().into(),
            database,
            quant: self.quant.map(Into::into).unwrap_or_default(),
            mzml_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            report_psms: self.report_psms.unwrap_or(1),
            max_peaks: self.max_peaks.unwrap_or(150),
            min_peaks: self.min_peaks.unwrap_or(15),
            min_matched_peaks: self.min_matched_peaks.unwrap_or(4),
            max_fragment_charge: self.max_fragment_charge,
            annotate_matches: self.annotate_matches.unwrap_or(false),
            precursor_charge: self.precursor_charge.unwrap_or((2, 4)),
            override_precursor_charge: self.override_precursor_charge.unwrap_or(false),
            isotope_errors: self.isotope_errors.unwrap_or((0, 0)),
            deisotope: self.deisotope.unwrap_or(true),
            chimera: self.chimera.unwrap_or(false),
            wide_window: self.wide_window.unwrap_or(false),
            predict_rt: self.predict_rt.unwrap_or(true),
            output_paths: Vec::new(),
            write_pin: self.write_pin.unwrap_or(false),
            bruker_config: self.bruker_config.unwrap_or_default(),
            write_report: self.write_report.unwrap_or(false),
            protein_grouping: self.protein_grouping.unwrap_or(true),
            protein_grouping_peptide_fdr: self.protein_grouping_peptide_fdr.unwrap_or(0.01),
            score_type,
        })
    }
}

#[cfg(test)]
mod test {

    use sage_core::{database::EnzymeBuilder, enzyme::EnzymeParameters};

    #[test]
    fn deserialize_enzyme_builder() -> Result<(), serde_json::Error> {
        let a: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
        }))?;
        let b: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "P",
        }))?;
        let c: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "",
        }))?;

        let a: EnzymeParameters = a.into();
        let b: EnzymeParameters = b.into();
        let c: EnzymeParameters = c.into();

        assert_eq!(a.enzyme.map(|e| e.skip_suffix), Some([false; 26]));
        {
            let mut expected = [false; 26];
            expected[(b'P' - b'A') as usize] = true;
            assert_eq!(b.enzyme.map(|e| e.skip_suffix), Some(expected));
        }
        assert_eq!(c.enzyme.map(|e| e.skip_suffix), Some([false; 26]));

        Ok(())
    }

    #[test]
    fn unknown_keys_are_reported() -> anyhow::Result<()> {
        use super::{Input, Search};

        let mut config = serde_json::to_value(Search::defaults()?)?;
        config["deisotop"] = true.into();
        config["database"]["enzyme"]["misssed_cleavages"] = 2.into();
        config["quant"]["lfq_setings"] = serde_json::json!({});

        let path = std::env::temp_dir().join("sage-unknown-keys.json");
        std::fs::write(&path, config.to_string())?;
        let (_, mut unknown) = Input::load_checked(path.to_str().unwrap())?;
        std::fs::remove_file(&path)?;

        unknown.sort();
        assert_eq!(
            unknown,
            [
                "database.enzyme.misssed_cleavages",
                "deisotop",
                "quant.lfq_setings"
            ]
        );
        Ok(())
    }

    #[test]
    fn templates_are_valid() -> anyhow::Result<()> {
        use super::{Input, Template};
        use sage_core::tmt::Isobaric;

        for tmt in [None, Some(Isobaric::Tmt11), Some(Isobaric::Tmt18)] {
            for (lfq, open, wide_window) in [
                (false, false, false),
                (true, false, false),
                (false, true, false),
                (false, false, true),
            ] {
                let template = Template {
                    tmt: tmt.clone(),
                    lfq,
                    open,
                    wide_window,
                    ..Default::default()
                };
                let (mut input, unknown) = Input::from_value_checked(template.build()?)?;
                assert!(unknown.is_empty(), "{unknown:?}");
                input.mzml_paths = Some(Vec::new());
                input.build()?;
            }
        }
        Ok(())
    }

    #[test]
    fn default_search_round_trips() -> anyhow::Result<()> {
        use super::{Input, Search};

        let defaults = serde_json::to_value(Search::defaults()?)?;
        let input: Input = serde_json::from_value(defaults.clone())?;
        assert_eq!(serde_json::to_value(input.build()?)?, defaults);
        Ok(())
    }
}
