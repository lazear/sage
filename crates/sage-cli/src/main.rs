use clap::{value_parser, Arg, Command, ValueHint};
use sage_cli::input::Input;
use sage_cli::runner::Runner;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("sage")
        .version(clap::crate_version!())
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .about("\u{1F52E} Sage \u{1F9D9} - Proteomics searching so fast it feels like magic!")
        .arg(
            Arg::new("parameters")
                .required(true)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Path to configuration parameters (JSON file)")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("mzml_paths")
                .num_args(1..)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Paths to mzML files to process. Overrides mzML files listed in the \
                     configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to FASTA database. Overrides the FASTA file \
                     specified in the configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("output_directory")
                .short('o')
                .long("output_directory")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path where search and quant results will be written. \
                     Overrides the directory specified in the configuration file.",
                )
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("batch-size")
                .long("batch-size")
                .value_parser(value_parser!(u16).range(1..))
                .help("Number of files to load and search in parallel (default = # of CPUs/2)")
                .value_hint(ValueHint::Other),
        )
        .arg(
            Arg::new("parquet")
                .long("parquet")
                .action(clap::ArgAction::SetTrue)
                .help("Write search output in parquet format instead of tsv"),
        )
        .arg(
            Arg::new("annotate-matches")
                .long("annotate-matches")
                .action(clap::ArgAction::SetTrue)
                .help("Write matched fragments output file."),
        )
        .arg(
            Arg::new("write-pin")
                .long("write-pin")
                .action(clap::ArgAction::SetTrue)
                .help("Write percolator-compatible `.pin` output files"),
        )
        .arg(
            Arg::new("write-report")
                .long("write-report")
                .action(clap::ArgAction::SetTrue)
                .help("Write `.html` report file"),
        )
        .arg(
            Arg::new("disable-telemetry")
                .long("disable-telemetry-i-dont-want-to-improve-sage")
                .action(clap::ArgAction::SetFalse)
                .help("Disable sending telemetry data"),
        )
        .arg(
            Arg::new("save-index")
                .long("save-index")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Save the built index to a parquet directory for reuse")
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("load-index")
                .long("load-index")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Load a pre-built index from parquet directory (skips FASTA processing)")
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("export-index")
                .long("export-index")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Export a user-friendly peptide index parquet file")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("validate-index")
                .long("validate-index")
                .action(clap::ArgAction::SetTrue)
                .help("Validate loaded index matches what would be built from FASTA"),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    let parallel = matches
        .get_one::<u16>("batch-size")
        .copied()
        .unwrap_or_else(|| num_cpus::get() as u16 / 2) as usize;

    let parquet = matches.get_one::<bool>("parquet").copied().unwrap_or(false);
    let send_telemetry = matches
        .get_one::<bool>("disable-telemetry")
        .copied()
        .unwrap_or(true);

    let input = Input::from_arguments(matches)?;

    let runner = input
        .build()
        .and_then(|parameters| Runner::new(parameters, parallel))?;

    let tel = runner.run(parallel, parquet)?;

    if send_telemetry {
        tel.send();
    }

    Ok(())
}
