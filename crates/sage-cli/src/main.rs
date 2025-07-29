use clap::{value_parser, Arg, Command, ValueHint};
use sage_cli::input::Input;
use sage_cli::runner::Runner;
use rayon::ThreadPoolBuilder;   // ← add this to increase stack size to fix stack overflow

fn main() -> anyhow::Result<()> {

    // ── enlarge Rayon worker-thread stacks ───────────────────────────
    ThreadPoolBuilder::new()
        .stack_size(64 * 1024 * 1024)   // 64 MiB per worker thread
        .build_global()
        .expect("configure Rayon pool");

    // ── the rest of your original startup code ───────────────────────

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
            Arg::new("disable-telemetry")
                .long("disable-telemetry-i-dont-want-to-improve-sage")
                .action(clap::ArgAction::SetFalse)
                .help("Disable sending telemetry data"),
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
