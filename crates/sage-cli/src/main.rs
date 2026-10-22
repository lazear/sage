use anyhow::Context;
use clap::{value_parser, Arg, Command, ValueHint};
use rayon::ThreadPoolBuilder;
use sage_cli::input::{Input, Template};
use sage_cli::runner::Runner;
use sage_core::tmt::Isobaric;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("sage")
        .version(clap::crate_version!())
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .about("\u{1F52E} Sage \u{1F9D9} - Proteomics searching so fast it feels like magic!")
        .after_help("For detailed documentation please visit https://sage.lazear.org/docs")
        .subcommand_negates_reqs(true)
        .args_conflicts_with_subcommands(true)
        .subcommand(
            Command::new("config")
                .about("Write a minimal configuration file to get started with")
                .args_conflicts_with_subcommands(true)
                .subcommand(
                    Command::new("lint")
                        .about("Lint a configuration file for errors without running a search")
                        .arg(
                            Arg::new("file")
                                .required(true)
                                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                                .help("Path to the configuration file to lint ('-' for stdin)")
                                .value_hint(ValueHint::FilePath),
                        )
                        .arg(
                            Arg::new("output")
                                .short('o')
                                .long("output")
                                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                                .help(
                                    "Write the fully resolved configuration, with every \
                                     default value filled in, to this path ('-' for stdout)",
                                )
                                .value_hint(ValueHint::FilePath),
                        )
                        .arg(
                            Arg::new("force")
                                .long("force")
                                .action(clap::ArgAction::SetTrue)
                                .help("Overwrite the output file if it already exists"),
                        ),
                )
                .arg(
                    Arg::new("mzml_paths")
                        .num_args(1..)
                        .value_parser(clap::builder::NonEmptyStringValueParser::new())
                        .help("Paths to mzML files to search")
                        .value_hint(ValueHint::FilePath),
                )
                .arg(
                    Arg::new("fasta")
                        .short('f')
                        .long("fasta")
                        .value_parser(clap::builder::NonEmptyStringValueParser::new())
                        .help("Path to FASTA database")
                        .value_hint(ValueHint::FilePath),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .value_parser(clap::builder::NonEmptyStringValueParser::new())
                        .help("Path to write the configuration file to [default: stdout]")
                        .value_hint(ValueHint::FilePath),
                )
                .arg(
                    Arg::new("force")
                        .long("force")
                        .action(clap::ArgAction::SetTrue)
                        .help("Overwrite the configuration file if it already exists"),
                )
                .arg(
                    Arg::new("tmt")
                        .long("tmt")
                        .value_name("PLEX")
                        .num_args(0..=1)
                        .require_equals(true)
                        .default_missing_value("tmt16")
                        .value_parser(["tmt6", "tmt10", "tmt11", "tmt16", "tmt18"])
                        .conflicts_with("lfq")
                        .help("Add TMT labels and isobaric quantification"),
                )
                .arg(
                    Arg::new("lfq")
                        .long("lfq")
                        .action(clap::ArgAction::SetTrue)
                        .help("Add label-free quantification"),
                )
                .arg(
                    Arg::new("open")
                        .long("open")
                        .action(clap::ArgAction::SetTrue)
                        .conflicts_with("wide-window")
                        .help("Use a wide precursor tolerance for open searching"),
                )
                .arg(
                    Arg::new("wide-window")
                        .long("wide-window")
                        .action(clap::ArgAction::SetTrue)
                        .help("Search wide-window acquisition, PRM or DIA data"),
                ),
        )
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
            Arg::new("stack-size")
                .long("stack-size")
                .value_parser(value_parser!(u32).range(1..))
                .help("Set Rayon worker thread stack size in MiB (default: 2 MiB)")
                .value_hint(ValueHint::Other),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    if let Some(("config", sub)) = matches.subcommand() {
        if let Some(("lint", lint)) = sub.subcommand() {
            let path = lint.get_one::<String>("file").expect("required");
            let name = match path.as_str() {
                "-" => "stdin".into(),
                path => format!("`{path}`"),
            };
            let (problems, resolved) = Input::lint(path)
                .with_context(|| format!("{name} is not a valid configuration"))?;
            if let Some(output) = lint.get_one::<String>("output") {
                let contents = serde_json::to_string_pretty(&resolved)?;
                write_new(output, &contents, "resolved configuration", lint)?;
            }
            anyhow::ensure!(problems == 0, "found {problems} problem(s) in {name}");
            log::info!("{name} is a valid configuration");
            return Ok(());
        }
        let path = sub.get_one::<String>("output").map_or("-", String::as_str);
        let template = Template {
            fasta: sub.get_one::<String>("fasta").cloned(),
            mzml_paths: sub
                .get_many::<String>("mzml_paths")
                .map(|paths| paths.cloned().collect())
                .unwrap_or_default(),
            tmt: sub
                .get_one::<String>("tmt")
                .map(|plex| match plex.as_str() {
                    "tmt6" => Isobaric::Tmt6,
                    "tmt10" => Isobaric::Tmt10,
                    "tmt11" => Isobaric::Tmt11,
                    "tmt16" => Isobaric::Tmt16,
                    _ => Isobaric::Tmt18,
                }),
            lfq: sub.get_flag("lfq"),
            open: sub.get_flag("open"),
            wide_window: sub.get_flag("wide-window"),
        };
        let contents = serde_json::to_string_pretty(&template.build()?)?;
        write_new(path, &contents, "configuration template", sub)?;
        return Ok(());
    }

    let stack_size_mib = matches.get_one::<u32>("stack-size").copied().unwrap_or(2);
    let stack_size_bytes = stack_size_mib as usize * 1024 * 1024;
    log::trace!(
        "setting Rayon worker thread stack size to {} MiB",
        stack_size_mib
    );
    ThreadPoolBuilder::new()
        .stack_size(stack_size_bytes)
        .build_global()
        .expect("configure Rayon pool");

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

/// Write `contents` to `path` (or stdout, if `path` is `-`), refusing to
/// overwrite an existing file unless `--force` was given
fn write_new(
    path: &str,
    contents: &str,
    what: &str,
    matches: &clap::ArgMatches,
) -> anyhow::Result<()> {
    if path == "-" {
        println!("{contents}");
        return Ok(());
    }
    if std::path::Path::new(path).exists() && !matches.get_flag("force") {
        anyhow::bail!("`{path}` already exists, use --force to overwrite");
    }
    std::fs::write(path, contents.to_owned() + "\n")?;
    log::info!("wrote {what} to `{path}`");
    Ok(())
}
