pub mod database;
pub mod enzyme;
pub mod fasta;
pub mod fdr;
pub mod heap;
pub mod ion_series;
pub mod isotopes;
pub mod lfq;
pub mod mass;
pub mod ml;
pub mod modification;
pub mod mzml;
pub mod peptide;
pub mod scoring;
pub mod spectrum;
pub mod tmt;

use sage_cloudpath::read_and_execute;
use tokio::io::AsyncReadExt;

#[derive(Debug)]
pub enum Error {
    CloudPathError(sage_cloudpath::Error),
    MzML(mzml::MzMLError),
    Json(serde_json::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CloudPathError(e) => e.fmt(f),
            Self::MzML(e) => e.fmt(f),
            Self::Json(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for Error {}

pub fn read_mzml<S: AsRef<str>>(
    s: S,
    signal_to_noise: Option<u8>,
) -> Result<Vec<mzml::Spectrum>, Error> {
    let res = read_and_execute(s, |bf| async move {
        Ok(mzml::MzMLReader::default()
            .set_signal_to_noise(signal_to_noise)
            .parse(bf)
            .await)
    });

    match res {
        Ok(t) => t.map_err(Error::MzML),
        Err(e) => Err(Error::CloudPathError(e)),
    }
}

pub fn read_fasta<S>(path: S, decoy_tag: S, generate_decoys: bool) -> Result<fasta::Fasta, Error>
where
    S: AsRef<str>,
{
    let res = read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(sage_cloudpath::Error::IO)?;
        Ok(fasta::Fasta::parse(
            contents,
            decoy_tag.as_ref(),
            generate_decoys,
        ))
    });

    match res {
        Ok(t) => Ok(t),
        Err(e) => Err(Error::CloudPathError(e)),
    }
}

pub fn read_json<S, T>(path: S) -> Result<T, Error>
where
    S: AsRef<str>,
    T: for<'de> serde::Deserialize<'de>,
{
    let res = read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(sage_cloudpath::Error::IO)?;
        Ok(serde_json::from_str(&contents))
    });

    match res {
        Ok(t) => t.map_err(Error::Json),
        Err(e) => Err(Error::CloudPathError(e)),
    }
}
