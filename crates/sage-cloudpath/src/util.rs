use crate::read_and_execute;
use sage_core::spectrum::RawSpectrum;
use tokio::io::AsyncReadExt;

#[derive(Debug)]
pub enum Error {
    CloudPathError(crate::Error),
    MzML(crate::mzml::MzMLError),
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
) -> Result<Vec<RawSpectrum>, Error> {
    let res = read_and_execute(s, |bf| async move {
        Ok(crate::mzml::MzMLReader::default()
            .set_signal_to_noise(signal_to_noise)
            .parse(bf)
            .await)
    });

    match res {
        Ok(t) => t.map_err(Error::MzML),
        Err(e) => Err(Error::CloudPathError(e)),
    }
}

pub fn read_fasta<S>(
    path: S,
    decoy_tag: S,
    generate_decoys: bool,
) -> Result<sage_core::fasta::Fasta, Error>
where
    S: AsRef<str>,
{
    let res = read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        Ok(sage_core::fasta::Fasta::parse(
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
            .map_err(crate::Error::IO)?;
        Ok(serde_json::from_str(&contents))
    });

    match res {
        Ok(t) => t.map_err(Error::Json),
        Err(e) => Err(Error::CloudPathError(e)),
    }
}
