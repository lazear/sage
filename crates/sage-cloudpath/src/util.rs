use crate::{read_and_execute, Error};
use sage_core::spectrum::RawSpectrum;
use serde::Serialize;
use tokio::io::AsyncReadExt;

pub fn read_mzml<S: AsRef<str>>(
    s: S,
    file_id: usize,
    signal_to_noise: Option<u8>,
) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(s, |bf| async move {
        Ok(crate::mzml::MzMLReader::with_file_id(file_id)
            .set_signal_to_noise(signal_to_noise)
            .parse(bf)
            .await?)
    })
}

pub fn read_tdf<S: AsRef<str>>(s: S, file_id: usize) -> Result<Vec<RawSpectrum>, Error> {
    let res = crate::tdf::TdfReader.parse(s, file_id);
    match res {
        Ok(t) => Ok(t),
        Err(e) => Err(Error::TDF(e)),
    }
}

#[cfg(feature = "mzdata")]
pub fn read_mzmlb<S: AsRef<str>>(
    s: S,
    file_id: usize,
    signal_to_noise: Option<u8>,
) -> Result<Vec<RawSpectrum>, Error> {
    let res = crate::mzmlb::MzMLbReader::with_file_id(file_id)
        .set_signal_to_noise(signal_to_noise)
        .parse(s.as_ref());
    match res {
        Ok(spectra) => Ok(spectra),
        Err(e) => Err(Error::IO(e)),
    }
}

pub fn read_mgf<S: AsRef<str>>(path: S, file_id: usize) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        let res = crate::mgf::MgfReader::with_file_id(file_id).parse(contents);
        match res {
            Ok(m) => Ok(m),
            Err(e) => Err(Error::MGF(e)),
        }
    })
}

pub fn read_fasta<S>(
    path: S,
    decoy_tag: S,
    generate_decoys: bool,
) -> Result<sage_core::fasta::Fasta, Error>
where
    S: AsRef<str>,
{
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        Ok(sage_core::fasta::Fasta::parse(
            contents,
            decoy_tag.as_ref(),
            generate_decoys,
        ))
    })
}

pub fn read_json<S, T>(path: S) -> Result<T, Error>
where
    S: AsRef<str>,
    T: for<'de> serde::Deserialize<'de>,
{
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents).await?;
        Ok(serde_json::from_str(&contents)?)
    })
}

/// Send telemetry data
pub fn send_data<T>(url: &str, data: &T) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
where
    T: Serialize,
{
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    rt.block_on(async {
        let client = reqwest::ClientBuilder::default().https_only(true).build()?;
        let res = client.post(url).json(data).send().await?;
        res.error_for_status()?;
        Ok(())
    })
}
