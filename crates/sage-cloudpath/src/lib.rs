use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::write::GzipEncoder;
use futures::TryStreamExt;
use tokio::io::{AsyncBufRead, AsyncRead, AsyncWriteExt, BufReader};

pub use url::Url;

pub mod mgf;
pub mod mzml;
pub mod tdf;
pub mod util;
pub use util::FileFormat;

#[cfg(feature = "parquet")]
pub mod parquet;

/// Convert a path string (local path or cloud URL) into a [`Url`].
pub fn to_url(s: &str) -> Result<Url, Error> {
    if let Ok(url) = Url::parse(s) {
        return Ok(url);
    }
    let path = std::path::Path::new(s);
    let canonical = path.canonicalize()?;
    Url::from_file_path(&canonical).map_err(|_| Error::InvalidUri)
}

/// Does the URL path end in "gz" or "gzip"?
fn gzip_heuristic(url: &Url) -> bool {
    let p = url.path();
    p.ends_with("gz") || p.ends_with("gzip")
}

/// Return the filename portion of a URL path. If the filename ends with `.tdf`,
/// return the parent directory name instead (Bruker `.d` convention).
pub fn filename(url: &Url) -> Option<&str> {
    let path = url.path();
    let name = path.rsplit('/').next().filter(|s| !s.is_empty());
    match name {
        Some(n) if n.ends_with("tdf") => {
            let mut iter = path.rsplit('/');
            iter.next();
            iter.next().filter(|s| !s.is_empty())
        }
        other => other,
    }
}

/// Open a streaming reader for the given URL.
async fn read_url(url: &Url) -> Result<Box<dyn AsyncBufRead + Unpin + Send>, Error> {
    let (store, obj_path) = object_store::parse_url(url).map_err(Error::ObjectStore)?;
    let result = store.get(&obj_path).await.map_err(Error::ObjectStore)?;
    let stream = result
        .into_stream()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e));
    let reader: BufReader<Box<dyn AsyncRead + Unpin + Send>> =
        BufReader::new(Box::new(tokio_util::io::StreamReader::new(stream)));

    if gzip_heuristic(url) {
        Ok(Box::new(BufReader::new(GzipDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

/// Write bytes to the given URL. Gzip-compresses if the path ends in `.gz`.
pub async fn write_bytes_async(url: &Url, bytes: Vec<u8>) -> Result<(), Error> {
    let bytes: Vec<u8> = if gzip_heuristic(url) {
        let inner = Vec::with_capacity(bytes.len() / 2);
        let mut wtr = GzipEncoder::new(inner);
        wtr.write_all(&bytes).await?;
        wtr.flush().await?;
        wtr.into_inner()
    } else {
        bytes
    };

    // Ensure parent directories exist for local paths
    if let Ok(local) = url.to_file_path() {
        if let Some(parent) = local.parent() {
            std::fs::create_dir_all(parent)?;
        }
    }

    let (store, obj_path) = object_store::parse_url(url).map_err(Error::ObjectStore)?;
    store
        .put(&obj_path, bytes::Bytes::from(bytes).into())
        .await
        .map_err(Error::ObjectStore)?;
    Ok(())
}

/// Synchronous wrapper around [`write_bytes_async`].
pub fn write_bytes_sync(url: &Url, bytes: Vec<u8>) -> Result<(), Error> {
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;
    rt.block_on(write_bytes_async(url, bytes))
}

/// Open a reader for `path`, then execute `func` on it. Handles both local and
/// remote paths, with transparent gzip decompression.
pub fn read_and_execute<P, F, Fut, T>(path: P, func: F) -> Result<T, Error>
where
    P: AsRef<str>,
    Fut: futures::Future<Output = Result<T, Error>>,
    F: FnOnce(Box<dyn AsyncBufRead + Unpin>) -> Fut,
{
    let url = to_url(path.as_ref())?;

    // Reject remote URLs that have no object key
    if url.scheme() != "file" {
        let key = url.path().strip_prefix('/').unwrap_or(url.path());
        if key.is_empty() {
            return Err(Error::InvalidUri);
        }
    }

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    rt.block_on(async {
        let reader = read_url(&url).await?;
        func(reader).await
    })
}

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("invalid uri")]
    InvalidUri,
    #[error("object store error: {0}")]
    ObjectStore(#[from] object_store::Error),
    #[error(transparent)]
    IO(#[from] tokio::io::Error),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
    #[error("MzML error: {0}")]
    MzML(#[from] mzml::MzMLError),
    #[error("TDF error: {0}")]
    TDF(#[from] timsrust::TimsRustError),
    #[error("MGF error: {0}")]
    MGF(#[from] mgf::MgfError),
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn filename_gcs() {
        let url = Url::parse("gs://my-bucket/path/to/file.mzML").unwrap();
        assert_eq!(filename(&url), Some("file.mzML"));
    }

    #[test]
    fn filename_azure() {
        let url = Url::parse("az://my-container/path/to/file.mzML").unwrap();
        assert_eq!(filename(&url), Some("file.mzML"));
    }

    #[test]
    fn invalid_remote_read() {
        assert!(read_and_execute("s3://my-bucket", |_| async move { Ok(()) }).is_err())
    }

    #[test]
    fn bruker_filenames() {
        let url = Url::parse("file:///data/20251005_sample_a.d/analysis.tdf").unwrap();
        assert_eq!(filename(&url), Some("20251005_sample_a.d"));

        let url = Url::parse("s3://bucket/baz/20251005_sample_a.d/analysis.tdf").unwrap();
        assert_eq!(filename(&url), Some("20251005_sample_a.d"));

        let url = Url::parse("file:///data/baz/20251005_sample_a.mzML").unwrap();
        assert_eq!(filename(&url), Some("20251005_sample_a.mzML"));
    }

    #[test]
    fn gzip_detection() {
        assert!(gzip_heuristic(&Url::parse("file:///file.mzML.gz").unwrap()));
        assert!(gzip_heuristic(
            &Url::parse("s3://bucket/file.mzML.gzip").unwrap()
        ));
        assert!(!gzip_heuristic(&Url::parse("file:///file.mzML").unwrap()));
    }
}
