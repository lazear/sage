use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::write::GzipEncoder;
use futures::TryStreamExt;
use object_store::{ObjectStore, ObjectStoreExt};
use tokio::io::{AsyncBufRead, AsyncRead, AsyncWriteExt, BufReader};

pub use url::Url;

pub mod mgf;
pub mod mzml;
pub mod tdf;
pub mod util;
pub use util::FileFormat;

#[cfg(feature = "parquet")]
pub mod parquet;

/// Schemes recognized by `object_store::parse_url_opts`. Anything outside
/// this set — most importantly Windows drive letters like `C:` which parse
/// as single-letter URL schemes — is treated as a local path.
const OBJECT_STORE_SCHEMES: &[&str] = &[
    "file", "memory", "s3", "s3a", "gs", "az", "adl", "azure", "abfs", "abfss", "http", "https",
];

/// Parse `s` as a URL, but only accept schemes that `object_store` knows how
/// to handle. Returns `None` for local paths (including Windows paths like
/// `C:\foo` that would otherwise parse as a URL with scheme `c`).
pub fn try_parse_url(s: &str) -> Option<Url> {
    Url::parse(s)
        .ok()
        .filter(|u| OBJECT_STORE_SCHEMES.contains(&u.scheme()))
}

/// Convert a path string (local path or cloud URL) into a [`Url`].
pub fn to_url(s: &str) -> Result<Url, Error> {
    if let Some(url) = try_parse_url(s) {
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

fn parse_url(url: &Url) -> Result<(Box<dyn ObjectStore>, object_store::path::Path), Error> {
    // AWS and Azure require lowercased config keys. By default, these aren't pulled from the env
    object_store::parse_url_opts(
        url,
        std::env::vars().map(|(k, v)| (k.to_ascii_lowercase(), v)),
    )
    .map_err(Error::ObjectStore)
}

/// Open a streaming reader for the given URL.
async fn read_url(url: &Url) -> Result<Box<dyn AsyncBufRead + Unpin + Send>, Error> {
    let (store, obj_path) = parse_url(url)?;
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

    let (store, obj_path) = parse_url(url)?;
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
    fn windows_drive_letter_is_not_a_url() {
        // `Url::parse("C:\\...")` succeeds with `c` as a single-letter
        // scheme, which object_store later rejects with "Unable to recognise
        // URL". `to_url` must treat such inputs as local paths — here the
        // path doesn't exist on this machine, so we expect an IO error
        // rather than a bogus `Ok(Url { scheme: "c", ... })`.
        let backslash = to_url(r"C:\Users\nonexistent\bar.json");
        assert!(
            matches!(backslash, Err(Error::IO(_))),
            "expected IO error for Windows path with backslashes, got {:?}",
            backslash
        );

        let forwardslash = to_url("C:/Users/nonexistent/bar.json");
        assert!(
            matches!(forwardslash, Err(Error::IO(_))),
            "expected IO error for Windows path with forward slashes, got {:?}",
            forwardslash
        );
    }

    #[test]
    fn cloud_urls_still_parse() {
        assert_eq!(to_url("s3://bucket/key").unwrap().scheme(), "s3");
        assert_eq!(to_url("gs://bucket/key").unwrap().scheme(), "gs");
        assert_eq!(to_url("az://container/key").unwrap().scheme(), "az");
        assert_eq!(to_url("https://example.com/key").unwrap().scheme(), "https");
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
