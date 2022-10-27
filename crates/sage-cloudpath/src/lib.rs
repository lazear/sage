use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::write::GzipEncoder;
use http::Uri;
use sage_core::mzml::{IOError, MzMlReader, Spectrum};
use std::path::PathBuf;
use std::str::FromStr;
use tokio::io::{AsyncBufRead, AsyncRead, AsyncReadExt, AsyncWriteExt, BufReader};

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum CloudPath {
    S3 { bucket: String, key: String },
    Local(PathBuf),
}

impl FromStr for CloudPath {
    type Err = IOError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let uri = match s.parse::<Uri>() {
            Ok(uri) => uri,
            Err(_) => return Ok(Self::Local(s.into())),
        };
        match uri.scheme_str() {
            Some("s3") => {
                let bucket = uri.authority().ok_or(IOError::InvalidUri)?.to_string();
                let key = uri
                    .path_and_query()
                    .and_then(|path| path.path().strip_prefix('/'))
                    .ok_or(IOError::InvalidUri)?
                    .into();
                Ok(Self::S3 { bucket, key })
            }
            Some(_) => Err(IOError::InvalidUri),
            None => Ok(Self::Local(s.into())),
        }
    }
}

impl std::fmt::Display for CloudPath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CloudPath::S3 { bucket, key } => {
                write!(f, "{}/{}", bucket, key)
            }
            CloudPath::Local(path) => write!(f, "{}", path.display()),
        }
    }
}

impl CloudPath {
    /// If `self` is a local path, recursively create directories as needed
    pub fn mkdir(&self) -> std::io::Result<()> {
        match self {
            CloudPath::S3 { .. } => Ok(()),
            CloudPath::Local(path) => std::fs::create_dir_all(path),
        }
    }

    /// Push a filename to the current path
    /// * If `self` is an S3 path, this adds to end of the object key, using '/' as a delimiter
    /// * If `self` is a local path, this behaves the same as [`PathBuf::push`]
    pub fn push<P: AsRef<str>>(&mut self, path: P) {
        match self {
            CloudPath::S3 { key, .. } => {
                if !key.is_empty() {
                    key.push('/');
                }

                key.push_str(path.as_ref());
            }
            CloudPath::Local(key) => key.push(path.as_ref()),
        }
    }

    /// Does the path end in "gz" or "gzip"?
    fn gzip_heuristic(&self) -> bool {
        match &self {
            Self::S3 { key, .. } => key.ends_with("gz") || key.ends_with("gzip"),
            Self::Local(path) => match path.extension() {
                Some(ext) => ext.to_ascii_lowercase() == "gz" || ext.to_ascii_lowercase() == "gzip",
                _ => false,
            },
        }
    }

    /// Open a BufReader stream to the object
    async fn mk_bufreader(&self) -> Result<BufReader<Box<dyn AsyncRead + Unpin>>, IOError> {
        match self {
            Self::S3 { bucket, key } => {
                let config = aws_config::load_from_env().await;
                let client = aws_sdk_s3::Client::new(&config);

                let object = client
                    .get_object()
                    .bucket(bucket)
                    .key(key)
                    .send()
                    .await
                    .map_err(|e| IOError::OtherError(e.into()))?;

                Ok(BufReader::new(Box::new(object.body.into_async_read())))
            }
            Self::Local(path) => Ok(BufReader::new(Box::new(
                tokio::fs::File::open(path)
                    .await
                    .map_err(IOError::IOError)?,
            ))),
        }
    }

    pub async fn read(&self) -> Result<Box<dyn AsyncBufRead + Unpin>, IOError> {
        let reader = self.mk_bufreader().await?;
        match self.gzip_heuristic() {
            true => {
                let gzip = GzipDecoder::new(reader);
                Ok(Box::new(BufReader::new(gzip)))
            }
            false => Ok(Box::new(reader)),
        }
    }

    pub async fn write_bytes(&self, bytes: Vec<u8>) -> Result<(), IOError> {
        let bytes = match self.gzip_heuristic() {
            true => {
                let inner = Vec::with_capacity(bytes.len() / 2);
                let mut wtr = GzipEncoder::new(inner);
                wtr.write_all(&bytes).await.map_err(IOError::IOError)?;
                wtr.flush().await.map_err(IOError::IOError)?;
                wtr.into_inner()
            }
            false => bytes,
        };
        match self {
            Self::Local(path) => {
                let mut file = tokio::fs::File::create(path)
                    .await
                    .map_err(IOError::IOError)?;
                file.write_all(&bytes).await.map_err(IOError::IOError)?;
                Ok(())
            }
            Self::S3 { bucket, key } => {
                let config = aws_config::load_from_env().await;
                let client = aws_sdk_s3::Client::new(&config);
                let bytes: bytes::Bytes = bytes.into();

                client
                    .put_object()
                    .bucket(bucket)
                    .body(bytes.into())
                    .key(key)
                    .send()
                    .await
                    .map_err(|e| IOError::OtherError(e.into()))?;
                Ok(())
            }
        }
    }

    pub fn write_bytes_sync(&self, bytes: Vec<u8>) -> Result<(), IOError> {
        let rt = tokio::runtime::Builder::new_current_thread()
            .enable_all()
            .build()
            .map_err(IOError::IOError)?;

        rt.block_on(async { self.write_bytes(bytes).await })
    }
}

pub fn read_bytes<S: AsRef<str>>(s: S) -> Result<Vec<u8>, IOError> {
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .map_err(IOError::IOError)?;

    let path = s.as_ref().parse::<CloudPath>()?;
    rt.block_on(async {
        let mut buffer = Vec::new();
        path.read()
            .await?
            .read_to_end(&mut buffer)
            .await
            .map_err(IOError::IOError)?;
        Ok(buffer)
    })
}

pub fn read_mzml<S: AsRef<str>>(s: S) -> Result<Vec<Spectrum>, IOError> {
    let buffer = read_bytes(s)?;
    let reader = std::io::BufReader::new(buffer.as_slice());
    MzMlReader::default().parse(reader)
}

#[cfg(test)]
mod test {
    use super::CloudPath;

    #[test]
    fn cloudpath_s3() {
        let a = CloudPath::S3 {
            bucket: "awscdkstack-bucket-foobar".into(),
            key: "/file.mzML.gz".into(),
        };
        assert_eq!(
            "s3://awscdkstack-bucket-foobar//file.mzML.gz"
                .parse::<CloudPath>()
                .unwrap(),
            a
        );

        let b = CloudPath::S3 {
            bucket: "awscdkstack-bucket-foobar".into(),
            key: "prefix/file.mzML.gz".into(),
        };
        assert_eq!(
            "s3://awscdkstack-bucket-foobar/prefix/file.mzML.gz"
                .parse::<CloudPath>()
                .unwrap(),
            b
        );

        let mut c = CloudPath::S3 {
            bucket: "awscdkstack-bucket-foobar".into(),
            key: Default::default(),
        };
        assert_eq!(
            "s3://awscdkstack-bucket-foobar/"
                .parse::<CloudPath>()
                .unwrap(),
            c
        );
        c.push("prefix");
        assert_eq!(
            "s3://awscdkstack-bucket-foobar/prefix"
                .parse::<CloudPath>()
                .unwrap(),
            c
        );
        c.push("test.mzml");
        assert_eq!(
            "s3://awscdkstack-bucket-foobar/prefix/test.mzml"
                .parse::<CloudPath>()
                .unwrap(),
            c
        );
    }
}
