use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::write::GzipEncoder;
use http::Uri;
use sage_core::mzml::{MzMLError, MzMLReader, Spectrum};
use std::path::PathBuf;
use std::str::FromStr;
use tokio::io::{AsyncBufRead, AsyncRead, AsyncWriteExt, BufReader};

static S3_CLIENT: once_cell::sync::OnceCell<aws_sdk_s3::Client> = once_cell::sync::OnceCell::new();

async fn s3_client() -> &'static aws_sdk_s3::Client {
    if S3_CLIENT.get().is_none() {
        let config = aws_config::load_from_env().await;
        let client = aws_sdk_s3::Client::new(&config);
        S3_CLIENT.get_or_init(|| client)
    } else {
        S3_CLIENT.get().expect("once_cell invariant")
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum CloudPath {
    S3 { bucket: String, key: String },
    Local(PathBuf),
}

impl FromStr for CloudPath {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let uri = match s.parse::<Uri>() {
            Ok(uri) => uri,
            Err(_) => return Ok(Self::Local(s.into())),
        };
        match uri.scheme_str() {
            Some("s3") => {
                let bucket = uri.authority().ok_or(Error::InvalidUri)?.to_string();
                let key = uri
                    .path_and_query()
                    .and_then(|path| path.path().strip_prefix('/'))
                    .ok_or(Error::InvalidUri)?
                    .into();
                Ok(Self::S3 { bucket, key })
            }
            Some(_) => Err(Error::InvalidUri),
            None => Ok(Self::Local(s.into())),
        }
    }
}

impl std::fmt::Display for CloudPath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CloudPath::S3 { bucket, key } => write!(f, "s3://{}/{}", bucket, key),
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
    async fn mk_bufreader(&self) -> Result<BufReader<Box<dyn AsyncRead + Unpin>>, Error> {
        match self {
            Self::S3 { bucket, key } => {
                let object = s3_client()
                    .await
                    .get_object()
                    .bucket(bucket)
                    .key(key)
                    .send()
                    .await
                    .map_err(|e| Error::S3Error(e.into()))?;

                Ok(BufReader::new(Box::new(object.body.into_async_read())))
            }
            Self::Local(path) => Ok(BufReader::new(Box::new(
                tokio::fs::File::open(path).await.map_err(Error::IOError)?,
            ))),
        }
    }

    pub async fn read(&self) -> Result<Box<dyn AsyncBufRead + Unpin>, Error> {
        let reader = self.mk_bufreader().await?;
        match self.gzip_heuristic() {
            true => {
                let gzip = GzipDecoder::new(reader);
                Ok(Box::new(BufReader::new(gzip)))
            }
            false => Ok(Box::new(reader)),
        }
    }

    pub async fn write_bytes(&self, bytes: Vec<u8>) -> Result<(), Error> {
        let bytes: Vec<u8> = match self.gzip_heuristic() {
            true => {
                let inner = Vec::with_capacity(bytes.len() / 2);
                let mut wtr = GzipEncoder::new(inner);
                wtr.write_all(&bytes).await.map_err(Error::IOError)?;
                wtr.flush().await.map_err(Error::IOError)?;
                wtr.into_inner()
            }
            false => bytes,
        };
        match self {
            Self::Local(path) => {
                let mut file = tokio::fs::File::create(path)
                    .await
                    .map_err(Error::IOError)?;
                file.write_all(&bytes).await.map_err(Error::IOError)?;
                Ok(())
            }
            Self::S3 { bucket, key } => Ok(multipart_upload(bucket, key, bytes)
                .await
                .map_err(Error::S3Error)?),
        }
    }

    pub fn write_bytes_sync(&self, bytes: Vec<u8>) -> Result<(), Error> {
        let rt = tokio::runtime::Builder::new_current_thread()
            .enable_all()
            .build()
            .map_err(Error::IOError)?;

        rt.block_on(async { self.write_bytes(bytes).await })
    }
}

/// For very large searches, the size of the results file can exceed the 5GB
/// limit for sending S3 PUT requests. Instead, we can use multipart uploads
async fn multipart_upload(
    bucket: &str,
    key: &str,
    bytes: Vec<u8>,
) -> Result<(), aws_sdk_s3::Error> {
    let upload = s3_client()
        .await
        .create_multipart_upload()
        .bucket(bucket)
        .key(key)
        .send()
        .await?;

    let upload_id = upload.upload_id().ok_or_else(|| {
        aws_sdk_s3::Error::NoSuchUpload(
            aws_sdk_s3::error::no_such_upload::Builder::default()
                .message("S3 CreateMultipartUpload did not return upload id!")
                .build(),
        )
    })?;

    let bytes: bytes::Bytes = bytes.into();

    // 256 MB part size
    let part_size = 256 * 1024 * 1024;
    let mut completed_parts = Vec::new();
    for (part_number, start) in (0..bytes.len()).step_by(part_size).enumerate() {
        let end = bytes.len().min(start + part_size);
        let body = bytes.slice(start..end);

        // Part numbers can be 1 to 10,000
        let part_number = 1 + part_number as i32;

        let resp = s3_client()
            .await
            .upload_part()
            .bucket(bucket)
            .key(key)
            .part_number(part_number)
            .upload_id(upload_id)
            .body(body.into())
            .send()
            .await?;

        completed_parts.push(
            aws_sdk_s3::model::completed_part::Builder::default()
                .part_number(part_number)
                .e_tag(resp.e_tag().unwrap_or_default())
                .build(),
        );
    }

    let completed_upload = aws_sdk_s3::model::completed_multipart_upload::Builder::default()
        .set_parts(Some(completed_parts))
        .build();

    s3_client()
        .await
        .complete_multipart_upload()
        .bucket(bucket)
        .key(key)
        .upload_id(upload_id)
        .multipart_upload(completed_upload)
        .send()
        .await?;
    Ok(())
}

pub fn read_mzml<S: AsRef<str>>(s: S) -> Result<Vec<Spectrum>, Error> {
    let path = s.as_ref().parse::<CloudPath>()?;

    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .map_err(Error::IOError)?;

    rt.block_on(async {
        let reader = path.read().await?;
        MzMLReader::default()
            .parse(reader)
            .await
            .map_err(Error::MzMLError)
    })
}

#[derive(Debug)]
pub enum Error {
    InvalidUri,
    MzMLError(MzMLError),
    S3Error(aws_sdk_s3::Error),
    IOError(tokio::io::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidUri => f.write_str("invalid URI"),
            Error::MzMLError(x) => x.fmt(f),
            Error::S3Error(x) => x.fmt(f),
            Error::IOError(x) => x.fmt(f),
        }
    }
}

impl std::error::Error for Error {}

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
