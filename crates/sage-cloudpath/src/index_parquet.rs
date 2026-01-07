//! Parquet serialization for IndexedDatabase

#![cfg(feature = "parquet")]

use std::sync::Arc;

use parquet::basic::ZstdLevel;
use parquet::data_type::{BoolType, ByteArray, ByteArrayType, FloatType, Int32Type};
use parquet::file::properties::WriterProperties;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::file::writer::SerializedFileWriter;
use parquet::record::RowAccessor;
use parquet::schema::types::Type;

use sage_core::database::{IndexedDatabase, PeptideIx, Theoretical};
use sage_core::enzyme::Position;
use sage_core::ion_series::Kind;
use sage_core::modification::ModificationSpecificity;
use sage_core::peptide::Peptide;

use crate::parquet::{ROW_GROUP_SIZE, ZSTD_COMPRESSION_LEVEL};
use crate::{CloudPath, Error};

fn pq_err(e: parquet::errors::ParquetError) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::Other, e.to_string())
}

#[derive(Debug)]
pub enum ValidationError {
    PeptideCountMismatch { expected: usize, actual: usize },
    FragmentCountMismatch { expected: usize, actual: usize },
    PeptideMismatch { index: usize, field: &'static str },
    FragmentMismatch { index: usize, field: &'static str },
    MetadataMismatch { field: &'static str },
}

impl std::fmt::Display for ValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PeptideCountMismatch { expected, actual } => {
                write!(f, "peptide count: expected {expected}, got {actual}")
            }
            Self::FragmentCountMismatch { expected, actual } => {
                write!(f, "fragment count: expected {expected}, got {actual}")
            }
            Self::PeptideMismatch { index, field } => {
                write!(f, "peptide[{index}].{field} mismatch")
            }
            Self::FragmentMismatch { index, field } => {
                write!(f, "fragment[{index}].{field} mismatch")
            }
            Self::MetadataMismatch { field } => write!(f, "metadata.{field} mismatch"),
        }
    }
}

impl std::error::Error for ValidationError {}

fn build_peptides_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        "message peptides {
            required int32 id; required boolean decoy; required byte_array seq;
            required byte_array mods; optional float nterm; optional float cterm;
            required float mono; required int32 mc; required boolean semi;
            required int32 pos; required byte_array proteins (utf8);
        }",
    )
}

fn build_fragments_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        "message fragments { required int32 pep_ix; required float mz; }",
    )
}

fn build_metadata_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        "message metadata { required byte_array key (utf8); required byte_array val (utf8); }",
    )
}

fn position_to_i32(p: Position) -> i32 {
    match p {
        Position::Nterm => 0,
        Position::Cterm => 1,
        Position::Full => 2,
        Position::Internal => 3,
    }
}

fn i32_to_position(v: i32) -> Position {
    match v {
        0 => Position::Nterm,
        1 => Position::Cterm,
        2 => Position::Full,
        _ => Position::Internal,
    }
}

fn kind_to_char(k: Kind) -> char {
    match k {
        Kind::A => 'a',
        Kind::B => 'b',
        Kind::C => 'c',
        Kind::X => 'x',
        Kind::Y => 'y',
        Kind::Z => 'z',
    }
}

fn char_to_kind(c: char) -> Kind {
    match c {
        'a' => Kind::A,
        'b' => Kind::B,
        'c' => Kind::C,
        'x' => Kind::X,
        'y' => Kind::Y,
        'z' => Kind::Z,
        _ => Kind::B,
    }
}

fn pack_mods(mods: &[f32]) -> Vec<u8> {
    mods.iter().flat_map(|m| m.to_le_bytes()).collect()
}

fn unpack_mods(b: &[u8]) -> Vec<f32> {
    b.chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect()
}

fn serialize_kinds(kinds: &[Kind]) -> String {
    kinds.iter().map(|k| kind_to_char(*k)).collect()
}

fn deserialize_kinds(s: &str) -> Vec<Kind> {
    s.chars().map(char_to_kind).collect()
}

fn serialize_potential_mods(mods: &[(ModificationSpecificity, f32)]) -> String {
    serde_json::to_string(
        &mods
            .iter()
            .map(|(s, m)| (format!("{s}"), *m))
            .collect::<Vec<_>>(),
    )
    .unwrap_or_default()
}

fn deserialize_potential_mods(s: &str) -> Vec<(ModificationSpecificity, f32)> {
    if s.is_empty() {
        return Vec::new();
    }
    serde_json::from_str::<Vec<(String, f32)>>(s)
        .map(|v| {
            v.into_iter()
                .filter_map(|(spec, m)| spec.parse().ok().map(|s| (s, m)))
                .collect()
        })
        .unwrap_or_default()
}

fn serialize_min_values(values: &[f32]) -> String {
    values
        .iter()
        .map(|v| format!("{:08x}", v.to_bits()))
        .collect::<Vec<_>>()
        .join("")
}

fn deserialize_min_values(s: &str) -> Vec<f32> {
    (0..s.len())
        .step_by(8)
        .filter_map(|i| {
            s.get(i..i + 8)
                .and_then(|h| u32::from_str_radix(h, 16).ok().map(f32::from_bits))
        })
        .collect()
}

fn writer_props() -> parquet::errors::Result<WriterProperties> {
    Ok(WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(
            ZSTD_COMPRESSION_LEVEL,
        )?))
        .build())
}

fn write_peptides(peps: &[Peptide], tag: &str, gen: bool) -> parquet::errors::Result<Vec<u8>> {
    let mut w = SerializedFileWriter::new(
        Vec::new(),
        build_peptides_schema()?.into(),
        writer_props()?.into(),
    )?;
    for chunk in peps.chunks(ROW_GROUP_SIZE) {
        let mut rg = w.next_row_group()?;
        macro_rules! col {
            ($v:expr, $t:ty) => {
                if let Some(mut c) = rg.next_column()? {
                    c.typed::<$t>().write_batch(&$v, None, None)?;
                    c.close()?;
                }
            };
        }
        col!(
            chunk
                .iter()
                .enumerate()
                .map(|(i, _)| i as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        col!(chunk.iter().map(|p| p.decoy).collect::<Vec<_>>(), BoolType);
        col!(
            chunk
                .iter()
                .map(|p| ByteArray::from(p.sequence.as_ref()))
                .collect::<Vec<_>>(),
            ByteArrayType
        );
        col!(
            chunk
                .iter()
                .map(|p| ByteArray::from(pack_mods(&p.modifications)))
                .collect::<Vec<_>>(),
            ByteArrayType
        );
        for getter in [|p: &Peptide| p.nterm, |p: &Peptide| p.cterm] {
            if let Some(mut c) = rg.next_column()? {
                let (mut vals, mut defs) = (Vec::new(), Vec::new());
                for p in chunk {
                    if let Some(v) = getter(p) {
                        vals.push(v);
                        defs.push(1);
                    } else {
                        defs.push(0);
                    }
                }
                c.typed::<FloatType>()
                    .write_batch(&vals, Some(&defs), None)?;
                c.close()?;
            }
        }
        col!(
            chunk.iter().map(|p| p.monoisotopic).collect::<Vec<_>>(),
            FloatType
        );
        col!(
            chunk
                .iter()
                .map(|p| p.missed_cleavages as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        col!(
            chunk.iter().map(|p| p.semi_enzymatic).collect::<Vec<_>>(),
            BoolType
        );
        col!(
            chunk
                .iter()
                .map(|p| position_to_i32(p.position))
                .collect::<Vec<_>>(),
            Int32Type
        );
        col!(
            chunk
                .iter()
                .map(|p| ByteArray::from(p.proteins(tag, gen).as_bytes()))
                .collect::<Vec<_>>(),
            ByteArrayType
        );
        rg.close()?;
    }
    w.into_inner()
}

fn write_fragments(frags: &[Theoretical]) -> parquet::errors::Result<Vec<u8>> {
    let mut w = SerializedFileWriter::new(
        Vec::new(),
        build_fragments_schema()?.into(),
        writer_props()?.into(),
    )?;
    for chunk in frags.chunks(ROW_GROUP_SIZE) {
        let mut rg = w.next_row_group()?;
        if let Some(mut c) = rg.next_column()? {
            c.typed::<Int32Type>().write_batch(
                &chunk
                    .iter()
                    .map(|f| f.peptide_index.0 as i32)
                    .collect::<Vec<_>>(),
                None,
                None,
            )?;
            c.close()?;
        }
        if let Some(mut c) = rg.next_column()? {
            c.typed::<FloatType>().write_batch(
                &chunk.iter().map(|f| f.fragment_mz).collect::<Vec<_>>(),
                None,
                None,
            )?;
            c.close()?;
        }
        rg.close()?;
    }
    w.into_inner()
}

fn write_metadata(db: &IndexedDatabase) -> parquet::errors::Result<Vec<u8>> {
    let mut w = SerializedFileWriter::new(
        Vec::new(),
        build_metadata_schema()?.into(),
        writer_props()?.into(),
    )?;
    let meta = [
        ("bucket_size", db.bucket_size.to_string()),
        ("generate_decoys", db.generate_decoys.to_string()),
        ("decoy_tag", db.decoy_tag.clone()),
        ("ion_kinds", serialize_kinds(&db.ion_kinds)),
        (
            "potential_mods",
            serialize_potential_mods(&db.potential_mods),
        ),
        ("min_values", serialize_min_values(&db.min_value)),
    ];
    let mut rg = w.next_row_group()?;
    if let Some(mut c) = rg.next_column()? {
        c.typed::<ByteArrayType>().write_batch(
            &meta
                .iter()
                .map(|(k, _)| ByteArray::from(k.as_bytes()))
                .collect::<Vec<_>>(),
            None,
            None,
        )?;
        c.close()?;
    }
    if let Some(mut c) = rg.next_column()? {
        c.typed::<ByteArrayType>().write_batch(
            &meta
                .iter()
                .map(|(_, v)| ByteArray::from(v.as_bytes()))
                .collect::<Vec<_>>(),
            None,
            None,
        )?;
        c.close()?;
    }
    rg.close()?;
    w.into_inner()
}

fn read_metadata(
    bytes: &[u8],
) -> parquet::errors::Result<std::collections::HashMap<String, String>> {
    let r = SerializedFileReader::new(bytes::Bytes::from(bytes.to_vec()))?;
    let mut m = std::collections::HashMap::new();
    for row in r.get_row_iter(None)? {
        let row = row?;
        m.insert(row.get_string(0)?.into(), row.get_string(1)?.into());
    }
    Ok(m)
}

fn read_peptides(bytes: &[u8], tag: &str, gen: bool) -> parquet::errors::Result<Vec<Peptide>> {
    let r = SerializedFileReader::new(bytes::Bytes::from(bytes.to_vec()))?;
    let mut peps = Vec::new();
    for row in r.get_row_iter(None)? {
        let row = row?;
        let decoy = row.get_bool(1)?;
        let proteins: Vec<Arc<str>> = row
            .get_string(10)?
            .split(';')
            .filter(|s| !s.is_empty())
            .map(|s| {
                Arc::from(if decoy && gen && s.starts_with(tag) {
                    &s[tag.len()..]
                } else {
                    s
                })
            })
            .collect();
        peps.push(Peptide {
            decoy,
            sequence: Arc::from(row.get_bytes(2)?.data()),
            modifications: unpack_mods(row.get_bytes(3)?.data()),
            nterm: row.get_float(4).ok(),
            cterm: row.get_float(5).ok(),
            monoisotopic: row.get_float(6)?,
            missed_cleavages: row.get_int(7)? as u8,
            semi_enzymatic: row.get_bool(8)?,
            position: i32_to_position(row.get_int(9)?),
            proteins,
        });
    }
    Ok(peps)
}

fn read_fragments(bytes: &[u8]) -> parquet::errors::Result<Vec<Theoretical>> {
    let r = SerializedFileReader::new(bytes::Bytes::from(bytes.to_vec()))?;
    let mut frags = Vec::new();
    for row in r.get_row_iter(None)? {
        let row = row?;
        frags.push(Theoretical {
            peptide_index: PeptideIx(row.get_int(0)? as u32),
            fragment_mz: row.get_float(1)?,
        });
    }
    Ok(frags)
}

fn read_file(path: &CloudPath) -> Result<Vec<u8>, Error> {
    match path {
        CloudPath::Local(p) => Ok(std::fs::read(p)?),
        CloudPath::S3 { .. } => Err(Error::IO(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "S3 not supported",
        ))),
    }
}

/// Serialize an IndexedDatabase to a parquet directory
pub fn serialize_index(db: &IndexedDatabase, dir: &CloudPath) -> Result<(), Error> {
    dir.mkdir()?;
    for (name, bytes) in [
        (
            "peptides.parquet",
            write_peptides(&db.peptides, &db.decoy_tag, db.generate_decoys).map_err(pq_err)?,
        ),
        (
            "fragments.parquet",
            write_fragments(&db.fragments).map_err(pq_err)?,
        ),
        ("metadata.parquet", write_metadata(db).map_err(pq_err)?),
    ] {
        let mut p = dir.clone();
        p.push(name);
        p.write_bytes_sync(bytes)?;
    }
    Ok(())
}

/// Deserialize an IndexedDatabase from a parquet directory
pub fn deserialize_index(dir: &CloudPath) -> Result<IndexedDatabase, Error> {
    let mut mp = dir.clone();
    mp.push("metadata.parquet");
    let meta = read_metadata(&read_file(&mp)?).map_err(pq_err)?;

    let bucket_size = meta
        .get("bucket_size")
        .and_then(|s| s.parse().ok())
        .unwrap_or(8192);
    let generate_decoys = meta
        .get("generate_decoys")
        .map(|s| s == "true")
        .unwrap_or(true);
    let decoy_tag = meta
        .get("decoy_tag")
        .cloned()
        .unwrap_or_else(|| "rev_".into());
    let ion_kinds = meta
        .get("ion_kinds")
        .map(|s| deserialize_kinds(s))
        .unwrap_or_else(|| vec![Kind::B, Kind::Y]);
    let potential_mods = meta
        .get("potential_mods")
        .map(|s| deserialize_potential_mods(s))
        .unwrap_or_default();
    let min_value = meta
        .get("min_values")
        .map(|s| deserialize_min_values(s))
        .unwrap_or_default();

    let mut pp = dir.clone();
    pp.push("peptides.parquet");
    let peptides = read_peptides(&read_file(&pp)?, &decoy_tag, generate_decoys).map_err(pq_err)?;

    let mut fp = dir.clone();
    fp.push("fragments.parquet");
    let fragments = read_fragments(&read_file(&fp)?).map_err(pq_err)?;

    Ok(IndexedDatabase {
        peptides,
        fragments,
        ion_kinds,
        min_value,
        potential_mods,
        bucket_size,
        generate_decoys,
        decoy_tag,
    })
}

/// Validate that two IndexedDatabases are identical (bit-exact for floats)
pub fn validate_index(a: &IndexedDatabase, b: &IndexedDatabase) -> Result<(), ValidationError> {
    macro_rules! check {
        ($c:expr, $e:expr) => {
            if $c {
                return Err($e);
            }
        };
    }

    check!(
        a.peptides.len() != b.peptides.len(),
        ValidationError::PeptideCountMismatch {
            expected: a.peptides.len(),
            actual: b.peptides.len()
        }
    );
    check!(
        a.fragments.len() != b.fragments.len(),
        ValidationError::FragmentCountMismatch {
            expected: a.fragments.len(),
            actual: b.fragments.len()
        }
    );
    check!(
        a.bucket_size != b.bucket_size,
        ValidationError::MetadataMismatch {
            field: "bucket_size"
        }
    );
    check!(
        a.generate_decoys != b.generate_decoys,
        ValidationError::MetadataMismatch {
            field: "generate_decoys"
        }
    );
    check!(
        a.decoy_tag != b.decoy_tag,
        ValidationError::MetadataMismatch { field: "decoy_tag" }
    );
    check!(
        a.ion_kinds != b.ion_kinds,
        ValidationError::MetadataMismatch { field: "ion_kinds" }
    );
    check!(
        a.min_value.len() != b.min_value.len(),
        ValidationError::MetadataMismatch { field: "min_value" }
    );
    for (x, y) in a.min_value.iter().zip(&b.min_value) {
        check!(
            x.to_bits() != y.to_bits(),
            ValidationError::MetadataMismatch { field: "min_value" }
        );
    }
    check!(
        a.potential_mods.len() != b.potential_mods.len(),
        ValidationError::MetadataMismatch {
            field: "potential_mods"
        }
    );
    for ((sa, ma), (sb, mb)) in a.potential_mods.iter().zip(&b.potential_mods) {
        check!(
            sa != sb || ma.to_bits() != mb.to_bits(),
            ValidationError::MetadataMismatch {
                field: "potential_mods"
            }
        );
    }

    for (i, (o, l)) in a.peptides.iter().zip(&b.peptides).enumerate() {
        check!(
            o.decoy != l.decoy,
            ValidationError::PeptideMismatch {
                index: i,
                field: "decoy"
            }
        );
        check!(
            o.sequence != l.sequence,
            ValidationError::PeptideMismatch {
                index: i,
                field: "sequence"
            }
        );
        check!(
            o.modifications != l.modifications,
            ValidationError::PeptideMismatch {
                index: i,
                field: "modifications"
            }
        );
        check!(
            o.nterm != l.nterm,
            ValidationError::PeptideMismatch {
                index: i,
                field: "nterm"
            }
        );
        check!(
            o.cterm != l.cterm,
            ValidationError::PeptideMismatch {
                index: i,
                field: "cterm"
            }
        );
        check!(
            o.monoisotopic.to_bits() != l.monoisotopic.to_bits(),
            ValidationError::PeptideMismatch {
                index: i,
                field: "monoisotopic"
            }
        );
        check!(
            o.missed_cleavages != l.missed_cleavages,
            ValidationError::PeptideMismatch {
                index: i,
                field: "missed_cleavages"
            }
        );
        check!(
            o.semi_enzymatic != l.semi_enzymatic,
            ValidationError::PeptideMismatch {
                index: i,
                field: "semi_enzymatic"
            }
        );
        check!(
            o.position != l.position,
            ValidationError::PeptideMismatch {
                index: i,
                field: "position"
            }
        );
        check!(
            o.proteins != l.proteins,
            ValidationError::PeptideMismatch {
                index: i,
                field: "proteins"
            }
        );
    }

    for (i, (o, l)) in a.fragments.iter().zip(&b.fragments).enumerate() {
        check!(
            o.peptide_index != l.peptide_index,
            ValidationError::FragmentMismatch {
                index: i,
                field: "peptide_index"
            }
        );
        check!(
            o.fragment_mz.to_bits() != l.fragment_mz.to_bits(),
            ValidationError::FragmentMismatch {
                index: i,
                field: "fragment_mz"
            }
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn test_db() -> IndexedDatabase {
        IndexedDatabase {
            peptides: vec![Peptide {
                decoy: false,
                sequence: Arc::from(b"PEPTIDE".as_slice()),
                modifications: vec![0.0; 7],
                nterm: Some(42.01),
                cterm: None,
                monoisotopic: std::f32::consts::PI,
                missed_cleavages: 1,
                semi_enzymatic: false,
                position: Position::Internal,
                proteins: vec![Arc::from("PROT1")],
            }],
            fragments: vec![Theoretical {
                peptide_index: PeptideIx(0),
                fragment_mz: std::f32::consts::E,
            }],
            ion_kinds: vec![Kind::B, Kind::Y],
            min_value: vec![100.0, f32::MIN_POSITIVE],
            potential_mods: vec![(ModificationSpecificity::Residue(b'M'), 15.994915)],
            bucket_size: 8192,
            generate_decoys: true,
            decoy_tag: "rev_".into(),
        }
    }

    #[test]
    fn round_trip() {
        let db = test_db();
        let tmp = TempDir::new().unwrap();
        let path = CloudPath::Local(tmp.path().into());
        serialize_index(&db, &path).unwrap();
        let loaded = deserialize_index(&path).unwrap();
        validate_index(&db, &loaded).unwrap();
    }

    #[test]
    fn empty_db() {
        let db = IndexedDatabase::default();
        let tmp = TempDir::new().unwrap();
        let path = CloudPath::Local(tmp.path().into());
        serialize_index(&db, &path).unwrap();
        let loaded = deserialize_index(&path).unwrap();
        validate_index(&db, &loaded).unwrap();
    }

    #[test]
    fn detects_mismatch() {
        let db1 = test_db();
        let mut db2 = test_db();
        db2.fragments[0].fragment_mz = 999.0;
        assert!(validate_index(&db1, &db2).is_err());
    }
}
