//! Parquet serialization for IndexedDatabase
//!
//! Provides functions to serialize and deserialize the fragment ion index
//! to parquet format for fast loading, as well as user-friendly export.

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

/// Convert parquet error to io error
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
                write!(
                    f,
                    "peptide count mismatch: expected {expected}, got {actual}"
                )
            }
            Self::FragmentCountMismatch { expected, actual } => {
                write!(
                    f,
                    "fragment count mismatch: expected {expected}, got {actual}"
                )
            }
            Self::PeptideMismatch { index, field } => {
                write!(f, "peptide mismatch at index {index}: field '{field}'")
            }
            Self::FragmentMismatch { index, field } => {
                write!(f, "fragment mismatch at index {index}: field '{field}'")
            }
            Self::MetadataMismatch { field } => write!(f, "metadata mismatch: field '{field}'"),
        }
    }
}

impl std::error::Error for ValidationError {}

// Schema definitions
fn build_peptides_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        r#"message peptides_schema {
            required int32 peptide_id; required boolean decoy; required byte_array sequence;
            required byte_array modifications; optional float nterm; optional float cterm;
            required float monoisotopic; required int32 missed_cleavages;
            required boolean semi_enzymatic; required int32 position;
            required byte_array proteins (utf8);
        }"#,
    )
}

fn build_fragments_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        r#"message fragments_schema { required int32 peptide_index; required float fragment_mz; }"#,
    )
}

fn build_metadata_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        r#"message metadata_schema { required byte_array key (utf8); required byte_array value (utf8); }"#,
    )
}

fn build_export_schema() -> parquet::errors::Result<Type> {
    parquet::schema::parser::parse_message_type(
        r#"message peptide_export_schema {
            required byte_array peptide (utf8); required byte_array stripped_sequence (utf8);
            required byte_array proteins (utf8); required boolean is_decoy;
            required float monoisotopic_mass; required int32 missed_cleavages;
            required boolean semi_enzymatic; required byte_array position (utf8);
            required int32 num_proteins; required int32 fragment_count;
        }"#,
    )
}

// Conversion helpers
fn position_to_i32(pos: Position) -> i32 {
    match pos {
        Position::Nterm => 0,
        Position::Cterm => 1,
        Position::Full => 2,
        Position::Internal => 3,
    }
}

fn i32_to_position(val: i32) -> Position {
    match val {
        0 => Position::Nterm,
        1 => Position::Cterm,
        2 => Position::Full,
        _ => Position::Internal,
    }
}

fn position_to_str(pos: Position) -> &'static str {
    match pos {
        Position::Nterm => "N-term",
        Position::Cterm => "C-term",
        Position::Full => "Full",
        Position::Internal => "Internal",
    }
}

fn kind_to_str(kind: Kind) -> &'static str {
    match kind {
        Kind::A => "a",
        Kind::B => "b",
        Kind::C => "c",
        Kind::X => "x",
        Kind::Y => "y",
        Kind::Z => "z",
    }
}

fn str_to_kind(s: &str) -> Kind {
    match s {
        "a" => Kind::A,
        "b" => Kind::B,
        "c" => Kind::C,
        "x" => Kind::X,
        "y" => Kind::Y,
        "z" => Kind::Z,
        _ => Kind::B,
    }
}

fn pack_mods(mods: &[f32]) -> Vec<u8> {
    mods.iter().flat_map(|m| m.to_le_bytes()).collect()
}

fn unpack_mods(bytes: &[u8]) -> Vec<f32> {
    bytes
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect()
}

fn serialize_ion_kinds(kinds: &[Kind]) -> String {
    kinds
        .iter()
        .map(|k| kind_to_str(*k))
        .collect::<Vec<_>>()
        .join(",")
}

fn deserialize_ion_kinds(s: &str) -> Vec<Kind> {
    s.split(',')
        .filter(|s| !s.is_empty())
        .map(str_to_kind)
        .collect()
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
        .map(|pairs| {
            pairs
                .into_iter()
                .filter_map(|(spec, mass)| spec.parse().ok().map(|s| (s, mass)))
                .collect()
        })
        .unwrap_or_default()
}

fn serialize_min_values(values: &[f32]) -> String {
    values
        .iter()
        .map(|v| format!("{:08x}", v.to_bits()))
        .collect::<Vec<_>>()
        .join(",")
}

fn deserialize_min_values(s: &str) -> Vec<f32> {
    s.split(',')
        .filter(|s| !s.is_empty())
        .filter_map(|s| u32::from_str_radix(s, 16).ok().map(f32::from_bits))
        .collect()
}

fn writer_props() -> parquet::errors::Result<WriterProperties> {
    Ok(WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(
            ZSTD_COMPRESSION_LEVEL,
        )?))
        .build())
}

// Serialization
fn write_peptides(
    peptides: &[Peptide],
    decoy_tag: &str,
    gen_decoys: bool,
) -> parquet::errors::Result<Vec<u8>> {
    let mut writer = SerializedFileWriter::new(
        Vec::new(),
        build_peptides_schema()?.into(),
        writer_props()?.into(),
    )?;

    for chunk in peptides.chunks(ROW_GROUP_SIZE) {
        let mut rg = writer.next_row_group()?;
        macro_rules! col {
            ($vals:expr, $ty:ty) => {{
                if let Some(mut c) = rg.next_column()? {
                    c.typed::<$ty>().write_batch(&$vals, None, None)?;
                    c.close()?;
                }
            }};
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

        // Optional nterm/cterm
        for getter in [|p: &Peptide| p.nterm, |p: &Peptide| p.cterm] {
            if let Some(mut c) = rg.next_column()? {
                let (mut vals, mut defs) = (Vec::new(), Vec::new());
                for p in chunk {
                    match getter(p) {
                        Some(v) => {
                            vals.push(v);
                            defs.push(1);
                        }
                        None => defs.push(0),
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
                .map(|p| ByteArray::from(p.proteins(decoy_tag, gen_decoys).as_bytes()))
                .collect::<Vec<_>>(),
            ByteArrayType
        );
        rg.close()?;
    }
    writer.into_inner()
}

fn write_fragments(fragments: &[Theoretical]) -> parquet::errors::Result<Vec<u8>> {
    let mut writer = SerializedFileWriter::new(
        Vec::new(),
        build_fragments_schema()?.into(),
        writer_props()?.into(),
    )?;
    for chunk in fragments.chunks(ROW_GROUP_SIZE) {
        let mut rg = writer.next_row_group()?;
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
    writer.into_inner()
}

fn write_metadata(db: &IndexedDatabase) -> parquet::errors::Result<Vec<u8>> {
    let mut writer = SerializedFileWriter::new(
        Vec::new(),
        build_metadata_schema()?.into(),
        writer_props()?.into(),
    )?;
    let meta = [
        ("version", env!("CARGO_PKG_VERSION").to_string()),
        ("bucket_size", db.bucket_size.to_string()),
        ("generate_decoys", db.generate_decoys.to_string()),
        ("decoy_tag", db.decoy_tag.clone()),
        ("ion_kinds", serialize_ion_kinds(&db.ion_kinds)),
        (
            "potential_mods",
            serialize_potential_mods(&db.potential_mods),
        ),
        ("min_values", serialize_min_values(&db.min_value)),
    ];
    let mut rg = writer.next_row_group()?;
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
    writer.into_inner()
}

// Deserialization
fn read_metadata(
    bytes: &[u8],
) -> parquet::errors::Result<std::collections::HashMap<String, String>> {
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes.to_vec()))?;
    let mut meta = std::collections::HashMap::new();
    for row in reader.get_row_iter(None)? {
        let row = row?;
        meta.insert(
            row.get_string(0)?.to_string(),
            row.get_string(1)?.to_string(),
        );
    }
    Ok(meta)
}

fn read_peptides(
    bytes: &[u8],
    decoy_tag: &str,
    gen_decoys: bool,
) -> parquet::errors::Result<Vec<Peptide>> {
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes.to_vec()))?;
    let mut peptides = Vec::new();
    for row in reader.get_row_iter(None)? {
        let row = row?;
        let decoy = row.get_bool(1)?;
        let proteins: Vec<Arc<str>> = row
            .get_string(10)?
            .split(';')
            .filter(|s| !s.is_empty())
            .map(|s| {
                Arc::from(if decoy && gen_decoys && s.starts_with(decoy_tag) {
                    &s[decoy_tag.len()..]
                } else {
                    s
                })
            })
            .collect();
        peptides.push(Peptide {
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
    Ok(peptides)
}

fn read_fragments(bytes: &[u8]) -> parquet::errors::Result<Vec<Theoretical>> {
    let reader = SerializedFileReader::new(bytes::Bytes::from(bytes.to_vec()))?;
    let mut fragments = Vec::new();
    for row in reader.get_row_iter(None)? {
        let row = row?;
        fragments.push(Theoretical {
            peptide_index: PeptideIx(row.get_int(0)? as u32),
            fragment_mz: row.get_float(1)?,
        });
    }
    Ok(fragments)
}

fn read_file_bytes(path: &CloudPath) -> Result<Vec<u8>, Error> {
    match path {
        CloudPath::Local(p) => Ok(std::fs::read(p)?),
        CloudPath::S3 { .. } => Err(Error::IO(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "S3 index loading not yet supported",
        ))),
    }
}

// Public API
pub fn serialize_index(db: &IndexedDatabase, output_dir: &CloudPath) -> Result<(), Error> {
    output_dir.mkdir()?;
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
        let mut path = output_dir.clone();
        path.push(name);
        path.write_bytes_sync(bytes)?;
    }
    Ok(())
}

pub fn deserialize_index(input_dir: &CloudPath) -> Result<IndexedDatabase, Error> {
    let mut meta_path = input_dir.clone();
    meta_path.push("metadata.parquet");
    let meta = read_metadata(&read_file_bytes(&meta_path)?).map_err(pq_err)?;

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
        .map(|s| deserialize_ion_kinds(s))
        .unwrap_or_else(|| vec![Kind::B, Kind::Y]);
    let potential_mods = meta
        .get("potential_mods")
        .map(|s| deserialize_potential_mods(s))
        .unwrap_or_default();
    let min_value = meta
        .get("min_values")
        .map(|s| deserialize_min_values(s))
        .unwrap_or_default();

    let mut pep_path = input_dir.clone();
    pep_path.push("peptides.parquet");
    let peptides =
        read_peptides(&read_file_bytes(&pep_path)?, &decoy_tag, generate_decoys).map_err(pq_err)?;

    let mut frag_path = input_dir.clone();
    frag_path.push("fragments.parquet");
    let fragments = read_fragments(&read_file_bytes(&frag_path)?).map_err(pq_err)?;

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

pub fn validate_index(
    original: &IndexedDatabase,
    loaded: &IndexedDatabase,
) -> Result<(), ValidationError> {
    macro_rules! check {
        ($cond:expr, $err:expr) => {
            if $cond {
                return Err($err);
            }
        };
    }

    check!(
        original.peptides.len() != loaded.peptides.len(),
        ValidationError::PeptideCountMismatch {
            expected: original.peptides.len(),
            actual: loaded.peptides.len()
        }
    );
    check!(
        original.fragments.len() != loaded.fragments.len(),
        ValidationError::FragmentCountMismatch {
            expected: original.fragments.len(),
            actual: loaded.fragments.len()
        }
    );

    check!(
        original.bucket_size != loaded.bucket_size,
        ValidationError::MetadataMismatch {
            field: "bucket_size"
        }
    );
    check!(
        original.generate_decoys != loaded.generate_decoys,
        ValidationError::MetadataMismatch {
            field: "generate_decoys"
        }
    );
    check!(
        original.decoy_tag != loaded.decoy_tag,
        ValidationError::MetadataMismatch { field: "decoy_tag" }
    );
    check!(
        original.ion_kinds != loaded.ion_kinds,
        ValidationError::MetadataMismatch { field: "ion_kinds" }
    );

    check!(
        original.min_value.len() != loaded.min_value.len(),
        ValidationError::MetadataMismatch { field: "min_value" }
    );
    for (a, b) in original.min_value.iter().zip(&loaded.min_value) {
        check!(
            a.to_bits() != b.to_bits(),
            ValidationError::MetadataMismatch { field: "min_value" }
        );
    }

    check!(
        original.potential_mods.len() != loaded.potential_mods.len(),
        ValidationError::MetadataMismatch {
            field: "potential_mods"
        }
    );
    for ((sa, ma), (sb, mb)) in original.potential_mods.iter().zip(&loaded.potential_mods) {
        check!(
            sa != sb || ma.to_bits() != mb.to_bits(),
            ValidationError::MetadataMismatch {
                field: "potential_mods"
            }
        );
    }

    for (i, (o, l)) in original.peptides.iter().zip(&loaded.peptides).enumerate() {
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

    for (i, (o, l)) in original.fragments.iter().zip(&loaded.fragments).enumerate() {
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

pub fn export_index(db: &IndexedDatabase, output_path: &CloudPath) -> Result<(), Error> {
    let mut frag_counts = vec![0i32; db.peptides.len()];
    for f in &db.fragments {
        if (f.peptide_index.0 as usize) < frag_counts.len() {
            frag_counts[f.peptide_index.0 as usize] += 1;
        }
    }

    let mut writer = SerializedFileWriter::new(
        Vec::new(),
        build_export_schema().map_err(pq_err)?.into(),
        writer_props().map_err(pq_err)?.into(),
    )
    .map_err(pq_err)?;

    for (chunk_idx, chunk) in db.peptides.chunks(ROW_GROUP_SIZE).enumerate() {
        let start = chunk_idx * ROW_GROUP_SIZE;
        let mut rg = writer.next_row_group().map_err(pq_err)?;
        macro_rules! col {
            ($vals:expr, $ty:ty) => {{
                if let Some(mut c) = rg.next_column().map_err(pq_err)? {
                    c.typed::<$ty>()
                        .write_batch(&$vals, None, None)
                        .map_err(pq_err)?;
                    c.close().map_err(pq_err)?;
                }
            }};
        }

        col!(
            chunk
                .iter()
                .map(|p| ByteArray::from(p.to_string().as_bytes()))
                .collect::<Vec<_>>(),
            ByteArrayType
        );
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
                .map(|p| ByteArray::from(p.proteins(&db.decoy_tag, db.generate_decoys).as_bytes()))
                .collect::<Vec<_>>(),
            ByteArrayType
        );
        col!(chunk.iter().map(|p| p.decoy).collect::<Vec<_>>(), BoolType);
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
                .map(|p| ByteArray::from(position_to_str(p.position).as_bytes()))
                .collect::<Vec<_>>(),
            ByteArrayType
        );
        col!(
            chunk
                .iter()
                .map(|p| p.proteins.len() as i32)
                .collect::<Vec<_>>(),
            Int32Type
        );
        col!(
            (0..chunk.len())
                .map(|i| frag_counts.get(start + i).copied().unwrap_or(0))
                .collect::<Vec<_>>(),
            Int32Type
        );
        rg.close().map_err(pq_err)?;
    }
    output_path.write_bytes_sync(writer.into_inner().map_err(pq_err)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn make_test_peptide(decoy: bool, seq: &[u8], proteins: Vec<&str>) -> Peptide {
        Peptide {
            decoy,
            sequence: Arc::from(seq),
            modifications: vec![0.0; seq.len()],
            nterm: None,
            cterm: None,
            monoisotopic: 1000.5,
            missed_cleavages: 1,
            semi_enzymatic: false,
            position: Position::Internal,
            proteins: proteins.into_iter().map(Arc::from).collect(),
        }
    }

    fn make_test_database() -> IndexedDatabase {
        let peptides = vec![
            make_test_peptide(false, b"PEPTIDE", vec!["sp|P12345|PROT1"]),
            make_test_peptide(
                false,
                b"SEQVENCE",
                vec!["sp|P12345|PROT1", "sp|P67890|PROT2"],
            ),
            make_test_peptide(true, b"EDITPEP", vec!["sp|P12345|PROT1"]),
        ];
        let fragments = vec![
            Theoretical {
                peptide_index: PeptideIx(0),
                fragment_mz: 500.25,
            },
            Theoretical {
                peptide_index: PeptideIx(0),
                fragment_mz: 600.30,
            },
            Theoretical {
                peptide_index: PeptideIx(1),
                fragment_mz: 450.15,
            },
            Theoretical {
                peptide_index: PeptideIx(2),
                fragment_mz: 550.20,
            },
        ];
        IndexedDatabase {
            peptides,
            fragments,
            ion_kinds: vec![Kind::B, Kind::Y],
            min_value: vec![100.0, 200.0, 300.0],
            potential_mods: vec![
                (ModificationSpecificity::Residue(b'M'), 15.994915),
                (ModificationSpecificity::Residue(b'C'), 57.021464),
            ],
            bucket_size: 8192,
            generate_decoys: true,
            decoy_tag: "rev_".to_string(),
        }
    }

    #[test]
    fn test_pack_unpack_modifications() {
        let mods = vec![1.5f32, -2.3f32, 0.0f32, 100.123f32];
        let unpacked = unpack_mods(&pack_mods(&mods));
        assert!(mods
            .iter()
            .zip(&unpacked)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[test]
    fn test_position_conversion() {
        for (i, pos) in [
            Position::Nterm,
            Position::Cterm,
            Position::Full,
            Position::Internal,
        ]
        .iter()
        .enumerate()
        {
            assert_eq!(position_to_i32(*pos), i as i32);
            assert_eq!(i32_to_position(i as i32), *pos);
        }
        assert_eq!(i32_to_position(99), Position::Internal);
    }

    #[test]
    fn test_ion_kinds_serialization() {
        let kinds = vec![Kind::A, Kind::B, Kind::Y];
        assert_eq!(deserialize_ion_kinds(&serialize_ion_kinds(&kinds)), kinds);
    }

    #[test]
    fn test_min_values_serialization() {
        let values = vec![
            1.0f32,
            2.5f32,
            100.123f32,
            f32::MIN_POSITIVE,
            f32::MAX,
            -0.0f32,
        ];
        let deser = deserialize_min_values(&serialize_min_values(&values));
        assert!(values
            .iter()
            .zip(&deser)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[test]
    fn test_potential_mods_serialization() {
        let mods = vec![
            (ModificationSpecificity::Residue(b'M'), 15.994915f32),
            (ModificationSpecificity::PeptideN(None), 42.010565f32),
            (ModificationSpecificity::Residue(b'C'), 57.021464f32),
        ];
        let deser = deserialize_potential_mods(&serialize_potential_mods(&mods));
        assert!(mods
            .iter()
            .zip(&deser)
            .all(|((sa, ma), (sb, mb))| sa == sb && ma.to_bits() == mb.to_bits()));
    }

    #[test]
    fn test_empty_serialization() {
        assert!(deserialize_min_values(&serialize_min_values(&[])).is_empty());
        assert!(deserialize_potential_mods(&serialize_potential_mods(&[])).is_empty());
    }

    #[test]
    fn test_full_round_trip() {
        let db = make_test_database();
        let tmp = TempDir::new().unwrap();
        let path = CloudPath::Local(tmp.path().to_path_buf());

        serialize_index(&db, &path).expect("serialize failed");
        let loaded = deserialize_index(&path).expect("deserialize failed");
        validate_index(&db, &loaded).expect("validation failed");
    }

    #[test]
    fn test_round_trip_with_optional_mods() {
        let mut db = make_test_database();
        db.peptides[0].nterm = Some(42.01);
        db.peptides[1].cterm = Some(-17.03);
        db.peptides[2].nterm = Some(28.0);
        db.peptides[2].cterm = Some(14.5);

        let tmp = TempDir::new().unwrap();
        let path = CloudPath::Local(tmp.path().to_path_buf());

        serialize_index(&db, &path).unwrap();
        let loaded = deserialize_index(&path).unwrap();
        validate_index(&db, &loaded).expect("validation failed with optional mods");
    }

    #[test]
    fn test_round_trip_preserves_float_precision() {
        let mut db = make_test_database();
        db.peptides[0].monoisotopic = std::f32::consts::PI;
        db.peptides[1].monoisotopic = f32::MIN_POSITIVE;
        db.fragments[0].fragment_mz = std::f32::consts::E;

        let tmp = TempDir::new().unwrap();
        let path = CloudPath::Local(tmp.path().to_path_buf());

        serialize_index(&db, &path).unwrap();
        let loaded = deserialize_index(&path).unwrap();

        assert_eq!(
            db.peptides[0].monoisotopic.to_bits(),
            loaded.peptides[0].monoisotopic.to_bits()
        );
        assert_eq!(
            db.peptides[1].monoisotopic.to_bits(),
            loaded.peptides[1].monoisotopic.to_bits()
        );
        assert_eq!(
            db.fragments[0].fragment_mz.to_bits(),
            loaded.fragments[0].fragment_mz.to_bits()
        );
    }

    #[test]
    fn test_empty_database_round_trip() {
        let db = IndexedDatabase::default();
        let tmp = TempDir::new().unwrap();
        let path = CloudPath::Local(tmp.path().to_path_buf());

        serialize_index(&db, &path).unwrap();
        let loaded = deserialize_index(&path).unwrap();
        validate_index(&db, &loaded).expect("empty database validation failed");
    }

    #[test]
    fn test_validation_detects_peptide_count_mismatch() {
        let db1 = make_test_database();
        let mut db2 = make_test_database();
        db2.peptides.pop();

        let err = validate_index(&db1, &db2).unwrap_err();
        assert!(matches!(err, ValidationError::PeptideCountMismatch { .. }));
    }

    #[test]
    fn test_validation_detects_fragment_mismatch() {
        let db1 = make_test_database();
        let mut db2 = make_test_database();
        db2.fragments[0].fragment_mz = 999.99;

        let err = validate_index(&db1, &db2).unwrap_err();
        assert!(matches!(
            err,
            ValidationError::FragmentMismatch {
                field: "fragment_mz",
                ..
            }
        ));
    }

    #[test]
    fn test_validation_detects_metadata_mismatch() {
        let db1 = make_test_database();
        let mut db2 = make_test_database();
        db2.bucket_size = 9999;

        let err = validate_index(&db1, &db2).unwrap_err();
        assert!(matches!(
            err,
            ValidationError::MetadataMismatch {
                field: "bucket_size"
            }
        ));
    }

    #[test]
    fn test_export_index() {
        let db = make_test_database();
        let tmp = TempDir::new().unwrap();
        let export_path = CloudPath::Local(tmp.path().join("export.parquet"));

        export_index(&db, &export_path).expect("export failed");

        // Verify file exists and is readable
        let bytes = std::fs::read(tmp.path().join("export.parquet")).unwrap();
        assert!(!bytes.is_empty());
    }

    #[test]
    fn test_validation_error_display() {
        let err = ValidationError::PeptideCountMismatch {
            expected: 100,
            actual: 50,
        };
        assert_eq!(
            err.to_string(),
            "peptide count mismatch: expected 100, got 50"
        );

        let err = ValidationError::PeptideMismatch {
            index: 42,
            field: "sequence",
        };
        assert_eq!(
            err.to_string(),
            "peptide mismatch at index 42: field 'sequence'"
        );

        let err = ValidationError::MetadataMismatch { field: "decoy_tag" };
        assert_eq!(err.to_string(), "metadata mismatch: field 'decoy_tag'");
    }
}
