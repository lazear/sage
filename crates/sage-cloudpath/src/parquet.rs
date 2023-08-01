#![cfg(feature = "parquet")]

use parquet::data_type::FloatType;
use parquet::{
    basic::ZstdLevel,
    data_type::{ByteArrayType, DataType, Int32Type},
    file::{properties::WriterProperties, writer::SerializedFileWriter},
    schema::types::Type,
};
use sage_core::database::IndexedDatabase;
use sage_core::scoring::Feature;

pub fn build_schema() -> anyhow::Result<Type> {
    let msg = r#"
        message schema {
            required byte_array filename (utf8);
            required byte_array scannr (utf8);
            required byte_array peptide (utf8);
            required byte_array proteins (utf8);
            required int32 num_proteins;
            required int32 rank;
            required int32 label;
            required float expmass;
            required float calcmass;
            required int32 charge;
            required int32 peptide_len;
            required int32 missed_cleavages;
            required float isotope_error;
            required float precursor_ppm;
            required float fragment_ppm;
            required float hyperscore;
            required float delta_next;
            required float delta_best;
            required float rt;
            required float aligned_rt;
            required float predicted_rt;
            required float delta_rt_model;
            required int32 matched_peaks;
            required int32 longest_b;
            required int32 longest_y;
            required float longest_y_pct;
            required float matched_intensity_pct;
            required int32 scored_candidates;
            required float poisson;
            required float sage_discriminant_score;
            required float posterior_error;
            required float spectrum_q;
            required float peptide_q;
            required float protein_q;
        }
    "#;
    let schema = parquet::schema::parser::parse_message_type(msg)?;
    Ok(schema)
}

pub fn serialize_parquet_to_bytes(
    features: &[Feature],
    filenames: &[String],
    database: &IndexedDatabase,
) -> anyhow::Result<Vec<u8>> {
    let schema = build_schema()?;

    let options = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::ZSTD(ZstdLevel::try_new(4)?))
        .build();

    let buf = Vec::new();
    // let buf = std::fs::File::create("test.parquet")?;
    let mut writer = SerializedFileWriter::new(buf, schema.into(), options.into())?;
    let mut rg = writer.next_row_group()?;

    macro_rules! write_col {
        ($field:ident, $ty:ident) => {
            if let Some(mut col) = rg.next_column()? {
                col.typed::<$ty>().write_batch(
                    &features
                        .iter()
                        .map(|f| f.$field as <$ty as DataType>::T)
                        .collect::<Vec<_>>(),
                    None,
                    None,
                )?;
                col.close()?;
            }
        };
        ($lambda:expr, $ty:ident) => {
            if let Some(mut col) = rg.next_column()? {
                col.typed::<$ty>().write_batch(
                    &features.iter().map($lambda).collect::<Vec<_>>(),
                    None,
                    None,
                )?;
                col.close()?;
            }
        };
    }

    write_col!(
        |f: &Feature| filenames[f.file_id].as_str().into(),
        ByteArrayType
    );
    write_col!(|f: &Feature| f.spec_id.as_str().into(), ByteArrayType);
    write_col!(
        |f: &Feature| database[f.peptide_idx].sequence.as_ref().into(),
        ByteArrayType
    );
    write_col!(
        |f: &Feature| database[f.peptide_idx]
            .proteins(&database.decoy_tag, database.generate_decoys)
            .as_str()
            .into(),
        ByteArrayType
    );
    write_col!(
        |f: &Feature| database[f.peptide_idx].proteins.len() as i32,
        Int32Type
    );
    write_col!(rank, Int32Type);
    write_col!(label, Int32Type);
    write_col!(expmass, FloatType);
    write_col!(calcmass, FloatType);
    write_col!(charge, Int32Type);
    write_col!(peptide_len, Int32Type);
    write_col!(missed_cleavages, Int32Type);
    write_col!(isotope_error, FloatType);
    write_col!(delta_mass, FloatType);
    write_col!(average_ppm, FloatType);
    write_col!(hyperscore, FloatType);
    write_col!(delta_next, FloatType);
    write_col!(delta_best, FloatType);
    write_col!(rt, FloatType);
    write_col!(aligned_rt, FloatType);
    write_col!(predicted_rt, FloatType);
    write_col!(delta_rt_model, FloatType);
    write_col!(matched_peaks, Int32Type);
    write_col!(longest_b, Int32Type);
    write_col!(longest_y, Int32Type);
    write_col!(longest_y_pct, FloatType);
    write_col!(matched_intensity_pct, FloatType);
    write_col!(scored_candidates, Int32Type);
    write_col!(poisson, FloatType);
    write_col!(discriminant_score, FloatType);
    write_col!(posterior_error, FloatType);
    write_col!(spectrum_q, FloatType);
    write_col!(peptide_q, FloatType);
    write_col!(protein_q, FloatType);

    rg.close()?;
    writer.into_inner().map_err(Into::into)
}
