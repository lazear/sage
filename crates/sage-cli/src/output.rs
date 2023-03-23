use std::collections::HashMap;

use rayon::prelude::*;
use sage_core::{database::PeptideIx, scoring::Feature, tmt::TmtQuant};

use crate::Runner;

impl Runner {
    pub fn serialize_feature(&self, feature: &Feature, filenames: &[String]) -> csv::ByteRecord {
        let mut record = csv::ByteRecord::new();
        record.push_field(feature.peptide.as_str().as_bytes());
        record.push_field(feature.proteins.as_str().as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.num_proteins).as_bytes());
        record.push_field(filenames[feature.file_id].as_bytes());
        record.push_field(feature.spec_id.as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.rank).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.label).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.expmass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.calcmass).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.charge).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.peptide_len).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format(feature.missed_cleavages)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.isotope_error).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_mass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.hyperscore).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_hyperscore)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.aligned_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_rt_model).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.matched_peaks).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_b).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_y).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.longest_y_pct).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.matched_intensity_pct)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.scored_candidates)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.poisson).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.discriminant_score)
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.posterior_error)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.spectrum_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.peptide_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.protein_q).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ms1_intensity).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.ms2_intensity).as_bytes());
        record
    }

    pub fn write_features(
        &self,
        features: Vec<Feature>,
        filenames: &[String],
    ) -> anyhow::Result<String> {
        let path = self.make_path("results.sage.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        let headers = csv::ByteRecord::from(vec![
            "peptide",
            "proteins",
            "num_proteins",
            "filename",
            "scannr",
            "rank",
            "label",
            "expmass",
            "calcmass",
            "charge",
            "peptide_len",
            "missed_cleavages",
            "isotope_error",
            "precursor_ppm",
            "fragment_ppm",
            "hyperscore",
            "delta_hyperscore",
            "rt",
            "aligned_rt",
            "predicted_rt",
            "delta_rt_model",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "matched_intensity_pct",
            "scored_candidates",
            "poisson",
            "sage_discriminant_score",
            "posterior_error",
            "spectrum_fdr",
            "peptide_fdr",
            "protein_fdr",
            "ms1_intensity",
            "ms2_intensity",
            // "ms1_apex_intensity",
        ]);

        wtr.write_byte_record(&headers)?;
        for record in features
            .into_par_iter()
            .map(|feat| self.serialize_feature(&feat, filenames))
            .collect::<Vec<_>>()
        {
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    pub fn write_tmt(&self, quant: Vec<TmtQuant>, filenames: &[String]) -> anyhow::Result<String> {
        let path = self.make_path("tmt.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = csv::ByteRecord::from(vec!["file", "scannr", "ion_injection_time"]);
        headers.extend(
            self.parameters
                .quant
                .tmt
                .as_ref()
                .map(|tmt| tmt.headers())
                .expect("TMT quant cannot be performed without setting this parameter"),
        );

        wtr.write_byte_record(&headers)?;

        let records = quant
            .into_par_iter()
            .map(|q| {
                let mut record = csv::ByteRecord::new();
                record.push_field(filenames[q.file_id].as_bytes());
                record.push_field(q.spec_id.as_bytes());
                record.push_field(ryu::Buffer::new().format(q.ion_injection_time).as_bytes());
                for peak in &q.peaks {
                    record.push_field(ryu::Buffer::new().format(*peak).as_bytes());
                }
                record
            })
            .collect::<Vec<csv::ByteRecord>>();

        for record in records {
            wtr.write_record(&record)?;
        }
        wtr.flush()?;

        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    pub fn write_lfq(
        &self,
        areas: HashMap<PeptideIx, Vec<f64>>,
        filenames: &[String],
    ) -> anyhow::Result<String> {
        let path = self.make_path("lfq.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = csv::ByteRecord::from(vec!["peptide"]);
        headers.extend(filenames);

        wtr.write_byte_record(&headers)?;

        let records = areas
            .into_par_iter()
            .map(|(peptide_ix, data)| {
                let mut record = csv::ByteRecord::new();
                record.push_field(self.database[peptide_ix].to_string().as_bytes());
                for x in data {
                    record.push_field(ryu::Buffer::new().format(x).as_bytes());
                }
                record
            })
            .collect::<Vec<csv::ByteRecord>>();

        for record in records {
            wtr.write_record(&record)?;
        }
        wtr.flush()?;

        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }
}
