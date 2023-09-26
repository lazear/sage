use std::collections::HashMap;

use rayon::prelude::*;
use sage_core::{
    lfq::{Peak, PrecursorId},
    scoring::Feature,
    tmt::TmtQuant,
};

use crate::Runner;

impl Runner {
    pub fn serialize_feature(&self, feature: &Feature, filenames: &[String]) -> csv::ByteRecord {
        let mut record = csv::ByteRecord::new();
        let peptide = &self.database[feature.peptide_idx];
        record.push_field(peptide.to_string().as_bytes());
        record.push_field(
            peptide
                .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(peptide.proteins.len())
                .as_bytes(),
        );
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
        record.push_field(ryu::Buffer::new().format(feature.delta_next).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.delta_best).as_bytes());
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
        features: &[Feature],
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
            "delta_next",
            "delta_best",
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
            "spectrum_q",
            "peptide_q",
            "protein_q",
            "ms1_intensity",
            "ms2_intensity",
        ]);

        wtr.write_byte_record(&headers)?;
        for record in features
            .into_par_iter()
            .map(|feat| self.serialize_feature(feat, filenames))
            .collect::<Vec<_>>()
        {
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    fn serialize_pin(
        &self,
        re: &regex::Regex,
        idx: usize,
        feature: &Feature,
        filenames: &[String],
    ) -> csv::ByteRecord {
        let scannr = re
            .captures_iter(&feature.spec_id)
            .last()
            .and_then(|cap| cap.get(1).map(|cap| cap.as_str()))
            .unwrap_or(&feature.spec_id);

        let mut record = csv::ByteRecord::new();
        let peptide = &self.database[feature.peptide_idx];
        record.push_field(itoa::Buffer::new().format(idx).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.label).as_bytes());
        record.push_field(scannr.as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.expmass).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.calcmass).as_bytes());
        record.push_field(filenames[feature.file_id].as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.rt).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.rank).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 2) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 3) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 4) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 5) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format((feature.charge == 6) as i32)
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(
                    (feature.charge < 2 || feature.charge > 6)
                        .then_some(feature.charge)
                        .unwrap_or(0),
                )
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.peptide_len).as_bytes());
        record.push_field(
            itoa::Buffer::new()
                .format(feature.missed_cleavages)
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.isotope_error).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_mass.ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.average_ppm).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.hyperscore.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_next.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_best.ln_1p())
                .as_bytes(),
        );
        record.push_field(ryu::Buffer::new().format(feature.aligned_rt).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.predicted_rt).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.delta_rt_model.clamp(0.001, 1.0).sqrt())
                .as_bytes(),
        );
        record.push_field(itoa::Buffer::new().format(feature.matched_peaks).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_b).as_bytes());
        record.push_field(itoa::Buffer::new().format(feature.longest_y).as_bytes());
        record.push_field(ryu::Buffer::new().format(feature.longest_y_pct).as_bytes());
        record.push_field(
            ryu::Buffer::new()
                .format(feature.matched_intensity_pct.ln_1p())
                .as_bytes(),
        );
        record.push_field(
            itoa::Buffer::new()
                .format(feature.scored_candidates)
                .as_bytes(),
        );
        record.push_field(
            ryu::Buffer::new()
                .format((-feature.poisson).ln_1p())
                .as_bytes(),
        );
        record.push_field(peptide.to_string().as_bytes());
        record.push_field(
            peptide
                .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                .as_bytes(),
        );
        record
    }

    pub fn write_pin(&self, features: &[Feature], filenames: &[String]) -> anyhow::Result<String> {
        let path = self.make_path("results.sage.pin");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);

        let headers = csv::ByteRecord::from(vec![
            "SpecId",
            "Label",
            "ScanNr",
            "ExpMass",
            "CalcMass",
            "FileName",
            "retentiontime",
            "rank",
            "z=2",
            "z=3",
            "z=4",
            "z=5",
            "z=6",
            "z=other",
            "peptide_len",
            "missed_cleavages",
            "isotope_error",
            "ln(precursor_ppm)",
            "fragment_ppm",
            "ln(hyperscore)",
            "ln(delta_next)",
            "ln(delta_best)",
            "aligned_rt",
            "predicted_rt",
            "sqrt(delta_rt_model)",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "ln(matched_intensity_pct)",
            "scored_candidates",
            "ln(-poisson)",
            "Peptide",
            "Proteins",
        ]);

        let re = regex::Regex::new(r"scan=(\d+)").expect("This is valid regex");

        wtr.write_byte_record(&headers)?;
        for record in features
            .into_par_iter()
            .enumerate()
            .map(|(idx, feat)| self.serialize_pin(&re, idx, feat, filenames))
            .collect::<Vec<_>>()
        {
            wtr.write_byte_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        path.write_bytes_sync(bytes)?;
        Ok(path.to_string())
    }

    pub fn write_tmt(&self, quant: &[TmtQuant], filenames: &[String]) -> anyhow::Result<String> {
        let path = self.make_path("tmt.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = csv::ByteRecord::from(vec!["filename", "scannr", "ion_injection_time"]);
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
        areas: HashMap<(PrecursorId, bool), (Peak, Vec<f64>), fnv::FnvBuildHasher>,
        filenames: &[String],
    ) -> anyhow::Result<String> {
        let path = self.make_path("lfq.tsv");

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = csv::ByteRecord::from(vec![
            "peptide",
            "charge",
            "proteins",
            "q_value",
            "score",
            "spectral_angle",
        ]);
        headers.extend(filenames);

        wtr.write_byte_record(&headers)?;

        let records = areas
            .into_par_iter()
            .filter_map(|((id, decoy), (peak, data))| {
                if decoy {
                    return None;
                };
                let mut record = csv::ByteRecord::new();
                let (peptide_ix, charge) = match id {
                    PrecursorId::Combined(x) => (x, None),
                    PrecursorId::Charged((x, charge)) => (x, Some(charge as i32)),
                };
                record.push_field(self.database[peptide_ix].to_string().as_bytes());
                record.push_field(itoa::Buffer::new().format(charge.unwrap_or(-1)).as_bytes());
                record.push_field(
                    self.database[peptide_ix]
                        .proteins(&self.database.decoy_tag, self.database.generate_decoys)
                        .as_bytes(),
                );
                record.push_field(ryu::Buffer::new().format(peak.q_value).as_bytes());
                record.push_field(ryu::Buffer::new().format(peak.score).as_bytes());
                record.push_field(ryu::Buffer::new().format(peak.spectral_angle).as_bytes());
                for x in data {
                    record.push_field(ryu::Buffer::new().format(x).as_bytes());
                }
                Some(record)
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
