use std::time::Duration;

use super::*;

#[derive(Serialize)]
struct SearchDatabase {
    id: &'static str,
    location: String,
    #[serde(rename = "FileFormat")]
    file_format: Param,
    #[serde(rename = "DatabaseName")]
    database_name: UParam,

    #[serde(rename = "cvParam")]
    params: Vec<CvParam>,
}

#[derive(Serialize)]
pub struct SourceFile {
    id: String,
    location: String,
    #[serde(rename = "FileFormat")]
    file_format: Param,
}

#[derive(Serialize)]
pub struct SpectraData {
    id: String,
    location: String,
    #[serde(rename = "FileFormat")]
    file_format: Param,
    #[serde(rename = "SpectrumIDFormat")]
    spectrum_id_format: Param,
}

#[derive(Serialize)]
#[serde(rename_all = "PascalCase")]
pub struct Inputs {
    source_file: SourceFile,
    search_database: SearchDatabase,
    spectra_data: Vec<SpectraData>,
}

#[derive(Serialize)]
pub struct PeptideEvidenceRef {
    #[serde(rename = "peptideEvidence_ref")]
    id: String,
}

/// An identification of a single (poly)peptide, resulting from querying an input spectra, along with
/// the set of confidence values for that identification. PeptideEvidence elements should be given
/// for all mappings of the corresponding Peptide sequence within protein sequences.
#[derive(Serialize)]
pub struct SpectrumIdentificationItem {
    id: String,
    #[serde(rename = "experimentalMassToCharge")]
    expmass: f32,
    #[serde(rename = "calculatedMassToCharge")]
    calcmass: f32,
    #[serde(rename = "chargeState")]
    charge: u8,
    #[serde(rename = "passThreshold")]
    pass_threshold: bool,
    peptide_ref: String,
    rank: usize,

    #[serde(rename = "PeptideEvidenceRef")]
    evidence: Vec<PeptideEvidenceRef>,
    #[serde(rename = "cvParam")]
    params: Vec<CvParam>,
    // May supply MS:1001405
}

/// All identifications made from searching one spectrum. For MS/MS data, there will be ranked
/// SpectrumIdentificationItems corresponding to possible different peptide IDs.
#[derive(Serialize)]
pub struct SpectrumIdentificationResult {
    id: String,
    #[serde(rename = "spectraData_ref")]
    spectrum_data_ref: String,
    #[serde(rename = "spectrumID")]
    spectrum_id: String,
    #[serde(rename = "SpectrumIdentificationItem")]
    sii: Vec<SpectrumIdentificationItem>,

    #[serde(rename = "cvParam")]
    params: Vec<CvParam>,
}

#[derive(Serialize)]
pub struct SpectrumIdentificationList {
    id: &'static str,
    #[serde(rename = "SpectrumIdentificationResult")]
    sir: Vec<SpectrumIdentificationResult>,

    #[serde(rename = "cvParam")]
    params: Vec<CvParam>,
}

#[derive(Serialize)]
pub struct AnalysisData {
    #[serde(rename = "SpectrumIdentificationList")]
    sil: SpectrumIdentificationList,
}

#[derive(Serialize)]
#[serde(rename_all = "PascalCase")]
pub struct DataCollection {
    inputs: Inputs,
    analysis_data: AnalysisData,
}

impl AnalysisData {
    pub fn new(features: &[sage_core::scoring::Feature], search_time: Duration) -> Self {
        let mut sir = Vec::new();
        for (idx, feat) in features.iter().enumerate() {
            let evidence = feat
                .proteins
                .split(';')
                .map(|protein| PeptideEvidenceRef {
                    id: format!("evi_{}_{}", feat.peptide, protein),
                })
                .collect();

            // MS:1001491 - percolator q-value

            let ir = SpectrumIdentificationResult {
                id: format!("SIR_{}", idx),
                spectrum_data_ref: format!("spectrum_data_{}", feat.file_id),
                spectrum_id: format!("controllerType=0 controllerNumber=1 scan={}", feat.scannr),
                sii: vec![SpectrumIdentificationItem {
                    id: format!("SII_{}_{}", idx, feat.peptide_idx.0),
                    expmass: feat.expmass,
                    calcmass: feat.calcmass,
                    charge: feat.charge,
                    pass_threshold: feat.q_value <= 0.01,
                    peptide_ref: format!("peptide_{}", feat.peptide),
                    rank: 1,
                    evidence,
                    params: vec![
                        CvParam {
                            cv_ref: "PSI-MS",
                            accession: "MS:1001331",
                            name: "X!Tandem hyperscore",
                            value: Some(feat.hyperscore.to_string()),
                            ..Default::default()
                        },
                        CvParam {
                            cv_ref: "PSI-MS",
                            accession: "MS:1002351",
                            name: "PSM-level local FDR",
                            value: Some(feat.posterior_error.to_string()),
                            ..Default::default()
                        },
                        CvParam {
                            cv_ref: "PSI-MS",
                            accession: "MS:1002351",
                            name: "PSM-level q-value",
                            value: Some(feat.q_value.to_string()),
                            ..Default::default()
                        },
                        CvParam {
                            cv_ref: "PSI-MS",
                            accession: "MS:1002500",
                            name: "peptide passes threshold",
                            value: Some((feat.peptide_q <= 0.01).to_string()),
                            ..Default::default()
                        },
                        CvParam {
                            cv_ref: "PSI-MS",
                            accession: "MS:1001868",
                            name: "distinct peptide-level q-value",
                            value: Some(feat.peptide_q.to_string()),
                            ..Default::default()
                        },
                        CvParam {
                            cv_ref: "PSI-MS",
                            accession: "MS:1002500",
                            name: "protein group passes threshold",
                            value: Some((feat.protein_q <= 0.01).to_string()),
                            ..Default::default()
                        },
                        CvParam {
                            cv_ref: "PSI-MS",
                            accession: "MS:1002373",
                            name: "protein group-level q-value",
                            value: Some(feat.protein_q.to_string()),
                            ..Default::default()
                        },
                    ],
                }],
                params: vec![CvParam {
                    cv_ref: "PSI-MS",
                    accession: "MS:1000894",
                    name: "retention time",
                    value: Some(feat.rt.to_string()),
                    unit_accession: Some("second"),
                    unit_cv_ref: Some("UO"),
                    ..Default::default()
                }],
            };

            sir.push(ir);
        }

        AnalysisData {
            sil: SpectrumIdentificationList {
                id: "SIL_0",
                sir,
                params: vec![CvParam {
                    cv_ref: "PSI-MS",
                    accession: "MS:1001036",
                    name: "search time taken",
                    value: Some(search_time.as_secs().to_string()),
                    ..Default::default()
                }],
            },
        }
    }
}

impl DataCollection {
    pub fn new(
        features: &[sage_core::scoring::Feature],
        files: &[String],
        search_time: Duration,
    ) -> Self {
        let search = SearchDatabase {
            id: "search_database",
            location: "bar".into(),
            file_format: CvParam {
                cv_ref: "PSI-MS",
                accession: "MS:1001348",
                name: "FASTA FORMAT",
                ..Default::default()
            }
            .into(),
            database_name: UserParam {
                name: "supplied fasta file",
                value: None,
            }
            .into(),
            params: vec![CvParam {
                cv_ref: "PSI-MS",
                accession: "MS:1001197",
                name: "DB composition target+decoy",
                ..Default::default()
            }
            .into()],
        };

        let spectra_data = files
            .iter()
            .enumerate()
            .map(|(idx, file)| SpectraData {
                id: format!("spectrum_data_{}", idx),
                location: file.clone(),
                file_format: CvParam {
                    cv_ref: "PSI-MS",
                    accession: "MS:1000584",
                    name: "mzML format",
                    ..Default::default()
                }
                .into(),
                spectrum_id_format: CvParam {
                    cv_ref: "PSI-MS",
                    accession: "MS:1000776",
                    name: "scan number only nativeID format",
                    ..Default::default()
                }
                .into(),
            })
            .collect();

        DataCollection {
            inputs: Inputs {
                source_file: SourceFile {
                    id: "config.json".into(),
                    location: "config.json".into(),
                    file_format: CvParam {
                        cv_ref: "PSI-MS",
                        accession: "MS:1001401",
                        name: "X!Tandem xml file",
                        ..Default::default()
                    }
                    .into(),
                },
                search_database: search,
                spectra_data,
            },
            analysis_data: AnalysisData::new(features, search_time),
        }
    }
}
