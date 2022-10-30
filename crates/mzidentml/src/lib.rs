use std::{io::Write, path, time::Duration};

use sage_core::{database::IndexedDatabase, mass::Tolerance};
use serde::Serialize;
mod analysis_protocol;
mod data_collection;
mod sequence_collection;
use analysis_protocol::AnalysisProtocolCollection;
use data_collection::DataCollection;
use sequence_collection::SequenceCollection;

const XSI: &'static str = "http://www.w3.org/2001/XMLSchema-instance";
const SCHEMA_LOCATION: &'static str =
    "http://psidev.info/psi/pi/mzIdentML/1.2 http://www.psidev.info/files/mzIdentML1.2.0.xsd";
const XMLNS: &'static str = "http://psidev.info/psi/pi/mzIdentML/1.2";

#[derive(Serialize, Default)]
#[serde(rename_all = "camelCase")]
pub struct CvParam {
    cv_ref: &'static str,
    accession: &'static str,
    name: &'static str,
    value: Option<String>,
    unit_cv_ref: Option<&'static str>,
    unit_accession: Option<&'static str>,
    unit_name: Option<&'static str>,
}

#[derive(Serialize)]
pub struct Param {
    #[serde(rename = "cvParam")]
    p: CvParam,
}

impl From<CvParam> for Param {
    fn from(p: CvParam) -> Self {
        Self { p }
    }
}

#[derive(Serialize)]
pub struct UserParam {
    name: &'static str,
    value: Option<String>,
}

#[derive(Serialize)]
pub struct UParam {
    #[serde(rename = "userParam")]
    p: UserParam,
}

impl From<UserParam> for UParam {
    fn from(p: UserParam) -> Self {
        Self { p }
    }
}

#[derive(Serialize)]
pub struct SearchDatabaseRef {
    #[serde(rename = "searchDatabase_ref")]
    id: &'static str,
}

#[derive(Serialize)]
pub struct InputSpectra {
    #[serde(rename = "spectraData_ref")]
    id: String,
}

#[derive(Serialize)]
pub struct SpectrumIdentification {
    id: &'static str,
    #[serde(rename = "spectrumIdentificationList_ref")]
    sil: &'static str,
    #[serde(rename = "spectrumIdentificationProtocol_ref")]
    sip: &'static str,

    #[serde(rename = "InputSpectra")]
    input: Vec<InputSpectra>,

    #[serde(rename = "SearchDatabaseRef")]
    sdr: SearchDatabaseRef,
}

#[derive(Serialize)]
pub struct AnalysisCollection {
    #[serde(rename = "SpectrumIdentification")]
    s: SpectrumIdentification,
}

#[derive(Serialize)]
pub struct Cv {
    id: &'static str,
    uri: &'static str,
    #[serde(rename = "fullName")]
    name: &'static str,
}

#[derive(Serialize)]
pub struct CvList {
    #[serde(rename = "cv")]
    cv_list: Vec<Cv>,
}

#[derive(Serialize)]
pub struct AnalysisSoftware {
    id: &'static str,
    name: &'static str,
    uri: &'static str,
    version: &'static str,

    #[serde(rename = "SoftwareName")]
    sw_name: Param,
}

#[derive(Serialize)]
pub struct AnalysisSoftwareList {
    #[serde(rename = "AnalysisSoftware")]
    software: AnalysisSoftware,
}

#[derive(Serialize)]
pub struct MzIdentML {
    version: &'static str,
    id: &'static str,

    #[serde(rename = "xmlns:xsi")]
    xsi: &'static str,
    #[serde(rename = "xsi:schemaLocation")]
    schema_location: &'static str,
    xmlns: &'static str,

    #[serde(rename = "cvList")]
    cv_list: CvList,

    #[serde(rename = "AnalysisSoftwareList")]
    software: AnalysisSoftwareList,

    #[serde(rename = "SequenceCollection")]
    sequence_collection: SequenceCollection,
    #[serde(rename = "AnalysisCollection")]
    analysis_collection: AnalysisCollection,
    #[serde(rename = "AnalysisProtocolCollection")]
    protocol_collection: AnalysisProtocolCollection,
    #[serde(rename = "DataCollection")]
    data_collection: DataCollection,
}

pub fn doit(
    paths: &[String],
    features: &[sage_core::scoring::Feature],
    db: &IndexedDatabase,
    fragment: Tolerance,
    precursor: Tolerance,
    search_time: Duration,
) {
    // let mut buffer = Vec::new();
    let mut buffer = std::fs::File::create("mzid.xml").unwrap();
    buffer
        .write(r#"<?xml version="1.0" encoding="UTF-8" ?>"#.as_bytes())
        .unwrap();
    buffer.write(b"\n").unwrap();
    let mut wtr = quick_xml::se::Serializer::new(&mut buffer);

    let mz = MzIdentML {
        version: "1.2.0",
        id: "Sage:MzIdentML-v1.2.0",
        xmlns: XMLNS,
        xsi: XSI,
        schema_location: SCHEMA_LOCATION,
        cv_list: CvList {
            cv_list: vec![
                // <cv id="PSI-MS" uri="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo" fullName="PSI-MS"/>
                Cv {
                    id: "PSI-MS",
                    uri: "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo",
                    name: "PSI-MS",
                },
                // <cv id="UNIMOD" uri="http://www.unimod.org/obo/unimod.obo" fullName="UNIMOD"/>
                Cv {
                    id: "UNIMOD",
                    uri: "http://www.unimod.org/obo/unimod.obo",
                    name: "UNIMOD",
                },
                Cv {
                    id: "UO",
                    uri: "https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo",
                    name: "UNIT-ONTOLOGY",
                }
            ],
        },
        software: AnalysisSoftwareList {
            software: AnalysisSoftware {
                id: "sage",
                name: "Sage",
                uri: "https://github.com/lazear/sage",
                version: "v0.5.0",
                sw_name: CvParam {
                    cv_ref: "PSI-MS",
                    accession: "MS:1001476",
                    name: "X!Tandem",
                    ..Default::default()
                }
                .into(),
            },
        },
        sequence_collection: SequenceCollection::new(features, db),
        analysis_collection: AnalysisCollection {
            s: SpectrumIdentification {
                id: "SI_0",
                sil: "SIL_0",
                sip: "SIP_0",
                input: paths
                    .iter()
                    .enumerate()
                    .map(|(idx, _)| InputSpectra {
                        id: format!("spectrum_data_{}", idx),
                    })
                    .collect(),
                sdr: SearchDatabaseRef {
                    id: "search_database",
                },
            },
        },
        protocol_collection: AnalysisProtocolCollection::new(fragment, precursor),
        data_collection: DataCollection::new(features, paths, search_time),
    };

    mz.serialize(&mut wtr).unwrap();
}
