use sage_core::mass::Tolerance;

use super::*;

#[derive(Serialize)]
pub struct Threshold {
    #[serde(rename = "cvParam")]
    cv: Vec<CvParam>,
}

#[derive(Serialize)]
pub struct Tol {
    #[serde(rename = "cvParam")]
    cv: Vec<CvParam>,
}

#[derive(Serialize)]
pub struct SpectrumIdProtcol {
    #[serde(rename = "analysisSoftware_ref")]
    sofware_ref: &'static str,
    id: &'static str,
    #[serde(rename = "SearchType")]
    search_type: Param,

    #[serde(rename = "FragmentTolerance")]
    fragment_tolerance: Tol,

    #[serde(rename = "ParentTolerance")]
    parent_tolerance: Tol,

    #[serde(rename = "Threshold")]
    threshold: Threshold,
}

fn tol_to_cv(tolerance: sage_core::mass::Tolerance) -> Tol {
    let cv = match tolerance {
        sage_core::mass::Tolerance::Ppm(lo, hi) => vec![
            CvParam {
                cv_ref: "PSI-MS",
                accession: "MS:1001412",
                name: "search tolerance plus value",
                value: Some(lo.to_string()),
                unit_accession: Some("UO:0000169".into()),
                unit_name: Some("parts per million".into()),
                unit_cv_ref: Some("UO".into()),
            },
            CvParam {
                cv_ref: "PSI-MS",
                accession: "MS:1001413",
                name: "search tolerance minus value",
                value: Some(hi.to_string()),
                unit_accession: Some("UO:0000169".into()),
                unit_name: Some("parts per million".into()),
                unit_cv_ref: Some("UO".into()),
            },
        ],
        sage_core::mass::Tolerance::Da(lo, hi) => vec![
            CvParam {
                cv_ref: "PSI-MS",
                accession: "MS:1001412",
                name: "search tolerance plus value",
                value: Some(lo.to_string()),
                unit_accession: Some("UO:0000221".into()),
                unit_name: Some("dalton".into()),
                unit_cv_ref: Some("UO".into()),
            },
            CvParam {
                cv_ref: "PSI-MS",
                accession: "MS:1001413",
                name: "search tolerance minus value",
                value: Some(hi.to_string()),
                unit_accession: Some("UO:0000221".into()),
                unit_name: Some("dalton".into()),
                unit_cv_ref: Some("UO".into()),
            },
        ],
    };
    Tol { cv }
}

#[derive(Serialize)]
pub struct AnalysisProtocolCollection {
    #[serde(rename = "SpectrumIdentificationProtocol")]
    sid: SpectrumIdProtcol,
}

impl AnalysisProtocolCollection {
    pub fn new(fragment: Tolerance, precursor: Tolerance) -> Self {
        let sid = SpectrumIdProtcol {
            sofware_ref: "sage",
            id: "SIP_0",
            search_type: CvParam {
                cv_ref: "PSI-MS",
                accession: "MS:1001083",
                name: "ms-ms search",
                ..Default::default()
            }
            .into(),
            fragment_tolerance: tol_to_cv(fragment),
            parent_tolerance: tol_to_cv(precursor),
            threshold: Threshold {
                cv: vec![
                    CvParam {
                        cv_ref: "PSI-MS",
                        accession: "MS:1002350",
                        name: "PSM-level global FDR",
                        value: Some("1.0".into()),
                        ..Default::default()
                    }
                    .into(),
                    CvParam {
                        cv_ref: "PSI-MS",
                        accession: "MS:1001194",
                        name: "quality estimation with decoy database",
                        ..Default::default()
                    }
                    .into(),
                    CvParam {
                        cv_ref: "PSI-MS",
                        accession: "MS:1002260",
                        name: "PSM:FDR threshold",
                        value: Some("1.0".into()),
                        ..Default::default()
                    }
                    .into(),
                    CvParam {
                        cv_ref: "PSI-MS",
                        accession: "MS:1001448",
                        name: "pep:FDR threshold",
                        value: Some("1.0".into()),
                        ..Default::default()
                    }
                    .into(),
                    CvParam {
                        cv_ref: "PSI-MS",
                        accession: "MS:1002910",
                        name: "proteoform-level global FDR threshold",
                        value: Some("1.0".into()),
                        ..Default::default()
                    }
                    .into(),
                ],
            },
        };
        Self { sid }
    }
}
