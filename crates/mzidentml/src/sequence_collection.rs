use std::collections::{HashMap, HashSet};

use sage_core::database::IndexedDatabase;

use super::*;

#[derive(Serialize)]
struct PeptideSequence(String);

#[derive(Serialize)]
pub struct Modification {
    #[serde(rename = "monoisotopicMassDelta")]
    mass_delta: f32,
    location: usize,
    #[serde(rename = "cvParam")]
    param: CvParam,
}

#[derive(Serialize)]
struct Peptide {
    id: String,
    #[serde(rename = "PeptideSequence")]
    peptide_sequence: PeptideSequence,
    #[serde(rename = "Modification")]
    modifications: Vec<Modification>,
}

#[derive(Serialize)]
struct PeptideEvidence {
    id: String,
    #[serde(rename = "isDecoy")]
    is_decoy: bool,
    peptide_ref: String,
    #[serde(rename = "dBSequence_ref")]
    db_sequence_ref: String,
    start: u32,
    end: u32,
}

#[derive(Serialize)]
struct DBSequence {
    accession: String,
    id: String,
    #[serde(rename = "searchDatabase_ref")]
    search_db_ref: &'static str,
}

#[derive(Serialize)]
pub struct SequenceCollection {
    #[serde(rename = "DBSequence")]
    db_sequence: Vec<DBSequence>,
    #[serde(rename = "Peptide")]
    peptide: Vec<Peptide>,
    #[serde(rename = "PeptideEvidence")]
    evidence: Vec<PeptideEvidence>,
}

impl SequenceCollection {
    pub fn new(features: &[sage_core::scoring::Feature], db: &IndexedDatabase) -> Self {
        let mut peptides = Vec::new();
        let mut evidence_set: HashSet<&str> = HashSet::new();
        let mut evidence = Vec::new();

        let mut proteins: HashMap<&str, DBSequence> = HashMap::new();
        for feat in features {
            if evidence_set.insert(&feat.peptide) {
                let mut modifications = Vec::new();
                let mut base = String::new();
                for (idx, seq) in db[feat.peptide_idx].sequence.iter().enumerate() {
                    match seq {
                        sage_core::mass::Residue::Just(c) => base.push(*c),
                        sage_core::mass::Residue::Mod(c, m) => {
                            base.push(*c);
                            modifications.push(Modification {
                                mass_delta: *m,
                                location: idx + 1,
                                param: CvParam {
                                    cv_ref: "PSI-MS",
                                    accession: "MS:1001460",
                                    name: "unknown modification",
                                    ..Default::default()
                                },
                            })
                        }
                    }
                }
                peptides.push(Peptide {
                    id: format!("peptide_{}", feat.peptide),
                    peptide_sequence: PeptideSequence(base),
                    modifications,
                });

                for protein in feat.proteins.split(';') {
                    let protein = protein.trim();
                    if !proteins.contains_key(protein) {
                        let cleaned = if protein.starts_with("rev_") {
                            let mut parts = protein
                                .trim_start_matches("rev_")
                                .split("|")
                                .map(String::from)
                                .collect::<Vec<_>>();
                            parts[1] = format!("rev_{}", parts[1]);
                            parts.join("|")
                        } else {
                            protein.into()
                        };

                        proteins.insert(
                            protein,
                            DBSequence {
                                accession: cleaned,
                                id: protein.into(),
                                search_db_ref: "search_database",
                            },
                        );
                    }

                    evidence.push(PeptideEvidence {
                        id: format!("evi_{}_{}", feat.peptide, protein),
                        is_decoy: feat.label == -1,
                        peptide_ref: format!("peptide_{}", feat.peptide),
                        db_sequence_ref: protein.into(),
                        start: db[feat.peptide_idx].start + 1,
                        end: db[feat.peptide_idx].end,
                    });
                }
            }
        }

        SequenceCollection {
            db_sequence: proteins.into_values().collect(),
            peptide: peptides,
            evidence: evidence,
        }
    }
}
