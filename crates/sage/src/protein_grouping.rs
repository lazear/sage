//! # Protein Grouping and Inference
//!
//! This module provides functionality for grouping proteins based on peptide evidence
//! and performing protein inference for proteomics data analysis. It supports grouping
//! proteins and inferring (almost) minimal covering sets and annotating features with
//! their corresponding protein groups.
//!
//! ## Main Features
//! - Groups proteins when peptide evidence for proteins is identical.
//! - Infers (almost) minimal protein-group covers using bipartite graph algorithms.
//! - Supports different protein inference strategies (e.g., "All", "Slim").
//! - Annotates features with protein group information in parallel for efficiency.
//! - Handles decoy proteins and supports deterministic group naming.
//!
//! ## Key Types
//! - `ProteinGroup`: Represents a group of proteins with identical peptide evidence.
//! - `ProteinGrouping`: Builds protein groups and manages peptide-to-protein mapping.
//! - `ProteinGroupMap`: Maps proteins to their groups for annotation.
//! - `ProteinInference`: Enum specifying the strategy for protein inference,
//! such as reporting all possible protein groups ("All") or inferring an (almost)
//! minimal covering set ("Slim").
//!
//! ## Usage
//! Use [`generate_proteingroups`] to annotate features with protein group information
//! based on the provided indexed database and peptide identifications.
//!
//! ## Parallelization
//! The module leverages Rayon for parallel processing of features to improve performance
//! on large datasets.
//!
//! ## Dependencies
//! Relies on FnvHashMap/FnvHashSet for efficient hashing, Rayon for parallelism, and
//! Itertools for convenient iterator utilities.
//!
//! ## References
// 1. Zhang, B., Chambers, M. C., & Tabb, D. L. (2007). Proteomic parsimony through
// bipartite graph analysis improves accuracy and transparency. Journal of proteome research,
// 6(9), 3549-3557. https://doi.org/10.1021/pr070230d
//!

use crate::database::{IndexedDatabase, PeptideIx};
use crate::peptide::Peptide;
use crate::scoring::Feature;
use fnv::{FnvHashMap, FnvHashSet};
use itertools::Itertools;
use log::info;
use rayon::prelude::*;
use std::hash::Hash;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;
use std::time::Instant;

const PROTEIN_INFERENCE: ProteinInference = ProteinInference::Slim;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct ProteinIx(u32);

impl ProteinIx {
    fn to_string(
        &self,
        protein_map: &Vec<(Arc<str>, bool)>,
        decoy_tag: &str,
        generate_decoys: bool,
    ) -> String {
        let (prot, decoy) = &protein_map[self.0 as usize];
        if *decoy & generate_decoys {
            format!("{}{}", decoy_tag, prot)
        } else {
            prot.to_string()
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Default, PartialOrd, Ord)]
struct ProteinGroup(Vec<ProteinIx>);

impl ProteinGroup {
    fn to_string(
        &self,
        protein_map: &Vec<(Arc<str>, bool)>,
        decoy_tag: &str,
        generate_decoys: bool,
    ) -> String {
        self.0
            .iter()
            .map(|prot_ix| prot_ix.to_string(protein_map, decoy_tag, generate_decoys))
            .sorted()
            .join("/")
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ProteinInference {
    All,
    #[default]
    Slim,
}

impl ProteinInference {
    fn infer(
        &self,
        connections: Vec<(u32, u32)>,
        peptide_count: usize,
        protein_count: usize,
    ) -> Vec<bool> {
        match self {
            ProteinInference::All => vec![true; protein_count],
            ProteinInference::Slim => {
                let mut graph = BipartiteGraph::new(connections, protein_count, peptide_count);
                while !graph.is_empty() {
                    graph.trim_connections();
                    if !graph.is_empty() {
                        graph.add_biggest_left_to_cover();
                    }
                }
                graph.left_cover
            }
        }
    }
}

struct BipartiteGraph {
    connections: Vec<(u32, u32)>,
    original_left: Vec<u32>,
    remaining_left: Vec<u32>,
    remaining_right: Vec<u32>,
    left_cover: Vec<bool>,
    right_cover: Vec<bool>,
}

impl BipartiteGraph {
    fn new(connections: Vec<(u32, u32)>, left_size: usize, right_size: usize) -> Self {
        let mut remaining_left = vec![0; left_size];
        let mut remaining_right = vec![0; right_size];
        connections.iter().for_each(|(l, r)| {
            remaining_left[*l as usize] += 1;
            remaining_right[*r as usize] += 1;
        });
        let left_cover = vec![false; left_size];
        let right_cover = vec![false; right_size];
        Self {
            connections,
            original_left: remaining_left.clone(),
            remaining_left,
            remaining_right,
            left_cover,
            right_cover,
        }
    }

    fn is_empty(&self) -> bool {
        self.connections.is_empty()
    }

    fn remove_edge(&mut self, left_index: usize, right_index: usize) {
        self.remaining_left[left_index] -= 1;
        self.remaining_right[right_index] -= 1;
    }

    fn trim_connections(&mut self) {
        let mut connections = self.connections.to_owned();
        let mut connection_count = 0;
        while connection_count != connections.len() {
            connection_count = connections.len();
            connections.iter().for_each(|(left_index, right_index)| {
                let left_index = *left_index as usize;
                let right_index = *right_index as usize;
                if self.remaining_right[right_index] == 1 {
                    self.left_cover[left_index] = true;
                }
            });
            connections.retain(|(left_index, right_index)| {
                let left_index = *left_index as usize;
                let right_index = *right_index as usize;
                if self.left_cover[left_index] {
                    self.right_cover[right_index] = true;
                    self.remove_edge(left_index, right_index);
                    false
                } else {
                    true
                }
            });
            connections.retain(|(left_index, right_index)| {
                let left_index = *left_index as usize;
                let right_index = *right_index as usize;
                if self.right_cover[right_index] {
                    self.remove_edge(left_index, right_index);
                    false
                } else {
                    true
                }
            });
        }
        self.connections = connections;
    }

    fn add_biggest_left_to_cover(&mut self) {
        match self
            .remaining_left
            .iter()
            .zip(&self.original_left)
            .enumerate()
            .max_by_key(|(_, i)| *i)
        {
            Some((index, _)) => self.left_cover[index] = true,
            _ => {}
        }
    }
}

#[derive(Debug, Default)]
struct ProteinGrouping {
    protein_groups: Vec<ProteinGroup>,
    connections: Vec<(u32, u32)>,
    prot_name_map: FnvHashMap<(Arc<str>, bool), ProteinIx>,
    protein_cover: Vec<bool>,
    peptide_count: usize,
}

impl ProteinGrouping {
    fn new(db: &IndexedDatabase, peps: FnvHashSet<PeptideIx>) -> Self {
        let mut mapping = Self::default();
        mapping.peptide_count = peps.len();
        let time = Instant::now();
        let meta_peptides = mapping.set_proteins_and_get_metapeptides(peps, &db);
        info!(
            "-  found {} meta peptides in {:?}ms",
            meta_peptides.len(),
            time.elapsed().as_millis()
        );
        let time = Instant::now();
        mapping.set_protein_groups_and_connections(meta_peptides);
        info!(
            "-  found {} protein groups in {:?}ms",
            mapping.protein_groups.len(),
            time.elapsed().as_millis()
        );
        mapping
    }

    fn set_proteins_and_get_metapeptides(
        &mut self,
        peps: FnvHashSet<PeptideIx>,
        db: &&IndexedDatabase,
    ) -> FnvHashSet<Vec<ProteinIx>> {
        peps.into_iter()
            .sorted() //Needed to ensure the same order for prot_name_map
            .map(|pep_id| {
                let db_peptide = &db[pep_id];
                let prot_ids = db_peptide
                    .proteins
                    .iter()
                    .map(|prot| {
                        let prot = (prot.clone(), db_peptide.decoy);
                        self.get_or_insert_prot(prot)
                    })
                    .sorted()
                    .collect::<Vec<_>>();
                prot_ids
            })
            .collect()
    }

    fn get_or_insert_prot(&mut self, prot: (Arc<str>, bool)) -> ProteinIx {
        match self.prot_name_map.get(&prot) {
            Some(&prot_id) => prot_id,
            None => {
                let prot_id = ProteinIx(self.prot_name_map.len() as u32);
                self.prot_name_map.insert(prot, prot_id);
                prot_id
            }
        }
    }

    fn set_protein_groups_and_connections<'a>(
        &mut self,
        meta_peptides: FnvHashSet<Vec<ProteinIx>>,
    ) {
        let mut protein_map = FnvHashMap::default();
        meta_peptides
            .into_iter()
            .sorted() //Needed to ensure the same order for protein_groups
            .enumerate()
            .for_each(|(i, prot_ids)| {
                prot_ids.into_iter().for_each(|prot_id| {
                    let entry = protein_map.entry(prot_id).or_insert_with(|| vec![]);
                    entry.push(i);
                });
            });
        let mut mapping = FnvHashMap::default();
        protein_map
            .into_iter()
            .for_each(|(prot_id, meta_peptides)| {
                let entry = mapping
                    .entry(meta_peptides)
                    .or_insert_with(|| ProteinGroup::default());
                entry.0.push(prot_id);
            });
        let mut protein_groups = vec![];
        let mut connections = vec![];
        mapping.into_iter().sorted().enumerate().for_each(
            |(meta_protein_index, (meta_peptides, meta_protein))| {
                protein_groups.push(meta_protein);
                meta_peptides.into_iter().for_each(|meta_peptide_index| {
                    connections.push((meta_protein_index as u32, meta_peptide_index as u32));
                });
            },
        );
        self.protein_groups = protein_groups;
        self.connections = connections;
    }
}

struct ProteinGroupMap {
    protein_groups: Vec<ProteinGroup>,
    proteins: Vec<(Arc<str>, bool)>,
    protein_group_map: FnvHashMap<(Arc<str>, bool), Vec<u32>>,
}

impl ProteinGroupMap {
    fn new(mapping: ProteinGrouping) -> Self {
        let proteins = mapping
            .prot_name_map
            .into_iter()
            .sorted_by_key(|(_, v)| v.0)
            .map(|(k, _)| k.clone())
            .collect::<Vec<_>>();
        let mut protein_group_map = FnvHashMap::default();
        mapping
            .protein_cover
            .into_iter()
            .enumerate()
            .filter(|(_, i)| *i)
            .for_each(|(i, _)| {
                let protein_group = &mapping.protein_groups[i];
                protein_group.0.iter().for_each(|&prot_id| {
                    let (name, decoy) = &proteins[prot_id.0 as usize];
                    let entry = protein_group_map
                        .entry((name.clone(), *decoy))
                        .or_insert_with(|| vec![]);
                    entry.push(i as u32);
                });
            });
        let protein_group_map = protein_group_map;
        Self {
            protein_groups: mapping.protein_groups,
            proteins,
            protein_group_map,
        }
    }

    fn get_protein_group_string(&self, peptide: &Peptide, db: &IndexedDatabase) -> Vec<String> {
        let proteingroups = peptide
            .proteins
            .iter()
            .map(|protein| (protein.clone(), peptide.decoy))
            .filter_map(|prot| self.protein_group_map.get(&prot))
            .flat_map(|protein_group_indices| {
                protein_group_indices
                    .iter()
                    .map(|&i| &self.protein_groups[i as usize])
            })
            .collect::<FnvHashSet<_>>();
        if proteingroups.is_empty() {
            vec![]
        } else {
            proteingroups
                .iter()
                .map(|p| p.to_string(&self.proteins, &db.decoy_tag, db.generate_decoys))
                .collect()
        }
    }
}

pub fn generate_proteingroups(
    db: &IndexedDatabase,
    features: &mut [Feature],
    protein_grouping: bool,
    confident_peptide_threshold: Option<f32>,
) {
    let time = Instant::now();
    if protein_grouping {
        if confident_peptide_threshold.is_some() {
            update_features_with_proteingroups(features, db, confident_peptide_threshold);
        }
        update_features_with_proteingroups(features, db, None);
    }
    features
        .par_iter_mut()
        .filter(|f| f.proteingroups.is_none())
        .for_each(|feat| {
            let pep = &db[feat.peptide_idx];
            let proteingroups =
                annotate_proteins(&pep.proteins, &db.decoy_tag, db.generate_decoys, pep.decoy)
                    .collect::<Vec<_>>();
            feat.proteingroups = Some(proteingroups.iter().sorted().join(";"));
            feat.num_proteingroups = proteingroups.len() as i32;
        });
    info!(
        "Grouped and inferred proteins in {:?}ms",
        time.elapsed().as_millis()
    );
}

fn update_features_with_proteingroups(
    features: &mut [Feature],
    db: &IndexedDatabase,
    confident_peptide_threshold: Option<f32>,
) {
    let time = Instant::now();
    info!("Protein grouping with {} features", features.len());
    let peps = if let Some(threshold) = confident_peptide_threshold {
        let t = threshold.max(0.0).min(1.0);
        features
            .par_iter()
            .filter(|f| (f.label != -1) && (f.peptide_q < t))
            .map(|feature| feature.peptide_idx)
            .collect::<FnvHashSet<_>>()
    } else {
        features
            .par_iter()
            .filter(|f| (f.label != -1))
            .map(|feature| feature.peptide_idx)
            .collect::<FnvHashSet<_>>()
    };
    info!(
        "-  found {} unique peptides in {:?}ms",
        peps.len(),
        time.elapsed().as_millis()
    );
    let mut mapping = ProteinGrouping::new(db, peps);
    let time = Instant::now();
    mapping.protein_cover = PROTEIN_INFERENCE.infer(
        mapping.connections.drain(..).collect(),
        mapping.peptide_count,
        mapping.protein_groups.len(),
    );
    info!(
        "-  found cover of {} metaproteins in {:?}ms",
        mapping.protein_cover.len(),
        time.elapsed().as_millis()
    );
    let time = Instant::now();
    let protein_map = ProteinGroupMap::new(mapping);
    let counter = AtomicUsize::new(0);
    features
        .par_iter_mut()
        .filter(|feat| feat.proteingroups.is_none())
        .for_each(|feat| {
            let pep = &db[feat.peptide_idx];
            let proteingroups = protein_map.get_protein_group_string(pep, db);
            if !proteingroups.is_empty() {
                feat.proteingroups = Some(proteingroups.iter().sorted().join(";"));
                feat.num_proteingroups = proteingroups.len() as i32;
                counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
        });
    info!(
        "-  annotated {} features in {:?}ms",
        counter.load(std::sync::atomic::Ordering::Relaxed),
        time.elapsed().as_millis()
    );
}

pub fn annotate_proteins<'a>(
    proteins: &'a [Arc<str>],
    decoy_tag: &'a str,
    generate_decoys: bool,
    decoy: bool,
) -> impl Iterator<Item = String> + 'a {
    proteins
        .iter()
        .map(move |s| {
            if decoy & generate_decoys {
                format!("{}{}", decoy_tag, s)
            } else {
                s.to_string()
            }
        })
        .sorted()
}

#[cfg(test)]
mod test {
    use crate::peptide;

    use super::*;
    use std::sync::Arc;

    fn get_data() -> (Vec<Vec<impl AsRef<str>>>, Vec<bool>, Vec<f32>) {
        let proteins = vec![
            vec!["protein_7"],
            vec!["protein_4", "protein_6", "protein_9"],
            vec!["protein_1"],
            vec!["protein_1", "protein_5"],
            vec!["protein_7"],
            vec!["protein_3", "protein_6"],
            vec!["protein_1"],
            vec!["protein_1", "protein_2", "protein_5", "protein_8"],
            vec!["protein_1"],
            vec!["protein_4", "protein_9"],
        ];
        let decoys = vec![false; proteins.len()];
        let q_vals = vec![0.0; proteins.len()];
        (proteins, decoys, q_vals)
    }

    fn build_db_and_features(
        proteins: &[Vec<impl AsRef<str>>],
        decoys: &[bool],
        q_vals: &[f32],
    ) -> (IndexedDatabase, Vec<Feature>) {
        let features: Vec<Feature> = (0..proteins.len())
            .map(|ix| Feature {
                peptide_idx: PeptideIx(ix as u32),
                label: if decoys[ix] { -1 } else { 1 },
                peptide_q: q_vals[ix],
                ..Default::default()
            })
            .collect();
        let db = IndexedDatabase {
            peptides: proteins
                .into_iter()
                .zip(decoys)
                .map(|(prots, decoy)| Peptide {
                    proteins: prots.into_iter().map(|s| Arc::from(s.as_ref())).collect(),
                    decoy: *decoy,
                    ..Default::default()
                })
                .collect(),
            decoy_tag: "DECOY_".to_string(),
            generate_decoys: false,
            ..IndexedDatabase::default()
        };
        (db, features)
    }

    #[test]
    fn test_protein_grouping_expected_groups() {
        let (proteins, decoys, q_vals) = get_data();
        // let protein_slices = proteins.iter().map(|v| v.as_slice()).collect::<Vec<_>>();
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_proteingroups(&db, &mut features, true, Some(0.01));
        let expected = vec![
            "protein_7",
            "protein_4/protein_9;protein_6",
            "protein_1",
            "protein_1",
            "protein_7",
            "protein_6",
            "protein_1",
            "protein_1",
            "protein_1",
            "protein_4/protein_9",
        ]
        .iter()
        .enumerate()
        .map(|(i, v)| (i + 1, v.to_string()))
        .collect::<Vec<_>>();
        let actual = features
            .into_iter()
            .enumerate()
            .map(|(i, v)| (i + 1, v.proteingroups.unwrap()))
            .collect::<Vec<_>>();
        // Compare actual and expected
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_eq!(a, e);
        }
    }
}
