// This module implements a protein grouping algorihm based on the IDPicker (1) algorithm with
// extensions from the "Picked Protein FDR approach" (2). The Python implementation (3) of
// CsoDIAq (4) has been used as template and for testing of the IDPicker approach. Most functions
// are based on IDPicker, only generate_proteingroups() represents the "rescued subset grouping
// (rsG)", approach of (2). Discarding of shared peptides and picked FDR are implemented as part
// of the core Sage FDR routines.
//
// 1. Zhang, B., Chambers, M. C., & Tabb, D. L. (2007). Proteomic parsimony through
// bipartite graph analysis improves accuracy and transparency. Journal of proteome research,
// 6(9), 3549-3557. https://doi.org/10.1021/pr070230d
//
// 2. The, M., Samaras, P., Kuster, B., & Wilhelm, M. (2022). Reanalysis of ProteomicsDB using
// an accurate, sensitive, and scalable false discovery rate estimation approach for protein
// groups. Molecular & Cellular Proteomics, 21(12), 100437. https://doi.org/10.1016/j.mcpro.2022.100437
//
// 3. https://github.com/dg310012/CsoDIAq/blob/68abaa713eb719b488967cb34a876a71657827bd/idpicker.py
//
// 4. Cranney, C. W., & Meyer, J. G. (2021). CsoDIAq software for direct infusion shotgun proteome
// analysis. Analytical Chemistry, 93(36), 12312-12319. https://doi.org/10.1021/acs.analchem.1c02021

use crate::database::{IndexedDatabase, PeptideIx};
use crate::scoring::Feature;
use itertools::Itertools;
use log::info;
use rayon::{prelude::*, vec};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct ProtId(pub u32);

#[derive(Debug, PartialEq, Eq, Hash)]
struct MetaPeptide {
    pub peps: Vec<PeptideIx>,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
struct MetaProtein {
    pub prots: Vec<ProtId>,
}

#[derive(Debug)]
struct Mapping {
    pub meta_peptides: Vec<MetaPeptide>,
    pub meta_proteins: Vec<MetaProtein>,
    pub connections: Vec<Vec<usize>>,
    pub prot_name_map: HashMap<(Arc<String>, bool), ProtId>,
}

impl Mapping {
    fn new(
        db: &IndexedDatabase,
        features: &mut [Feature],
        predicate: fn(&&Feature) -> bool,
    ) -> Self {
        let time = Instant::now();
        let mut prot_name_map: HashMap<(Arc<String>, bool), ProtId> = HashMap::new();
        let meta_peptides = find_meta_peptides(&features, &db, &mut prot_name_map, predicate);
        let mapping = find_mapping(meta_peptides, prot_name_map);
        info!(
            "-  found {} metaproteins in {:?}ms",
            mapping.meta_proteins.len(),
            time.elapsed().as_millis()
        );
        mapping
    }

    fn find_greedy_protein_cover(&self) -> Vec<usize> {
        let time = Instant::now();
        let mut remaining = self.meta_peptides.len();
        let mut usable_peptides = vec![1; remaining];
        let mut greedy_cover = vec![];
        let mut connections = self.connections.clone();
        while remaining > 0 {
            let current_connection_counts = connections
                .iter()
                .map(|c| c.iter().map(|&i| usable_peptides[i]).sum())
                .collect::<Vec<_>>();
            let (largest_prot_index, count) = current_connection_counts
                .iter()
                .enumerate()
                .max_by_key(|(_, count)| *count)
                .unwrap();
            greedy_cover.push(largest_prot_index);
            remaining -= count;
            connections[largest_prot_index]
                .drain(..)
                .for_each(|pep_index| {
                    usable_peptides[pep_index] = 0;
                });
        }
        info!(
            "-  found metaprotein cover of {} in {:?}ms",
            greedy_cover.len(),
            time.elapsed().as_millis()
        );
        greedy_cover
    }

    fn get_protein_map(
        &self,
        meta_protein_cover: Vec<usize>,
    ) -> HashMap<(Arc<String>, bool), String> {
        let reversed_name_map = self
            .prot_name_map
            .iter()
            .map(|(k, v)| (v, k))
            .collect::<HashMap<_, _>>();
        let mut protein_map = HashMap::new();
        meta_protein_cover.iter().for_each(|&i| {
            let meta_protein = &self.meta_proteins[i];
            let prot_ids = meta_protein.prots.clone();
            let group_name = prot_ids
                .iter()
                .map(|&prot_id| {
                    let (name, _) = reversed_name_map[&prot_id];
                    name
                })
                .join("/");
            prot_ids.iter().for_each(|&prot_id| {
                let (name, decoy) = reversed_name_map[&prot_id];
                protein_map.insert((name.clone(), *decoy), group_name.clone());
            });
        });
        protein_map
    }
}

pub fn generate_proteingroups(db: &IndexedDatabase, features: &mut [Feature]) {
    let mapping = Mapping::new(db, features, |f| (f.label != -1) && (f.peptide_q < 0.01));
    let protein_cover = mapping.find_greedy_protein_cover();
    let protein_map_pg1 = mapping.get_protein_map(protein_cover);
    let mapping = Mapping::new(db, features, |f| f.label != -1);
    let protein_cover = mapping.find_greedy_protein_cover();
    let protein_map_pg2 = mapping.get_protein_map(protein_cover);
    features.par_iter_mut().for_each(|feat| {
        let pep = &db[feat.peptide_idx];
        let proteingroups: HashSet<_> = pep
            .proteins
            .iter()
            .map(|protein| {
                let prot = (protein.clone(), pep.decoy);
                match protein_map_pg1.get(&prot) {
                    Some(protein_group) => protein_group.clone(),
                    None => {
                        match protein_map_pg2.get(&prot) {
                            Some(protein_group) => protein_group.clone(),
                            None => {
                                // If the protein is not found in the map, return the protein itself
                                protein.to_string()
                            }
                        }
                    }
                }
            })
            .collect();
        feat.proteingroups = Some(proteingroups.iter().sorted().join(";"));
        feat.num_proteingroups = proteingroups.len() as i32;
    });
}

fn find_meta_peptides(
    features: &[Feature],
    db: &&IndexedDatabase,
    prot_name_map: &mut HashMap<(Arc<String>, bool), ProtId>,
    predicate: fn(&&Feature) -> bool,
) -> HashMap<Vec<ProtId>, MetaPeptide> {
    let found_pep_ids: HashSet<PeptideIx> = HashSet::new();
    let mut mapping = HashMap::new();
    features.iter().filter(predicate).for_each(|feature| {
        let pep_id = feature.peptide_idx;
        let db_peptide = &db[pep_id];
        if !found_pep_ids.contains(&pep_id) {
            {
                let prot_ids = db_peptide
                    .proteins
                    .iter()
                    .map(|prot| {
                        let prot = (prot.clone(), db_peptide.decoy);
                        get_or_insert_prot(prot_name_map, prot)
                    })
                    .sorted()
                    .collect::<Vec<_>>();
                let entry = mapping
                    .entry(prot_ids)
                    .or_insert_with(|| MetaPeptide { peps: vec![] });
                entry.peps.push(pep_id);
            }
        }
    });
    mapping
}

fn get_or_insert_prot(
    prot_name_map: &mut HashMap<(Arc<String>, bool), ProtId>,
    prot: (Arc<String>, bool),
) -> ProtId {
    match prot_name_map.get(&prot) {
        Some(&prot_id) => prot_id,
        None => {
            let prot_id = ProtId(prot_name_map.len() as u32);
            prot_name_map.insert(prot, prot_id);
            prot_id
        }
    }
}

fn find_mapping<'a>(
    meta_peptides: HashMap<Vec<ProtId>, MetaPeptide>,
    prot_name_map: HashMap<(Arc<String>, bool), ProtId>,
) -> Mapping {
    let mut protein_map = HashMap::new();
    let meta_peptides = meta_peptides
        .into_iter()
        .enumerate()
        .map(|(i, (prot_ids, meta_peptide))| {
            prot_ids.iter().for_each(|&prot_id| {
                let entry = protein_map.entry(prot_id).or_insert_with(|| vec![]);
                entry.push(i);
            });
            meta_peptide
        })
        .collect::<Vec<_>>();
    let mut mapping = HashMap::new();
    protein_map
        .into_iter()
        .for_each(|(prot_id, meta_peptides)| {
            let entry = mapping
                .entry(meta_peptides)
                .or_insert_with(|| MetaProtein { prots: vec![] });
            entry.prots.push(prot_id);
        });
    let mut meta_proteins = vec![];
    let mut connections = vec![];
    mapping
        .into_iter()
        .for_each(|(meta_peptides, meta_protein)| {
            meta_proteins.push(meta_protein);
            connections.push(meta_peptides);
        });
    let mapping = Mapping {
        prot_name_map,
        meta_peptides,
        meta_proteins,
        connections,
    };
    mapping
}

// fn get_proteingroups(pep_proteins: Vec<(PeptideIx, ProtId)>) -> FxHashMap<String, Vec<String>> {}

// #[cfg(test)]
// mod test {
//     use crate::idpicker::get_proteingroups;
//     use rustc_hash::FxHashMap;
//     use std::sync::Arc;

//     #[test]
//     fn test_valid_proteins() {
//         let test_data = vec![
//             ("peptide_1".to_string(), "protein_7".to_string()),
//             ("peptide_2".to_string(), "protein_4".to_string()),
//             ("peptide_2".to_string(), "protein_6".to_string()),
//             ("peptide_2".to_string(), "protein_9".to_string()),
//             ("peptide_3".to_string(), "protein_1".to_string()),
//             ("peptide_4".to_string(), "protein_1".to_string()),
//             ("peptide_4".to_string(), "protein_5".to_string()),
//             ("peptide_5".to_string(), "protein_7".to_string()),
//             ("peptide_6".to_string(), "protein_3".to_string()),
//             ("peptide_6".to_string(), "protein_6".to_string()),
//             ("peptide_7".to_string(), "protein_1".to_string()),
//             ("peptide_8".to_string(), "protein_1".to_string()),
//             ("peptide_8".to_string(), "protein_2".to_string()),
//             ("peptide_8".to_string(), "protein_5".to_string()),
//             ("peptide_8".to_string(), "protein_8".to_string()),
//             ("peptide_9".to_string(), "protein_1".to_string()),
//             ("peptide_10".to_string(), "protein_4".to_string()),
//             ("peptide_10".to_string(), "protein_9".to_string()),
//         ];

//         let data: Vec<_> = test_data
//             .into_iter()
//             .map(|(k, v)| (Arc::new(k.to_string()), Arc::new(v.to_string())))
//             .collect();

//         let proteingroups = get_proteingroups(data);

//         let expected_proteingroups: FxHashMap<String, Vec<String>> = vec![
//             ("protein_9", vec!["protein_9", "protein_4"]),
//             ("protein_6", vec!["protein_6"]),
//             ("protein_1", vec!["protein_1"]),
//             ("protein_4", vec!["protein_9", "protein_4"]),
//             ("protein_7", vec!["protein_7"]),
//         ]
//         .into_iter()
//         .map(|(k, v)| {
//             (
//                 k.to_string(),
//                 v.into_iter().map(|s| s.to_string()).collect::<Vec<_>>(),
//             )
//         })
//         .collect();

//         assert_eq!(expected_proteingroups, proteingroups);
//     }
// }
