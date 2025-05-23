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
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct ProteinIx(pub u32);

#[derive(Debug, PartialEq, Eq, Hash)]
struct MetaPeptide {
    pub peps: Vec<PeptideIx>,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
struct MetaProtein {
    pub prots: Vec<ProteinIx>,
}

#[derive(Debug, Default)]
struct ProteinMapping {
    pub meta_peptides: Vec<MetaPeptide>,
    pub meta_proteins: Vec<MetaProtein>,
    pub connections: Vec<(usize, usize)>,
    pub prot_name_map: HashMap<(Arc<String>, bool), ProteinIx>,
}

impl ProteinMapping {
    fn new(
        db: &IndexedDatabase,
        features: &mut [Feature],
        predicate: fn(&&Feature) -> bool,
    ) -> Self {
        let time = Instant::now();
        let mut mapping = Self::default();
        let meta_peptides = mapping.set_proteins_and_get_metapeptides(&features, &db, predicate);
        mapping.find_mapping(meta_peptides);
        info!(
            "-  found {} metaproteins in {:?}ms",
            mapping.meta_proteins.len(),
            time.elapsed().as_millis()
        );
        mapping
    }

    fn set_proteins_and_get_metapeptides(
        &mut self,
        features: &[Feature],
        db: &&IndexedDatabase,
        predicate: fn(&&Feature) -> bool,
    ) -> HashMap<Vec<ProteinIx>, MetaPeptide> {
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
                            self.get_or_insert_prot(prot)
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

    fn get_or_insert_prot(&mut self, prot: (Arc<String>, bool)) -> ProteinIx {
        match self.prot_name_map.get(&prot) {
            Some(&prot_id) => prot_id,
            None => {
                let prot_id = ProteinIx(self.prot_name_map.len() as u32);
                self.prot_name_map.insert(prot, prot_id);
                prot_id
            }
        }
    }

    fn find_mapping<'a>(&mut self, meta_peptides: HashMap<Vec<ProteinIx>, MetaPeptide>) {
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
        mapping.into_iter().enumerate().for_each(
            |(meta_protein_index, (meta_peptides, meta_protein))| {
                meta_proteins.push(meta_protein);
                meta_peptides.iter().for_each(|&meta_peptide_index| {
                    connections.push((meta_protein_index, meta_peptide_index));
                });
            },
        );
        self.meta_peptides = meta_peptides;
        self.meta_proteins = meta_proteins;
        self.connections = connections;
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
    let mapping = ProteinMapping::new(db, features, |f| (f.label != -1) && (f.peptide_q < 0.01));
    let protein_cover = find_greedy_protein_cover(mapping.connections.clone());
    let protein_map_pg1 = mapping.get_protein_map(protein_cover);
    let mapping = ProteinMapping::new(db, features, |f| f.label != -1);
    let protein_cover = find_greedy_protein_cover(mapping.connections.clone());
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

fn find_greedy_protein_cover(mut connections: Vec<(usize, usize)>) -> Vec<usize> {
    let time = Instant::now();
    let mut greedy_cover = HashSet::new();
    let mut remaining_pep_counts = HashMap::new();
    connections.iter().for_each(|(_, pep_index)| {
        let entry = remaining_pep_counts.entry(*pep_index).or_insert(0);
        *entry += 1;
    });
    let mut connection_count = 0;
    loop {
        while connection_count != connections.len() {
            connection_count = connections.len();
            connections.retain(|(prot_index, pep_index)| {
                if remaining_pep_counts[pep_index] == 1 {
                    greedy_cover.insert(*prot_index);
                    remaining_pep_counts.remove(pep_index);
                    false
                } else {
                    true
                }
            });
            connections.retain(|(prot_index, pep_index)| {
                if greedy_cover.contains(prot_index) {
                    let entry = remaining_pep_counts.get_mut(pep_index).unwrap();
                    *entry -= 1;
                    false
                } else {
                    true
                }
            })
        }
        if connections.len() == 0 {
            break;
        }
        let mut remaining_prot_counts = HashMap::new();
        connections.iter().for_each(|(prot_index, _)| {
            let entry = remaining_prot_counts.entry(*prot_index).or_insert(0);
            *entry += 1;
        });
        let &max_index = remaining_prot_counts.iter().max_by_key(|x| x.1).unwrap().0;
        greedy_cover.insert(max_index);
        connections.retain(|(prot_index, pep_index)| {
            if *prot_index == max_index {
                let entry = remaining_pep_counts.get_mut(pep_index).unwrap();
                *entry -= 1;
                false
            } else {
                true
            }
        })
    }
    info!(
        "-  found cover of {} metaproteins in {:?}ms",
        greedy_cover.len(),
        time.elapsed().as_millis()
    );
    greedy_cover.into_iter().collect()
}

// pub fn reduce_cluster(data: Vec<(&Vec<u32>, &Vec<Arc<String>>)>) -> Vec<(String, Vec<String>)> {
//     // println!("data {:?}",data);
//     let mut proteins = Vec::new();
//     let mut used_peptides: HashSet<&Vec<u32>> = HashSet::new();
//     let total_peptides: HashSet<_> = data.iter().map(|x| x.0).collect();

//     let mut overall_score_dict = HashMap::new();
//     data.iter().for_each(|(_, proteins)| {
//         *overall_score_dict.entry(proteins).or_insert(0) += 1;
//     });

//     while used_peptides.len() != total_peptides.len() {
//         let mut prot_dict: HashMap<&Vec<Arc<String>>, HashSet<&Vec<u32>>> = HashMap::new();

//         data.iter()
//             .filter(|(peptides, _)| !used_peptides.contains(peptides))
//             .for_each(|(peptides, proteins)| {
//                 prot_dict.entry(proteins).or_default().insert(peptides);
//             });

//         let score_dict: Vec<(&Vec<Arc<String>>, (usize, i32))> = prot_dict
//             .iter()
//             .map(|(key, peptides)| {
//                 (
//                     *key,
//                     (peptides.len(), *overall_score_dict.get(key).unwrap_or(&0)),
//                 )
//             })
//             .collect::<Vec<_>>();

//         // Sort proteins based on (peptide count, overall score)
//         let sorted_keys = score_dict
//             .into_iter()
//             .sorted_by(|a, b| b.1.cmp(&a.1))
//             .collect::<Vec<_>>();

//         if let Some((top_key, _)) = sorted_keys.first() {
//             let arry_protein = top_key
//                 .into_iter()
//                 .map(|x| x.to_string())
//                 .collect::<Vec<_>>();
//             proteins.extend(
//                 top_key
//                     .iter()
//                     .map(|x| (x.to_string(), arry_protein.clone()))
//                     .collect::<Vec<_>>(),
//             );

//             if let Some(peptides) = prot_dict.get(top_key).cloned() {
//                 used_peptides.extend(peptides);
//             }
//         }
//     }

//     proteins
// }

// fn get_proteingroups(pep_proteins: Vec<(PeptideIx, ProteinIx)>) -> FxHashMap<String, Vec<String>> {}

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
