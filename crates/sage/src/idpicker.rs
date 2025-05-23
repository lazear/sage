use crate::database::{IndexedDatabase, PeptideIx};
use crate::peptide::Peptide;
use crate::scoring::Feature;
use itertools::Itertools;
use log::info;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::sync::Arc;
use std::time::Instant;

const PROTEIN_COVER_TYPE: ProteinCoverType = ProteinCoverType::Slim;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct ProteinIx(u32);

#[derive(Debug, PartialEq, Eq, Hash, Clone, Default)]
struct ProteinGroup(Vec<ProteinIx>);

impl ProteinGroup {
    fn to_string(
        &self,
        protein_map: &Vec<(Arc<String>, bool)>,
        decoy_tag: &str,
        generate_decoys: bool,
    ) -> String {
        self.0
            .iter()
            .map(|prot_ix| {
                let (prot, decoy) = &protein_map[prot_ix.0 as usize];
                if *decoy & generate_decoys {
                    format!("{}{}", decoy_tag, prot)
                } else {
                    prot.to_string()
                }
            })
            .sorted()
            .join("/")
    }
}

#[derive(Debug, Default)]
struct ProteinMapping {
    protein_groups: Vec<ProteinGroup>,
    connections: Vec<(usize, usize)>,
    prot_name_map: HashMap<(Arc<String>, bool), ProteinIx>,
    protein_cover: Vec<usize>,
    peptide_count: usize,
}

pub enum ProteinCoverType {
    All,
    Slim,
}

impl ProteinMapping {
    fn new(
        db: &IndexedDatabase,
        features: &mut [Feature],
        predicate: fn(&&Feature) -> bool,
        protein_cover_type: ProteinCoverType,
    ) -> Self {
        let mut mapping = Self::default();
        let time = Instant::now();
        let peps = Self::get_peps(features, predicate);
        mapping.peptide_count = peps.len();
        info!(
            "-  found {} unique peptides in {:?}ms",
            peps.len(),
            time.elapsed().as_millis()
        );
        let time = Instant::now();
        let meta_peptides = mapping.set_proteins_and_get_metapeptides(peps, &db);
        info!(
            "-  found {} meta_peptides in {:?}ms",
            meta_peptides.len(),
            time.elapsed().as_millis()
        );
        let time = Instant::now();
        mapping.set_protein_groups_and_connections(meta_peptides);
        info!(
            "-  found {} protein_groups in {:?}ms",
            mapping.protein_groups.len(),
            time.elapsed().as_millis()
        );
        let time = Instant::now();
        match protein_cover_type {
            ProteinCoverType::All => {
                mapping.protein_cover = (0..mapping.protein_groups.len()).collect()
            }
            ProteinCoverType::Slim => mapping.find_greedy_protein_cover(),
        }
        info!(
            "-  found cover of {} metaproteins in {:?}ms",
            mapping.protein_cover.len(),
            time.elapsed().as_millis()
        );
        mapping
    }

    fn get_peps(features: &[Feature], predicate: fn(&&Feature) -> bool) -> HashSet<PeptideIx> {
        features
            .par_iter()
            .filter(predicate)
            .map(|feature| feature.peptide_idx)
            .collect::<HashSet<_>>()
    }

    fn set_proteins_and_get_metapeptides(
        &mut self,
        peps: HashSet<PeptideIx>,
        db: &&IndexedDatabase,
    ) -> HashSet<Vec<ProteinIx>> {
        peps.into_iter()
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

    fn set_protein_groups_and_connections<'a>(&mut self, meta_peptides: HashSet<Vec<ProteinIx>>) {
        let mut protein_map = HashMap::new();
        meta_peptides
            .into_iter()
            .enumerate()
            .for_each(|(i, prot_ids)| {
                prot_ids.into_iter().for_each(|prot_id| {
                    let entry = protein_map.entry(prot_id).or_insert_with(|| vec![]);
                    entry.push(i);
                });
            });
        let mut mapping = HashMap::new();
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
        mapping.into_iter().enumerate().for_each(
            |(meta_protein_index, (meta_peptides, meta_protein))| {
                protein_groups.push(meta_protein);
                meta_peptides.iter().for_each(|&meta_peptide_index| {
                    connections.push((meta_protein_index, meta_peptide_index));
                });
            },
        );
        self.protein_groups = protein_groups;
        self.connections = connections;
    }

    fn find_greedy_protein_cover(&mut self) {
        if self.connections.is_empty() {
            return;
        }
        let mut greedy_cover = vec![false; self.protein_groups.len()];
        let mut remaining_pep_counts = vec![0; self.peptide_count];
        self.connections.iter().for_each(|(_, pep_index)| {
            remaining_pep_counts[*pep_index] += 1;
        });
        let mut connection_count = 0;
        loop {
            while connection_count != self.connections.len() {
                connection_count = self.connections.len();
                self.connections.retain(|(prot_index, pep_index)| {
                    if remaining_pep_counts[*pep_index] == 1 {
                        greedy_cover[*prot_index] = true;
                        remaining_pep_counts[*pep_index] -= 1;
                        false
                    } else {
                        true
                    }
                });
                self.connections.retain(|(prot_index, pep_index)| {
                    if greedy_cover[*prot_index] {
                        remaining_pep_counts[*pep_index] -= 1;
                        false
                    } else {
                        true
                    }
                })
            }
            if self.connections.len() == 0 {
                break;
            }
            let mut remaining_prot_counts = HashMap::new();
            self.connections.iter().for_each(|(prot_index, _)| {
                let entry = remaining_prot_counts.entry(*prot_index).or_insert(0);
                *entry += 1;
            });
            let &max_index = remaining_prot_counts.iter().max_by_key(|x| x.1).unwrap().0;
            greedy_cover[max_index] = true;
            self.connections.retain(|(prot_index, pep_index)| {
                if *prot_index == max_index {
                    remaining_pep_counts[*pep_index] -= 1;
                    false
                } else {
                    true
                }
            })
        }
        self.protein_cover = greedy_cover
            .into_iter()
            .enumerate()
            .filter(|(_, x)| *x)
            .map(|(i, _)| i)
            .collect();
    }
}

struct ProteinGroupMap {
    protein_groups: Vec<ProteinGroup>,
    proteins: Vec<(Arc<String>, bool)>,
    protein_group_map: HashMap<(Arc<String>, bool), Vec<usize>>,
}

impl ProteinGroupMap {
    fn new(mapping: ProteinMapping) -> Self {
        let proteins = mapping
            .prot_name_map
            .into_iter()
            .sorted_by_key(|(_, v)| v.0)
            .map(|(k, _)| k.clone())
            .collect::<Vec<_>>();
        let mut protein_group_map = HashMap::new();
        mapping.protein_cover.into_iter().for_each(|i| {
            let protein_group = &mapping.protein_groups[i];
            protein_group.0.iter().for_each(|&prot_id| {
                let (name, decoy) = &proteins[prot_id.0 as usize];
                let entry = protein_group_map
                    .entry((name.clone(), *decoy))
                    .or_insert_with(|| vec![]);
                entry.push(i);
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
                    .map(|&i| &self.protein_groups[i])
            })
            .collect::<HashSet<_>>();
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

pub fn generate_proteingroups(db: &IndexedDatabase, features: &mut [Feature]) {
    let time = Instant::now();
    info!("Protein grouping with {} features", features.len());
    let protein_map_pg1 = {
        let map = ProteinMapping::new(
            db,
            features,
            |f| (f.label != -1) && (f.peptide_q < 0.01),
            PROTEIN_COVER_TYPE,
        );
        ProteinGroupMap::new(map)
    };
    let protein_map_pg2 = {
        let map = ProteinMapping::new(db, features, |f| f.label != -1, PROTEIN_COVER_TYPE);
        ProteinGroupMap::new(map)
    };
    features.par_iter_mut().for_each(|feat| {
        let pep = &db[feat.peptide_idx];
        let mut proteingroups = protein_map_pg1.get_protein_group_string(pep, db);
        if proteingroups.is_empty() {
            proteingroups = protein_map_pg2.get_protein_group_string(pep, db);
        }
        if proteingroups.is_empty() {
            proteingroups =
                proteins(&pep.proteins, &db.decoy_tag, db.generate_decoys, pep.decoy).collect()
        }
        feat.proteingroups = Some(proteingroups.iter().sorted().join(";"));
        feat.num_proteingroups = proteingroups.len() as i32;
    });
    info!(
        "Grouped and inferred proteins in {:?}ms",
        time.elapsed().as_millis()
    );
}

pub fn proteins<'a>(
    proteins: &'a [Arc<String>],
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
