// This module implements a protein grouping algorihm based on the IDPicker (1) algorithm with 
// extensions from the "Picked Protein FDR approach" (2). The Python implementation (3) of 
// CsoDIAq (4) has been used as template and for testing of the IDPicker approach. This function
// only implements IDPicker, the "rescued subset grouping (rsG)", discarding of shared peptides, 
// picked FDR are implemented as part of the core Sage FDR routines.
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

use crate::database::IndexedDatabase;
use crate::scoring::Feature;
use itertools::Itertools;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::{IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use rayon::prelude::{IntoParallelRefIterator, ParallelSliceMut};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;
use std::time::Instant;
use log::info;

pub fn generate_proteingroups(db: &IndexedDatabase, features: &mut [Feature]) {
    let gp1_filter = |label: i32, peptide_q: f32| label != -1 && peptide_q < 0.01;

    let gp2_filter = |label: i32, _: f32| label != -1;

    let pep_proteins_pg1 = get_peptide_protein_map(&db, features, gp1_filter);

    let pep_proteins_pg2 = get_peptide_protein_map(&db, features, gp2_filter);

    let protein_map_pg1 = get_proteingroups(pep_proteins_pg1);

    let protein_map_pg2 = get_proteingroups(pep_proteins_pg2);

    features.par_iter_mut().for_each(|feat| {
        let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);

        let array_proteins = proteins.split(";").collect::<Vec<_>>();

        let proteingroups: HashSet<_> = array_proteins
            .iter()
            .map(|&each_protein| {
                // Check if the protein exists in the IDpicker map
                if protein_map_pg1.contains_key(each_protein) {
                    // tier-1 proteins groups
                    protein_map_pg1.get(each_protein).unwrap().join("/")
                } else if protein_map_pg2.contains_key(each_protein) {
                    // tier-2 proteins groups
                    protein_map_pg2.get(each_protein).unwrap().join("/")
                } else {
                    // If the protein is not found in either tires, return an error
                    // panic!("{}",format!("Protein {} not found in protein group!!", each_protein));
                    each_protein.to_string()
                }
            })
            .collect();

        feat.proteingroups = Some(proteingroups.iter().sorted().join(";"));
        feat.num_proteingroups = proteingroups.len() as i32;
    });
}

fn get_peptide_protein_map(
    db: &&IndexedDatabase,
    features: &mut [Feature],
    gp1_filter: fn(i32, f32) -> bool,
) -> Vec<(Arc<String>, Arc<String>)> {
    features
        .iter()
        .filter(|feature: &&Feature| gp1_filter(feature.label, feature.peptide_q))
        .flat_map(|feature| {
            let idx = feature.peptide_idx;
            let proteins = &db[feature.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
            let peptide_proteins = proteins
                .split(";")
                .map(|s| (Arc::new(idx.0.to_string()), Arc::new(s.to_string())))
                .collect_vec();
            peptide_proteins
        })
        .collect_vec()
}

fn get_proteingroups(
    pep_proteins: Vec<(Arc<String>, Arc<String>)>,
) -> FxHashMap<String, Vec<String>> {
    info!("- number of records: {}", pep_proteins.len());
    let start = Instant::now();

    let start_meth = Instant::now();
    let group_by_peptides = get_group_by_peptides(pep_proteins);
    let end_time = Instant::now() - start_meth;
    info!("-  get_group_by_peptides: {:8} ms", end_time.as_millis());

    let start_meth = Instant::now();
    let grouped_peptides_and_grouped_proteins =
        get_grouped_peptides_and_proteins(group_by_peptides);
    let end_time = Instant::now() - start_meth;
    info!("-  get_grouped_peptides_and_proteins: {:8} ms", end_time.as_millis());

    let start_meth = Instant::now();
    let final_cluster = separate_into_clusters(grouped_peptides_and_grouped_proteins);
    let end_time = Instant::now() - start_meth;
    info!("-  separate_into_clusters: {:8} ms", end_time.as_millis());

    let final_proteins: Vec<(String, Vec<String>)> = final_cluster
        .into_par_iter()
        .flat_map(|cm| reduce_cluster(cm.1))
        .collect();
    let protein_map: FxHashMap<_, _> = final_proteins.into_par_iter().collect();

    let protein_map_process_time = Instant::now() - start;
    info!("- proteingroup: {:8} ms", protein_map_process_time.as_millis());
    protein_map
}

pub fn get_group_by_peptides(
    data: Vec<(Arc<String>, Arc<String>)>,
) -> FxHashSet<(Vec<Arc<String>>, Arc<String>)> {
    // Phase 1: Group proteins by peptide
    let group_by_peptides: FxHashMap<Arc<String>, Vec<Arc<String>>> = {
        data.par_iter()
            .fold(
                || FxHashMap::default(),
                |mut acc: FxHashMap<Arc<String>, Vec<Arc<String>>>, (peptide, protein)| {
                    acc.entry(peptide.clone()).or_default().push(protein.clone());
                    acc
                },
            )
            .reduce(FxHashMap::default, |mut acc, map| {
                for (key, val) in map {
                    acc.entry(key).or_default().extend(val);
                }
                acc
            })
    };
    

    // Phase 2: Group peptides by identical protein sets
    let group_by_protein: FxHashMap<Arc<[Arc<String>]>, Vec<Arc<String>>> = {
        group_by_peptides
            .into_par_iter()
            .map(|(peptide, mut proteins)| {
                proteins.sort(); // sort to make protein sets comparable
                let proteins_arc: Arc<[Arc<String>]> = Arc::from(proteins);
                (proteins_arc, peptide)
            })
            .fold(
                || FxHashMap::default(),
                |mut acc: FxHashMap<Arc<[Arc<String>]>, Vec<Arc<String>>>, (proteins, peptide)| {
                    acc.entry(proteins).or_default().push(peptide);
                    acc
                },
            )
            .reduce(FxHashMap::default, |mut acc, map| {
                for (key, val) in map {
                    acc.entry(key).or_default().extend(val);
                }
                acc
            })
    };
    
    // Phase 3: Flatten into peptide → peptide group map
    let peptide_to_group: FxHashMap<Arc<String>, Arc<[Arc<String>]>> = {
        group_by_protein
            .par_iter()
            .flat_map(|(protein_group, peptides)| {
                peptides
                    .par_iter()
                    .map(|p| (p.clone(), protein_group.clone()))
            })
            .collect()
    };
    
    // Phase 4: Construct final output set
    let result: FxHashSet<_> = data
        .into_par_iter()
        .map(|(peptide, protein)| {
            let group = peptide_to_group
                .get(&peptide)
                .expect("peptide group must exist");
            (group.to_vec(), protein)
        })
        .collect();
    
    result
}

pub fn get_grouped_peptides_and_proteins(
    grouped_peptides: FxHashSet<(Vec<Arc<String>>, Arc<String>)>,
) -> FxHashSet<(Vec<Arc<String>>, Vec<Arc<String>>)> {
    let group_by_proteins = grouped_peptides
        .par_iter()
        .fold(
            || FxHashMap::default(),
            |mut acc: FxHashMap<Arc<String>, Vec<Vec<Arc<String>>>>, x| {
                acc.entry(x.1.clone()).or_default().push(x.0.clone());
                acc
            },
        )
        .reduce(
            || FxHashMap::default(),
            |mut acc, map| {
                for (key, value) in map {
                    acc.entry(key).or_default().extend(value);
                }
                acc
            },
        );

    let group_by_peptide_group = group_by_proteins
        .par_iter()
        .fold(
            || FxHashMap::default(),
            |mut acc: FxHashMap<Vec<Vec<Arc<String>>>, Vec<Arc<String>>>, (key, val)| {
                let mut sorted_val = val.clone();
                sorted_val.sort();
                acc.entry(sorted_val).or_default().push(key.clone());
                acc
            },
        )
        .reduce(
            || FxHashMap::default(),
            |mut acc, map| {
                for (key, value) in map {
                    acc.entry(key).or_default().extend(value);
                }
                acc
            },
        );

    let map_group_proteins: FxHashMap<_, _> = group_by_peptide_group
        .values()
        .cloned()
        .collect::<Vec<_>>()
        .into_par_iter()
        .flat_map(|x1| {
            x1.iter()
                .map(|x2| (x2.clone(), x1.clone()))
                .collect::<Vec<_>>()
        })
        .collect();

    let grouped_peptides_and_proteins: FxHashSet<_> = grouped_peptides
        .par_iter()
        .map(|(peptide, protein)| (peptide.clone(), map_group_proteins[protein].clone()))
        .collect();

    grouped_peptides_and_proteins
}

pub fn separate_into_clusters(
    data: FxHashSet<(Vec<Arc<String>>, Vec<Arc<String>>)>,
) -> Vec<(usize, Vec<(Vec<Arc<String>>, Vec<Arc<String>>)>)> {
    // Group peptides by their protein groups first
    let group_peptides_proteins = group_peptides_by_proteins(&data);

    // Create sorted peptide and protein groups
    let (g_peps, g_proteins) = create_sorted_groups(&group_peptides_proteins);

    // Create protein sets for efficient lookup
    let g_proteins_sets = create_protein_sets(&g_proteins);

    // Build clusters using protein connectivity
    let mut cluster = build_clusters(&g_peps, &g_proteins, &g_proteins_sets);
    cluster.par_sort_unstable();
    cluster.dedup();

    // Create final cluster mapping
    let group_peptides_proteins = create_sorted_peptide_protein_map(group_peptides_proteins);
    let mut cluster_mapping = create_cluster_mapping(data, &cluster, &group_peptides_proteins);

    // Merge clusters with same index
    merge_clusters(&mut cluster_mapping);

    cluster_mapping
}

fn group_peptides_by_proteins(
    data: &FxHashSet<(Vec<Arc<String>>, Vec<Arc<String>>)>
) -> FxHashMap<Vec<Arc<String>>, Vec<Vec<Arc<String>>>> {
    let data_len = data.len();
    data.par_iter()
        .fold(
            || FxHashMap::with_capacity_and_hasher(data_len, Default::default()),
            |mut acc: FxHashMap<Vec<Arc<String>>, Vec<Vec<Arc<String>>>>, x| {
                acc.entry(x.0.clone()).or_default().push(x.1.clone());
                acc
            },
        )
        .reduce(
            || FxHashMap::with_capacity_and_hasher(data_len, Default::default()),
            |mut acc, map| {
                for (key, value) in map {
                    acc.entry(key).or_default().extend(value);
                }
                acc
            },
        )
}

fn create_sorted_groups(
    group_peptides_proteins: &FxHashMap<Vec<Arc<String>>, Vec<Vec<Arc<String>>>>
) -> (Vec<Vec<Arc<String>>>, Vec<Vec<Vec<Arc<String>>>>) {
    let gpp_len = group_peptides_proteins.len();
    let mut g_peps = Vec::with_capacity(gpp_len);
    let mut g_proteins = Vec::with_capacity(gpp_len);
    
    for (pep, prot) in group_peptides_proteins.iter() {
        let mut pep = pep.clone();
        pep.sort();
        g_peps.push(pep);
        g_proteins.push(prot.clone());
    }
    
    (g_peps, g_proteins)
}

fn create_protein_sets(
    g_proteins: &[Vec<Vec<Arc<String>>>]
) -> Vec<FxHashSet<Arc<String>>> {
    g_proteins
        .par_iter()
        .map(|proteins| {
            let capacity = proteins.iter().map(|v| v.len()).sum();
            let mut set = FxHashSet::with_capacity_and_hasher(capacity, Default::default());
            proteins.iter().flatten().for_each(|x| { set.insert(x.clone()); });
            set
        })
        .collect()
}

fn build_clusters<'a>(
    g_peps: &'a[Vec<Arc<String>>],
    g_proteins: &'a[Vec<Vec<Arc<String>>>],
    g_proteins_sets: &'a[FxHashSet<Arc<String>>]
) -> Vec<Vec<&'a Vec<Arc<String>>>> {
    g_peps
        .par_iter()
        .enumerate()
        .map(|(i, _)| {
            let mut mini_cluster = HashSet::with_capacity(16);
            let mut in_proteins_set = g_proteins_sets[i].clone();
            let mut new_proteins = HashSet::with_capacity(in_proteins_set.len());

            loop {
                new_proteins.clear();
                let prev_len = mini_cluster.len();
                
                check_protein_group(g_proteins, g_proteins_sets, &in_proteins_set)
                    .for_each(|pg| {
                        mini_cluster.insert(&g_peps[pg.0]);
                        pg.1.iter().flatten().for_each(|x| {
                            new_proteins.insert(Arc::clone(x));
                        });
                    });

                if new_proteins.is_empty() || prev_len == mini_cluster.len() {
                    break;
                }
                
                in_proteins_set.extend(new_proteins.drain());
            }
            
            let mut result = Vec::with_capacity(mini_cluster.len());
            result.extend(mini_cluster);
            result.sort_unstable();
            result
        })
        .collect()
}

fn create_sorted_peptide_protein_map(
    group_peptides_proteins: FxHashMap<Vec<Arc<String>>, Vec<Vec<Arc<String>>>>
) -> FxHashMap<Vec<Arc<String>>, Vec<Vec<Arc<String>>>> {
    group_peptides_proteins
        .into_par_iter()
        .map(|(mut k, v)| {
            k.sort();
            (k, v)
        })
        .collect()
}

fn create_cluster_mapping(
    data: FxHashSet<(Vec<Arc<String>>, Vec<Arc<String>>)>,
    cluster: &[Vec<&Vec<Arc<String>>>],
    group_peptides_proteins: &FxHashMap<Vec<Arc<String>>, Vec<Vec<Arc<String>>>>
) -> Vec<(usize, Vec<(Vec<Arc<String>>, Vec<Arc<String>>)>)> {
    data.into_par_iter()
        .map(|(mut peps, _)| {
            peps.sort();
            let cluster_idx = cluster
                .iter()
                .enumerate()
                .find(|x| x.1.contains(&&peps))
                .unwrap()
                .0;
                
            let proteins = group_peptides_proteins.get(&peps).unwrap();
            let protein_matches = proteins
                .iter()
                .map(|protein| (peps.clone(), protein.clone()))
                .collect::<Vec<_>>();

            (cluster_idx, protein_matches)
        })
        .collect()
}

fn merge_clusters(
    cluster_mapping: &mut Vec<(usize, Vec<(Vec<Arc<String>>, Vec<Arc<String>>)>)>
) {
    cluster_mapping.par_sort_unstable();
    cluster_mapping.dedup();

    cluster_mapping.dedup_by(|remove, keep| {
        if remove.0 == keep.0 {
            keep.1.extend(remove.1.iter().cloned());
            true
        } else {
            false
        }
    });
}

fn check_protein_group<'a>(
    g_proteins: &'a [Vec<Vec<Arc<String>>>],
    g_proteins_set: &'a [FxHashSet<Arc<String>>],
    in_proteins_set: &'a FxHashSet<Arc<String>>,
) -> impl Iterator<Item = (usize, &'a Vec<Vec<Arc<String>>>)> + 'a {
    g_proteins
        .iter()
        .zip(g_proteins_set)
        .enumerate()
        .filter_map(|(i, (proteins, proteins_set))| {
            if in_proteins_set.iter().any(|x| proteins_set.contains(x)) {
                Some((i, proteins))
            } else {
                None
            }
        })
}

pub fn reduce_cluster(
    data: Vec<(Vec<Arc<String>>, Vec<Arc<String>>)>,
) -> Vec<(String, Vec<String>)> {
    let mut proteins = Vec::new();
    let mut used_peptides: HashSet<Vec<Arc<String>>> = HashSet::new();
    let total_peptides: HashSet<_> = data.iter().map(|(p, _)| p).collect();

    // Compute overall protein scores
    let mut overall_score_dict = HashMap::new();

    data.iter().for_each(|(_, proteins)| {
        *overall_score_dict.entry(proteins.clone()).or_insert(0) += 1;
    });

    while used_peptides.len() != total_peptides.len() {
        let mut prot_dict: HashMap<Vec<Arc<String>>, HashSet<Vec<Arc<String>>>> = HashMap::new();

        data.iter()
            .filter(|(peptides, _)| !used_peptides.contains(peptides))
            .for_each(|(peptides, proteins)| {
                prot_dict
                    .entry(proteins.clone())
                    .or_default()
                    .insert(peptides.clone());
            });

        let score_dict: BTreeMap<Vec<Arc<String>>, (usize, usize)> = prot_dict
            .iter()
            .map(|(key, peptides)| {
                (
                    key.clone(),
                    (peptides.len(), *overall_score_dict.get(key).unwrap_or(&0)),
                )
            })
            .collect();

        // Sort proteins based on (peptide count, overall score)
        let sorted_keys = score_dict
            .into_iter()
            .sorted_by(|a, b| b.1.cmp(&a.1))
            .collect::<Vec<_>>();

        if let Some((top_key, _)) = sorted_keys.first() {
            // let protein = top_key[0].clone();
            // let protein_str = format!("{}/{}", top_key.len(), top_key.iter().join("/"));
            let arry_protein = top_key
                .into_iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>();
            // println!("{},{}", protein, protein_str);
            // proteins.push((protein.to_string(), protein_str));
            proteins.extend(
                top_key
                    .iter()
                    .map(|x| (x.to_string(), arry_protein.clone()))
                    .collect::<Vec<_>>(),
            );

            if let Some(peptides) = prot_dict.get(top_key).cloned() {
                used_peptides.extend(peptides);
            }
        }
    }

    proteins
}

#[cfg(test)]
mod test {
    use crate::idpicker::get_proteingroups;
    use std::sync::Arc;

    #[test]
    fn test_valid_proteins() {
        let test_data = vec![
            ("peptide_1".to_string(), "protein_7".to_string()),
            ("peptide_2".to_string(), "protein_4".to_string()),
            ("peptide_2".to_string(), "protein_6".to_string()),
            ("peptide_2".to_string(), "protein_9".to_string()),
            ("peptide_3".to_string(), "protein_1".to_string()),
            ("peptide_4".to_string(), "protein_1".to_string()),
            ("peptide_4".to_string(), "protein_5".to_string()),
            ("peptide_5".to_string(), "protein_7".to_string()),
            ("peptide_6".to_string(), "protein_3".to_string()),
            ("peptide_6".to_string(), "protein_6".to_string()),
            ("peptide_7".to_string(), "protein_1".to_string()),
            ("peptide_8".to_string(), "protein_1".to_string()),
            ("peptide_8".to_string(), "protein_2".to_string()),
            ("peptide_8".to_string(), "protein_5".to_string()),
            ("peptide_8".to_string(), "protein_8".to_string()),
            ("peptide_9".to_string(), "protein_1".to_string()),
            ("peptide_10".to_string(), "protein_4".to_string()),
            ("peptide_10".to_string(), "protein_9".to_string()),
        ];

        let data: Vec<_> = test_data
            .into_iter()
            .map(|(k, v)| (Arc::new(k.to_string()), Arc::new(v.to_string())))
            .collect();

        let protein_map = get_proteingroups(data);

        for (i, j) in protein_map.iter() {
            println!("Protein  ->{:?} ,  Protein Group-> {:?} ", i, j);
        }

        assert_eq!(5, protein_map.len());
    }
}
