use crate::database::IndexedDatabase;
use crate::scoring::Feature;
use itertools::Itertools;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::{IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use rayon::prelude::{IntoParallelRefIterator, ParallelSliceMut};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

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
    let group_by_peptides = get_group_by_peptides(pep_proteins);
    let grouped_peptides_and_grouped_proteins =
        get_grouped_peptides_and_proteins(group_by_peptides);
    let final_cluster = separate_into_clusters(grouped_peptides_and_grouped_proteins);
    let final_proteins: Vec<(String, Vec<String>)> = final_cluster
        .into_par_iter()
        .flat_map(|cm| reduce_cluster(cm.1))
        .collect();
    let protein_map: FxHashMap<_, _> = final_proteins.into_par_iter().collect();
    protein_map
}

pub fn get_group_by_peptides(
    data: Vec<(Arc<String>, Arc<String>)>,
) -> FxHashSet<(Vec<Arc<String>>, Arc<String>)> {
    let group_by_peptides = &data
        .par_iter()
        .fold(
            || FxHashMap::default(),
            |mut acc: FxHashMap<Arc<String>, Vec<Arc<String>>>, x| {
                acc.entry(x.0.clone()).or_default().push(x.1.clone());
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

    let group_by_protein = group_by_peptides
        .par_iter()
        .fold(
            || FxHashMap::default(),
            |mut acc: FxHashMap<Vec<Arc<String>>, Vec<Arc<String>>>, (key, val)| {
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

    let map_group_peptides: FxHashMap<_, _> = group_by_protein
        .values()
        .cloned()
        .collect::<Vec<_>>()
        .into_par_iter()
        .flat_map(|x1| {
            x1.clone()
                .into_par_iter()
                .map(|x2| (x2, x1.clone()))
                .collect::<Vec<_>>()
        })
        .collect();

    let grouped_peptides: FxHashSet<_> = data
        .into_par_iter()
        .map(|(peptide, protein)| (map_group_peptides[&peptide].clone(), protein))
        .collect();

    grouped_peptides
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
    let group_peptides_proteins = data
        .par_iter()
        .fold(
            || FxHashMap::default(),
            |mut acc: FxHashMap<Vec<Arc<String>>, Vec<Vec<Arc<String>>>>, x| {
                acc.entry(x.0.clone()).or_default().push(x.1.clone());
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

    let (mut g_peps, g_proteins): (Vec<_>, Vec<_>) =
        group_peptides_proteins.clone().into_par_iter().unzip();

    g_peps.iter_mut().for_each(|x| {
        x.sort();
    });
    let g_proteins_sets = g_proteins
        .par_iter()
        .map(|a| {
            let set: HashSet<_> = a.iter().flatten().map(|x| Arc::clone(x)).collect();
            set
        })
        .collect::<Vec<_>>();
    let mut cluster = g_peps
        .par_iter()
        .enumerate()
        .map(|(i, _)| {
            let mut mini_cluster = HashSet::new();
            let mut in_proteins_set = g_proteins_sets[i].clone();
            let mut cluster_len = 1; // Initialize to 1 to enter the loop
            while cluster_len != mini_cluster.len() {
                cluster_len = mini_cluster.len();
                let x = check_protein_group(&g_proteins, &g_proteins_sets, &in_proteins_set)
                    .flat_map(|pg| {
                        mini_cluster.insert(&g_peps[pg.0]);
                        pg.1.iter().flatten().map(|x| Arc::clone(x))
                    })
                    .collect::<HashSet<_>>();
                in_proteins_set.extend(x);
            }
            let mut mini_cluster = mini_cluster.into_iter().collect::<Vec<_>>();
            mini_cluster.sort_unstable();
            mini_cluster
        })
        .collect::<Vec<_>>();

    cluster.par_sort_unstable();
    cluster.dedup();

    let group_peptides_proteins = group_peptides_proteins
        .into_par_iter()
        .map(|(k, v)| {
            let mut k = k;
            k.sort();
            (k, v)
        })
        .collect::<FxHashMap<_, _>>();

    let mut cluster_mapping = data
        .into_par_iter()
        .map(|x1| {
            let mut peps = x1.0.clone();
            peps.sort();
            let e_cluster = cluster
                .iter()
                .enumerate()
                .find(|x| x.1.contains(&&peps))
                .unwrap();
            let proteins = group_peptides_proteins.get(&peps).unwrap().clone();

            let protein_match = proteins
                .into_iter()
                .map(|x1| (peps.clone(), x1))
                .collect::<Vec<_>>();

            (e_cluster.0, protein_match)
        })
        .collect::<Vec<_>>();

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

    cluster_mapping
}

fn check_protein_group<'a>(
    g_proteins: &'a Vec<Vec<Vec<Arc<String>>>>,
    g_proteins_set: &'a Vec<HashSet<Arc<String>>>,
    in_proteins_set: &'a HashSet<Arc<String>>,
) -> impl Iterator<Item = (usize, &'a Vec<Vec<Arc<String>>>)> + 'a {
    g_proteins
        .iter()
        .zip(g_proteins_set)
        .enumerate()
        .filter_map(move |(i, (proteins, proteins_set))| {
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
