use crate::database::IndexedDatabase;
use crate::scoring::Feature;
use itertools::Itertools;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefIterator;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;

pub fn picked_protein(db: &IndexedDatabase, features: &mut [Feature]) {
    let pep_proteins = features
        .iter()
        .filter(|feat| feat.label != -1 && feat.peptide_q < 0.01)
        .flat_map(|feature| {
            let idx = feature.peptide_idx;
            let proteins = &db[feature.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
            let peptide_proteins = proteins
                .split(";")
                .map(|s| (idx.0.to_string(), s.to_string()))
                .collect_vec();

            /* for pep_pro in &peptide_proteins {
                println!("{},{}", pep_pro.0, pep_pro.1);
            }*/

            peptide_proteins
        })
        .collect_vec();

    let protein_map = get_protein_map(pep_proteins);

    // println!("----------------------------------");
    /*for (i, j) in protein_map.iter() {
        println!("{},{}", i, j);
    }*/

    features.par_iter_mut().for_each(|feat| {
        let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
        let any_protein = match proteins
            .split(";")
            .next()
            .map(|s| protein_map.get(s))
            .unwrap()
        {
            None => Some(proteins),
            Some(protein) => Some(protein.to_string()),
        };
        feat.idpicker_proteingroups = any_protein;
    });
}

fn get_protein_map(pep_proteins: Vec<(String, String)>) -> HashMap<String, String> {
    let (peptides, proteins): (Vec<_>, Vec<_>) = pep_proteins.into_iter().unzip();

    let grp_peptides = group_node_with_same_edge(&peptides, &proteins);

    let (grouped_peptides, grouped_proteins): (Vec<_>, Vec<_>) = grp_peptides.into_iter().unzip();
    let grp_proteins = group_node_with_protein(&grouped_peptides, &grouped_proteins);

    let final_cluster: HashSet<_> = separate_into_clusters(grp_proteins).into_iter().collect();

    let final_proteins: Vec<(String, String)> =
        final_cluster.into_iter().flat_map(reduce_cluster).collect();

    let _protein_map: HashMap<_, _> = final_proteins.into_iter().collect();

    _protein_map
}

fn group_node_with_same_edge<'a>(
    l1: &'a Vec<String>,
    l2: &'a Vec<String>,
) -> Vec<(Vec<&'a str>, &'a str)> {
    let mut initial_map: HashMap<&str, Vec<&str>> = HashMap::new();

    // Populate initial_map
    for (key, value) in l1.iter().zip(l2.iter()) {
        initial_map.entry(key).or_default().push(value);
    }

    // Create reverse_map to group by values
    let mut reverse_map: BTreeMap<Vec<&str>, Vec<&str>> = BTreeMap::new();
    for (key, mut values) in initial_map {
        values.sort();
        reverse_map.entry(values).or_default().push(key);
    }

    // Create groups HashMap
    let groups: HashMap<&str, &Vec<&str>> = reverse_map
        .values()
        .flat_map(|value| value.iter().map(move |item| (*item, value)))
        .collect();

    // Build l_set vector with groups
    let mut l_set: Vec<HashSet<&str>> = l1
        .iter()
        .map(|key| {
            let mut set: HashSet<&str> = HashSet::new();
            set.insert(key);
            set.extend(groups[&key[..]]);
            set
        })
        .collect();

    // Create result set
    let result: HashSet<_> = l_set
        .into_iter()
        .zip(l2.iter())
        .map(|(set, psp)| {
            let mut sorted_set: Vec<_> = set.into_iter().collect();
            sorted_set.sort();
            (sorted_set, psp.as_str())
        })
        .collect();

    result.into_iter().collect()
}

fn group_node_with_protein<'a>(
    l2: &'a Vec<Vec<&str>>,
    l1: &'a Vec<&str>,
) -> Vec<(&'a Vec<&'a str>, Vec<&'a str>)> {
    let mut initial_map: HashMap<&str, Vec<&Vec<&str>>> = HashMap::new();

    // Populate initial_map
    for (key, value) in l1.iter().zip(l2.iter()) {
        initial_map.entry(key).or_default().push(value);
    }

    // Reverse map: group proteins by shared peptide sets
    let mut reverse_map: BTreeMap<Vec<&Vec<&str>>, Vec<&str>> = BTreeMap::new();
    for (key, mut values) in initial_map {
        values.sort();
        reverse_map.entry(values).or_default().push(key);
    }

    // Flatten reverse_map into a HashMap for quick lookup
    let groups: HashMap<&str, &Vec<&str>> = reverse_map
        .values()
        .flat_map(|value| value.iter().map(move |item| (*item, value)))
        .collect();

    // Construct l_set: Group proteins based on shared peptides
    let mut l_set: Vec<HashSet<&str>> = l1
        .into_iter()
        .map(|key| {
            let mut set: HashSet<&str> = HashSet::new();
            set.insert(key);
            set.extend(groups[&key[..]]);
            set
        })
        .collect();

    // Sort and create unique result set
    let result: HashSet<_> = l_set
        .into_iter()
        .map(|peps| {
            let mut sorted_set: Vec<_> = peps.into_iter().collect();
            sorted_set.sort();
            sorted_set
        })
        .zip(l2.iter())
        .map(|(x, y)| (y, x))
        .collect();

    let result: Vec<_> = result.into_iter().collect();

    result
}

fn separate_into_clusters<'a>(
    data: Vec<(&'a Vec<&'a str>, Vec<&'a str>)>,
) -> Vec<Vec<(Vec<&'a str>, Vec<&'a str>)>> {
    let mut node_dict: HashMap<&str, HashSet<&str>> = HashMap::new();

    // Build adjacency map of peptides to proteins
    data.iter().for_each(|(peptides, proteins)| {
        peptides.iter().for_each(|peptide| {
            node_dict.entry(peptide).or_default().extend(proteins);
        });
    });

    let mut used_keys = HashSet::new();
    let mut clusters = Vec::new();

    while used_keys.len() < node_dict.len() {
        let root = node_dict
            .keys()
            .find(|k| !used_keys.contains(*k))
            .unwrap()
            .clone();
        let mut proteins = node_dict[&root].clone();
        let mut cluster = HashSet::from([root.clone()]);

        loop {
            let previous_size = proteins.len();

            let keys: Vec<_> = node_dict
                .keys()
                .filter(|key| !cluster.contains(*key) && !node_dict[*key].is_disjoint(&proteins))
                .collect();

            keys.iter().for_each(|key| {
                cluster.insert(key.clone());
                proteins.extend(node_dict[*key].iter().cloned());
            });

            if proteins.len() == previous_size {
                break;
            }
        }

        clusters.push(cluster.clone());
        used_keys.extend(cluster);
    }

    // Group data into clusters
    let mut final_clusters = vec![Vec::new(); clusters.len()];

    data.into_iter().for_each(|(peptides, proteins)| {
        clusters
            .iter()
            .enumerate()
            .filter(|(_, cluster)| peptides.iter().any(|p| cluster.contains(p)))
            .for_each(|(i, _)| final_clusters[i].push((peptides.clone(), proteins.clone())));
    });

    final_clusters
}

fn reduce_cluster<'a>(data: Vec<(Vec<&'a str>, Vec<&'a str>)>) -> Vec<(String, String)> {
    let mut proteins = Vec::new();
    let mut used_peptides = HashSet::new();
    let total_peptides: HashSet<_> = data.iter().map(|(p, _)| p).collect();

    // Compute overall protein scores
    let mut overall_score_dict = HashMap::new();
    data.iter().for_each(|(_, proteins)| {
        *overall_score_dict.entry(proteins.clone()).or_insert(0) += 1;
    });

    while used_peptides.len() < total_peptides.len() {
        let mut prot_dict: HashMap<Vec<&str>, HashSet<Vec<&str>>> = HashMap::new();

        data.iter()
            .filter(|(peptides, _)| !used_peptides.contains(peptides))
            .for_each(|(peptides, proteins)| {
                prot_dict
                    .entry(proteins.clone())
                    .or_default()
                    .insert(peptides.clone());
            });

        // Compute scores for sorting
        let mut score_dict: BTreeMap<Vec<&str>, (usize, usize)> = prot_dict
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
            let protein = top_key[0].clone();
            let protein_str = format!("{}/{}", top_key.len(), top_key.iter().join("/"));
            // println!("{},{}", protein, protein_str);
            // proteins.push((protein, protein_str));
            proteins.extend(
                top_key
                    .iter()
                    .map(|x| (x.to_string(), protein_str.clone()))
                    .collect::<Vec<_>>(),
            );

            if let Some(peptides) = prot_dict.get(top_key).cloned() {
                used_peptides.extend(peptides);
            }
        }
    }

    proteins
}

fn group_node() {
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

    let (peptides, proteins) = test_data.into_iter().unzip();

    let data1 = group_node_with_same_edge(&peptides, &proteins);

    for i in 0..data1.len() {
        println!("PEP ->{:?} ,  PRO-> {:?} ", data1[i].0, data1[i].1);
    }
    let (l1, l2): (Vec<_>, Vec<_>) = data1.into_iter().unzip();

    let data2 = group_node_with_protein(&l1, &l2);
    println!("-------------------------------------------------");
    for index in 0..data2.len() {
        println!("PRO ->{:?} ,  PEP-> {:?} ", data2[index].0, data2[index].1);
    }

    let mut final_proteins = vec![];
    let final_cluster = separate_into_clusters(data2);
    let mut final_cluster: HashSet<_> = final_cluster.into_iter().collect();
    for cluster in final_cluster {
        println!("cluster -> {:?}", cluster);
        let rc = reduce_cluster(cluster);
        println!("rc -> {:?}", rc);
        final_proteins.extend(rc);
    }

    let map: HashMap<_, _> = final_proteins.into_iter().collect();
    for (i, j) in map.iter() {
        println!("Protein  ->{:?} ,  Protein Group-> {:?} ", i, j);
    }
    /*for i in 0..final_proteins.len() {
        println!("Protein  ->{:?} ,  Protein Group-> {:?} ", final_proteins[i].0, final_proteins[i].1);
    }*/
}

#[cfg(test)]
mod test {
    use crate::idpicker::{
        group_node, group_node_with_protein, group_node_with_same_edge, reduce_cluster,
        separate_into_clusters,
    };
    use std::collections::{HashMap, HashSet};

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

        // group_node();
        let (peptides, proteins): (Vec<_>, Vec<_>) = test_data.into_iter().unzip();

        let data1 = group_node_with_same_edge(&peptides, &proteins);

        let (grouped_peptides, grouped_proteins): (Vec<_>, Vec<_>) = data1.into_iter().unzip();
        let data2 = group_node_with_protein(&grouped_peptides, &grouped_proteins);

        let final_cluster: HashSet<_> = separate_into_clusters(data2).into_iter().collect();

        let final_proteins: Vec<(String, String)> =
            final_cluster.into_iter().flat_map(reduce_cluster).collect();

        let _protein_map: HashMap<_, _> = final_proteins.into_iter().collect();

        for (i, j) in _protein_map.iter() {
            println!("Protein  ->{:?} ,  Protein Group-> {:?} ", i, j);
        }

        assert_eq!(5, _protein_map.len());
    }
}
