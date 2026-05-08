//! Protein grouping and inference using IDPicker-style bipartite graph analysis
//!
//! Proteins with identical peptide evidence are collapsed into groups, then a
//! greedy set cover finds an (almost) minimal set of protein groups that
//! explains all observed peptides.
//!
//! Reference: Zhang, B., Chambers, M. C., & Tabb, D. L. (2007). Proteomic
//! parsimony through bipartite graph analysis improves accuracy and transparency.
//! J. Proteome Res., 6(9), 3549-3557. https://doi.org/10.1021/pr070230d

use crate::database::{IndexedDatabase, PeptideIx};
use crate::scoring::Feature;
use fnv::{FnvHashMap, FnvHashSet};
use itertools::Itertools;
use log::info;
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant;

/// Compact protein identifier used during grouping
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct ProteinIx(u32);

impl ProteinIx {
    fn format(
        &self,
        proteins: &[(Arc<str>, bool)],
        decoy_tag: &str,
        generate_decoys: bool,
    ) -> String {
        let (name, decoy) = &proteins[self.0 as usize];
        if *decoy && generate_decoys {
            format!("{}{}", decoy_tag, name)
        } else {
            name.to_string()
        }
    }
}

/// A group of proteins with identical peptide evidence
#[derive(Debug, Default, PartialEq, Eq, Hash, Clone, PartialOrd, Ord)]
struct ProteinGroup(Vec<ProteinIx>);

impl ProteinGroup {
    fn format(
        &self,
        proteins: &[(Arc<str>, bool)],
        decoy_tag: &str,
        generate_decoys: bool,
    ) -> String {
        self.0
            .iter()
            .map(|ix| ix.format(proteins, decoy_tag, generate_decoys))
            .sorted()
            .join("/")
    }
}

/// Bipartite graph for greedy set cover of proteins (left) to peptides (right)
struct BipartiteGraph {
    edges: Vec<(u32, u32)>,
    /// Degree of each left node at construction time (tiebreaker)
    original_degree: Vec<u32>,
    /// Current degree of each left node
    left_degree: Vec<u32>,
    /// Current degree of each right node
    right_degree: Vec<u32>,
    /// Whether each left node is in the cover
    left_cover: Vec<bool>,
    /// Whether each right node is covered
    right_cover: Vec<bool>,
}

impl BipartiteGraph {
    fn new(edges: Vec<(u32, u32)>, left_count: usize, right_count: usize) -> Self {
        let mut left_degree = vec![0u32; left_count];
        let mut right_degree = vec![0u32; right_count];
        for &(l, r) in &edges {
            left_degree[l as usize] += 1;
            right_degree[r as usize] += 1;
        }
        Self {
            edges,
            original_degree: left_degree.clone(),
            left_degree,
            right_degree,
            left_cover: vec![false; left_count],
            right_cover: vec![false; right_count],
        }
    }

    /// Compute a greedy minimal cover of left nodes that explains all right nodes
    fn into_cover(mut self) -> Vec<bool> {
        while !self.edges.is_empty() {
            self.trim();
            if !self.edges.is_empty() {
                self.add_largest_to_cover();
            }
        }
        self.left_cover
    }

    /// Iteratively add left nodes connected to degree-1 right nodes (unique
    /// peptides force their protein into the cover), then remove all edges
    /// incident to covered nodes.
    fn trim(&mut self) {
        let mut prev_len = 0;
        while prev_len != self.edges.len() {
            prev_len = self.edges.len();

            // Any right node with degree 1 forces its left neighbor into the cover
            for &(l, r) in &self.edges {
                if self.right_degree[r as usize] == 1 {
                    self.left_cover[l as usize] = true;
                }
            }

            // Remove edges where the left node is now covered
            self.edges.retain(|&(l, r)| {
                if self.left_cover[l as usize] {
                    self.right_cover[r as usize] = true;
                    self.left_degree[l as usize] -= 1;
                    self.right_degree[r as usize] -= 1;
                    false
                } else {
                    true
                }
            });

            // Remove edges where the right node is already covered
            self.edges.retain(|&(l, r)| {
                if self.right_cover[r as usize] {
                    self.left_degree[l as usize] -= 1;
                    self.right_degree[r as usize] -= 1;
                    false
                } else {
                    true
                }
            });
        }
    }

    /// Add the left node with the most remaining connections to the cover.
    /// Ties broken by original degree (prefer proteins with more total evidence).
    fn add_largest_to_cover(&mut self) {
        if let Some((idx, _)) = self
            .left_degree
            .iter()
            .zip(&self.original_degree)
            .enumerate()
            .max_by_key(|(_, (remaining, original))| (*remaining, *original))
        {
            self.left_cover[idx] = true;
        }
    }
}

/// Assigns proteins to integer indices and groups them by shared peptide evidence
struct ProteinGrouper {
    /// Map from (protein_name, is_decoy) -> ProteinIx
    protein_index: FnvHashMap<(Arc<str>, bool), ProteinIx>,
    /// Protein groups discovered by collapsing identical evidence
    groups: Vec<ProteinGroup>,
    /// Edges from group index -> meta-peptide index
    edges: Vec<(u32, u32)>,
    /// Number of distinct peptides
    peptide_count: usize,
}

impl ProteinGrouper {
    fn build(db: &IndexedDatabase, peptides: FnvHashSet<PeptideIx>) -> Self {
        let mut protein_index: FnvHashMap<(Arc<str>, bool), ProteinIx> = FnvHashMap::default();

        // Map each peptide to a sorted vector of ProteinIx ("meta-peptide"),
        // deduplicating peptides that map to identical protein sets
        let meta_peptides: FnvHashSet<Vec<ProteinIx>> = peptides
            .into_iter()
            .sorted()
            .map(|pep_ix| {
                let peptide = &db[pep_ix];
                peptide
                    .proteins
                    .iter()
                    .map(|name| {
                        let key = (name.clone(), peptide.decoy);
                        let next_id = ProteinIx(protein_index.len() as u32);
                        *protein_index.entry(key).or_insert(next_id)
                    })
                    .sorted()
                    .collect()
            })
            .collect();

        info!("-  found {} meta peptides", meta_peptides.len(),);

        // Invert: group proteins that share identical meta-peptide sets
        let mut prot_to_metapeps: FnvHashMap<ProteinIx, Vec<usize>> = FnvHashMap::default();
        for (i, meta_pep) in meta_peptides.iter().sorted().enumerate() {
            for &prot_ix in meta_pep {
                prot_to_metapeps.entry(prot_ix).or_default().push(i);
            }
        }

        // Proteins with identical meta-peptide vectors form a group
        let mut evidence_to_group: FnvHashMap<Vec<usize>, ProteinGroup> = FnvHashMap::default();
        for (prot_ix, meta_peps) in prot_to_metapeps {
            evidence_to_group
                .entry(meta_peps)
                .or_default()
                .0
                .push(prot_ix);
        }

        let mut groups = Vec::new();
        let mut edges = Vec::new();
        for (group_idx, (meta_peps, group)) in evidence_to_group.into_iter().sorted().enumerate() {
            groups.push(group);
            for meta_pep_idx in meta_peps {
                edges.push((group_idx as u32, meta_pep_idx as u32));
            }
        }

        info!("-  found {} protein groups", groups.len());

        Self {
            protein_index,
            groups,
            edges,
            peptide_count: meta_peptides.len(),
        }
    }

    /// Run set cover inference and produce a lookup map for annotation
    fn into_group_map(self) -> ProteinGroupLookup {
        let group_count = self.groups.len();
        let cover = BipartiteGraph::new(self.edges, group_count, self.peptide_count).into_cover();

        // Build protein name list ordered by ProteinIx
        let proteins: Vec<(Arc<str>, bool)> = self
            .protein_index
            .into_iter()
            .sorted_by_key(|(_, ix)| ix.0)
            .map(|(key, _)| key)
            .collect();

        // Map each protein to its covered group indices
        let mut protein_to_groups: FnvHashMap<(Arc<str>, bool), Vec<u32>> = FnvHashMap::default();
        for (i, in_cover) in cover.into_iter().enumerate() {
            if !in_cover {
                continue;
            }
            for &prot_ix in &self.groups[i].0 {
                let (name, decoy) = &proteins[prot_ix.0 as usize];
                protein_to_groups
                    .entry((name.clone(), *decoy))
                    .or_default()
                    .push(i as u32);
            }
        }

        ProteinGroupLookup {
            groups: self.groups,
            proteins,
            protein_to_groups,
        }
    }
}

/// Maps individual proteins to their group strings for feature annotation
struct ProteinGroupLookup {
    groups: Vec<ProteinGroup>,
    proteins: Vec<(Arc<str>, bool)>,
    protein_to_groups: FnvHashMap<(Arc<str>, bool), Vec<u32>>,
}

impl ProteinGroupLookup {
    /// Get the sorted, semicolon-delimited protein group string for a peptide
    fn group_string(
        &self,
        peptide: &crate::peptide::Peptide,
        db: &IndexedDatabase,
    ) -> Option<String> {
        let group_set: FnvHashSet<&ProteinGroup> = peptide
            .proteins
            .iter()
            .filter_map(|name| self.protein_to_groups.get(&(name.clone(), peptide.decoy)))
            .flat_map(|indices| indices.iter().map(|&i| &self.groups[i as usize]))
            .collect();

        if group_set.is_empty() {
            return None;
        }

        Some(
            group_set
                .into_iter()
                .map(|g| g.format(&self.proteins, &db.decoy_tag, db.generate_decoys))
                .sorted()
                .join(";"),
        )
    }
}

/// Annotate features with protein group information.
///
/// When `protein_grouping` is enabled, proteins are grouped by shared peptide
/// evidence and a minimal cover is inferred. If `confident_peptide_threshold`
/// is provided, an initial grouping pass is run on only confident peptides,
/// followed by a second pass on all peptides.
///
/// Features not assigned to a group are annotated with their raw protein list.
pub fn generate_protein_groups(
    db: &IndexedDatabase,
    features: &mut [Feature],
    protein_grouping: bool,
    confident_peptide_threshold: Option<f32>,
) {
    let time = Instant::now();
    if protein_grouping {
        if confident_peptide_threshold.is_some() {
            annotate_features(features, db, confident_peptide_threshold);
        }
        annotate_features(features, db, None);
    }

    // Fallback: features not assigned by grouping get their raw protein list
    features
        .par_iter_mut()
        .filter(|f| f.protein_groups.is_none())
        .for_each(|feat| {
            let pep = &db[feat.peptide_idx];
            feat.protein_groups = Some(pep.proteins(&db.decoy_tag, db.generate_decoys));
            feat.num_protein_groups = pep.proteins.len() as u32;
        });
    info!(
        "Grouped and inferred proteins in {:?}ms",
        time.elapsed().as_millis()
    );
}

fn annotate_features(
    features: &mut [Feature],
    db: &IndexedDatabase,
    confident_peptide_threshold: Option<f32>,
) {
    let time = Instant::now();
    let threshold = confident_peptide_threshold.unwrap_or(1.0).clamp(0.0, 1.0);

    let peptides: FnvHashSet<PeptideIx> = features
        .par_iter()
        .filter(|f| f.label != -1 && f.peptide_q < threshold)
        .map(|f| f.peptide_idx)
        .collect();

    info!(
        "Protein grouping: {} unique peptides (threshold={}) in {:?}ms",
        peptides.len(),
        threshold,
        time.elapsed().as_millis()
    );

    let grouper = ProteinGrouper::build(db, peptides);
    let lookup = grouper.into_group_map();

    let annotated: u32 = features
        .par_iter_mut()
        .filter(|f| f.protein_groups.is_none())
        .map(|feat| {
            let pep = &db[feat.peptide_idx];
            match lookup.group_string(pep, db) {
                Some(groups) => {
                    feat.num_protein_groups = groups.matches(';').count() as u32 + 1;
                    feat.protein_groups = Some(groups);
                    1u32
                }
                None => 0,
            }
        })
        .sum();

    info!(
        "-  annotated {} features in {:?}ms",
        annotated,
        time.elapsed().as_millis()
    );
}

#[cfg(test)]
mod test {
    use super::*;
    use std::sync::Arc;

    use crate::peptide::Peptide;

    fn get_data() -> (Vec<Vec<&'static str>>, Vec<bool>, Vec<f32>) {
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
        proteins: &[Vec<&str>],
        decoys: &[bool],
        q_vals: &[f32],
    ) -> (IndexedDatabase, Vec<Feature>) {
        build_db_and_features_with_scores(proteins, decoys, q_vals, None)
    }

    fn build_db_and_features_with_scores(
        proteins: &[Vec<&str>],
        decoys: &[bool],
        q_vals: &[f32],
        scores: Option<&[f32]>,
    ) -> (IndexedDatabase, Vec<Feature>) {
        let features: Vec<Feature> = (0..proteins.len())
            .map(|ix| Feature {
                peptide_idx: PeptideIx(ix as u32),
                label: if decoys[ix] { -1 } else { 1 },
                peptide_q: q_vals[ix],
                discriminant_score: scores.map(|s| s[ix]).unwrap_or(0.0),
                ..Default::default()
            })
            .collect();
        let db = IndexedDatabase {
            peptides: proteins
                .iter()
                .zip(decoys)
                .map(|(prots, &decoy)| Peptide {
                    proteins: prots.iter().map(|&s| Arc::from(s)).collect(),
                    decoy,
                    ..Default::default()
                })
                .collect(),
            decoy_tag: "rev_".to_string(),
            generate_decoys: false,
            ..IndexedDatabase::default()
        };
        (db, features)
    }

    #[test]
    fn test_protein_grouping_expected_groups() {
        let (proteins, decoys, q_vals) = get_data();
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));
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
        ];
        let actual: Vec<_> = features
            .iter()
            .map(|v| v.protein_groups.as_ref().unwrap().as_str())
            .collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_bipartite_cover_unique_peptides() {
        // Three proteins each with a single unique peptide — all in cover
        let edges = vec![(0, 0), (1, 1), (2, 2)];
        let cover = BipartiteGraph::new(edges, 3, 3).into_cover();
        assert_eq!(cover, vec![true, true, true]);
    }

    #[test]
    fn test_bipartite_cover_subset_protein() {
        // protein 0 -> peptides {0, 1, 2}  (superset)
        // protein 1 -> peptides {0, 1}     (subset)
        let edges = vec![(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)];
        let cover = BipartiteGraph::new(edges, 2, 3).into_cover();
        assert!(cover[0], "superset protein should be covered");
        assert!(!cover[1], "subset protein should not be covered");
    }

    #[test]
    fn test_bipartite_cover_shared_peptide() {
        // protein 0 -> peptides {0, 1}
        // protein 1 -> peptides {1, 2}
        // Both should be in cover because each has a unique peptide
        let edges = vec![(0, 0), (0, 1), (1, 1), (1, 2)];
        let cover = BipartiteGraph::new(edges, 2, 3).into_cover();
        assert_eq!(cover, vec![true, true]);
    }

    #[test]
    fn test_bipartite_cover_empty() {
        let cover = BipartiteGraph::new(vec![], 0, 0).into_cover();
        assert!(cover.is_empty());
    }

    #[test]
    fn test_bipartite_cover_single() {
        let cover = BipartiteGraph::new(vec![(0, 0)], 1, 1).into_cover();
        assert_eq!(cover, vec![true]);
    }

    #[test]
    fn test_decoy_features_excluded_from_grouping() {
        let proteins = vec![vec!["protA"], vec!["protA"], vec!["protB"]];
        let decoys = vec![false, true, false];
        let q_vals = vec![0.0, 0.0, 0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));

        for feat in &features {
            assert!(
                feat.protein_groups.is_some(),
                "every feature should be annotated"
            );
        }
        assert_eq!(features[1].protein_groups.as_deref(), Some("protA"));
    }

    #[test]
    fn test_decoy_features_with_generate_decoys() {
        let proteins = vec![vec!["protA"], vec!["protA"]];
        let decoys = vec![false, true];
        let q_vals = vec![0.0, 0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);

        // Override generate_decoys for this test
        let db = IndexedDatabase {
            generate_decoys: true,
            ..db
        };

        generate_protein_groups(&db, &mut features, false, None);
        assert_eq!(features[0].protein_groups.as_deref(), Some("protA"));
        assert_eq!(features[1].protein_groups.as_deref(), Some("rev_protA"));
    }

    #[test]
    fn test_grouping_disabled_falls_back_to_annotate() {
        let proteins = vec![vec!["protA", "protB"], vec!["protC"]];
        let decoys = vec![false, false];
        let q_vals = vec![0.0, 0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, false, None);

        assert_eq!(features[0].protein_groups.as_deref(), Some("protA;protB"));
        assert_eq!(features[0].num_protein_groups, 2);
        assert_eq!(features[1].protein_groups.as_deref(), Some("protC"));
        assert_eq!(features[1].num_protein_groups, 1);
    }

    #[test]
    fn test_single_protein_single_peptide() {
        let proteins = vec![vec!["protA"]];
        let decoys = vec![false];
        let q_vals = vec![0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));

        assert_eq!(features[0].protein_groups.as_deref(), Some("protA"));
        assert_eq!(features[0].num_protein_groups, 1);
    }

    #[test]
    fn test_all_shared_peptides() {
        // Every peptide maps to the same two proteins — they should form a single group
        let proteins = vec![
            vec!["protA", "protB"],
            vec!["protA", "protB"],
            vec!["protA", "protB"],
        ];
        let decoys = vec![false, false, false];
        let q_vals = vec![0.0, 0.0, 0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));

        let group = features[0].protein_groups.as_deref().unwrap();
        assert!(group.contains("protA"));
        assert!(group.contains("protB"));
        for feat in &features {
            assert_eq!(feat.protein_groups.as_deref(), Some(group));
            assert_eq!(feat.num_protein_groups, 1);
        }
    }

    #[test]
    fn test_peptide_fdr_threshold_filtering() {
        let proteins = vec![vec!["protA"], vec!["protB"]];
        let decoys = vec![false, false];
        let q_vals = vec![0.001, 0.5];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));

        assert!(features[0].protein_groups.is_some());
        assert!(features[1].protein_groups.is_some());
    }

    #[test]
    fn test_all_decoy_features() {
        let proteins = vec![vec!["protA"], vec!["protB"]];
        let decoys = vec![true, true];
        let q_vals = vec![0.0, 0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));

        for feat in &features {
            assert!(feat.protein_groups.is_some());
        }
    }

    #[test]
    fn test_proteins_with_identical_evidence_are_grouped() {
        let proteins = vec![
            vec!["protA", "protB"],
            vec!["protA", "protB"],
            vec!["protC"],
        ];
        let decoys = vec![false, false, false];
        let q_vals = vec![0.0, 0.0, 0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));

        let group_01 = features[0].protein_groups.as_deref().unwrap();
        assert!(group_01.contains("protA") && group_01.contains("protB"));
        assert_eq!(features[0].protein_groups, features[1].protein_groups);
        assert_eq!(features[2].protein_groups.as_deref(), Some("protC"));
    }

    #[test]
    fn test_num_protein_groups_counts_distinct_groups() {
        // protA and protB each have unique evidence; peptide 2 is shared across groups
        let proteins = vec![vec!["protA"], vec!["protB"], vec!["protA", "protB"]];
        let decoys = vec![false, false, false];
        let q_vals = vec![0.0, 0.0, 0.0];
        let (db, mut features) = build_db_and_features(&proteins, &decoys, &q_vals);
        generate_protein_groups(&db, &mut features, true, Some(0.01));

        assert_eq!(features[0].num_protein_groups, 1);
        assert_eq!(features[1].num_protein_groups, 1);
        assert_eq!(features[2].num_protein_groups, 2);
    }
}
