//! False discovery rate control using double-competition (picked-peptide &
//! picked-protein) approaches
//!
//! Lin et al., https://pubmed.ncbi.nlm.nih.gov/36166314/
//! Savitski et al., https://pubmed.ncbi.nlm.nih.gov/25987413/

use rayon::slice::ParallelSliceMut;

use crate::{
    database::{IndexedDatabase, PeptideIx},
    mass::Residue,
    peptide::Peptide,
    scoring::Percolator,
};
use std::{collections::HashMap, hash::Hash};

#[derive(Copy, Clone, Debug)]
struct Competition<Ix: Default> {
    forward: f32,
    foward_ix: Ix,
    reverse: f32,
    reverse_ix: Ix,
    q_value: f32,
}

impl<Ix: Default> Default for Competition<Ix> {
    fn default() -> Self {
        Self {
            forward: f32::MIN,
            reverse: f32::MIN,
            foward_ix: Ix::default(),
            reverse_ix: Ix::default(),
            q_value: 1.0,
        }
    }
}

impl<Ix: Default> Competition<Ix> {
    fn score(&self) -> f32 {
        self.forward.max(self.reverse)
    }

    fn is_decoy(&self) -> bool {
        self.reverse >= self.forward
    }

    fn assign_q_value(scores: &mut [Competition<Ix>]) -> usize {
        scores.sort_unstable_by(|a, b| b.score().total_cmp(&a.score()));

        let mut decoy = 1;
        let mut target = 0;
        for score in scores.iter_mut() {
            match score.is_decoy() {
                true => decoy += 1,
                false => target += 1,
            }
            score.q_value = decoy as f32 / target as f32
        }
        let mut q_min = 1.0f32;
        let mut passing = 0;
        for score in scores.iter_mut().rev() {
            q_min = q_min.min(score.q_value);
            score.q_value = q_min;
            if q_min <= 0.01 {
                passing += 1;
            }
        }
        passing
    }
}

pub fn picked_peptide(db: &IndexedDatabase, features: &mut [Percolator]) -> usize {
    let mut map: HashMap<String, Competition<PeptideIx>> = HashMap::new();
    for feat in features.iter() {
        let peptide = &db[feat.peptide_idx];
        let fwd = peptide.pseudo_forward();
        let key = fwd.as_ref().unwrap_or(peptide).to_string();

        let entry = map.entry(key).or_default();
        match peptide.decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = feat.peptide_idx;
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = feat.peptide_idx;
            }
        }
    }
    let mut scores = map.into_values().collect::<Vec<_>>();
    let passing = Competition::assign_q_value(&mut scores);

    let scores = scores
        .into_iter()
        .flat_map(|score| {
            [
                (score.foward_ix, score.q_value),
                (score.reverse_ix, score.q_value),
            ]
        })
        .collect::<HashMap<_, _>>();
    for feat in features.iter_mut() {
        feat.peptide_q = scores[&feat.peptide_idx];
    }

    passing
}

pub fn picked_protein(db: &IndexedDatabase, features: &mut [Percolator]) -> usize {
    let mut map: HashMap<&str, Competition<String>> = HashMap::new();
    for feat in features.iter() {
        let peptide = &db[feat.peptide_idx];
        let entry = map.entry(&feat.proteins).or_default();
        match peptide.decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = feat.proteins.clone();
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = feat.proteins.clone();
            }
        }
    }
    let mut scores = map.into_values().collect::<Vec<_>>();
    let passing = Competition::assign_q_value(&mut scores);

    let scores = scores
        .into_iter()
        .flat_map(|score| {
            [
                (score.foward_ix, score.q_value),
                (score.reverse_ix, score.q_value),
            ]
        })
        .collect::<HashMap<_, _>>();
    for feat in features.iter_mut() {
        feat.peptide_q = scores[&feat.proteins];
    }

    passing
}
