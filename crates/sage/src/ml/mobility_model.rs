//! Retention time prediction using linear regression
//!
//! See Klammer et al., Anal. Chem. 2007, 79, 16, 6111–6118
//! https://doi.org/10.1021/ac070262k

use super::{gauss::Gauss, matrix::Matrix};
use crate::database::IndexedDatabase;
use crate::mass::VALID_AA;
use crate::peptide::Peptide;
use crate::scoring::Feature;
use rayon::prelude::*;

/// Try to fit a retention time prediction model
// Add `decoy_free: bool` to the function signature
pub fn predict(db: &IndexedDatabase, features: &mut [Feature], decoy_free: bool) -> Option<()> {
    // Training LR might fail - not enough values, or r-squared is < 0.7
    // Pass the flag down to the `fit` function
    let lr = match MobilityModel::fit(db, features, decoy_free) {
        Some(lr) => lr,
        None => {
            log::warn!("Mobility model failed to train");
            return None;
        }
    };
    features.par_iter_mut().for_each(|feat| {
        // LR can sometimes predict crazy values - clamp predicted RT
        let ims = lr.predict_peptide(db, feat);
        let bounded = ims.clamp(0.0, 2.0) as f32;
        feat.predicted_ims = bounded;

        feat.delta_ims_model = (feat.ims - bounded).abs();
    });
    Some(())
}
pub struct MobilityModel {
    beta: Vec<f64>,
    map: [usize; 26],
    pub r2: f64,
}

const BULKY_AA_IDXS: [usize; 6] = [
    b'L' as usize - b'A' as usize,
    b'V' as usize - b'A' as usize,
    b'I' as usize - b'A' as usize,
    b'F' as usize - b'A' as usize,
    b'W' as usize - b'A' as usize,
    b'Y' as usize - b'A' as usize,
];

const UNCHARGED_POLAR_AA_IDXS: [usize; 4] = [
    b'S' as usize - b'A' as usize,
    b'T' as usize - b'A' as usize,
    b'N' as usize - b'A' as usize,
    b'Q' as usize - b'A' as usize,
];

const POSITIVE_AA_IDXS: [usize; 3] = [
    b'R' as usize - b'A' as usize,
    b'K' as usize - b'A' as usize,
    b'H' as usize - b'A' as usize,
];

const NEGATIVE_AA_IDXS: [usize; 2] = [b'D' as usize - b'A' as usize, b'E' as usize - b'A' as usize];

const TINY_AA_IDXS: [usize; 3] = [
    b'G' as usize - b'A' as usize,
    b'A' as usize - b'A' as usize,
    b'S' as usize - b'A' as usize,
];

const BRANCHED_AA_IDXS: [usize; 3] = [
    b'L' as usize - b'A' as usize,
    b'I' as usize - b'A' as usize,
    b'V' as usize - b'A' as usize,
];

const FEATURES: usize = VALID_AA.len() * 4 + 12;
const PCT_FEATURES_START: usize = VALID_AA.len();
const N_TERMINAL: usize = VALID_AA.len() * 2;
const C_TERMINAL: usize = VALID_AA.len() * 3;
const NUM_BRANCHED: usize = FEATURES - 12;
const NUM_TINY: usize = FEATURES - 11;
const NUM_UC_POLAR: usize = FEATURES - 10;
const NUM_BULKY: usize = FEATURES - 9;
const NUM_POSITIVE: usize = FEATURES - 8;
const NUM_NEGATIVE: usize = FEATURES - 7;
const INV_PEPTIDE_CHARGE: usize = FEATURES - 6;
const PEPTIDE_CHARGE: usize = FEATURES - 5;
const PEPTIDE_MZ: usize = FEATURES - 4;
const PEPTIDE_LEN: usize = FEATURES - 3;
const PEPTIDE_MASS: usize = FEATURES - 2;
const INTERCEPT: usize = FEATURES - 1;

// IN THEORY we could have only one model for both RT and IM
// And the RT should ignore the charge state "at training time"
impl MobilityModel {
    /// One-hot encoding of peptide sequences into feature vector
    /// Note that this currently does not take into account any modifications
    fn embed(peptide: &Peptide, charge: &u8, map: &[usize; 26]) -> [f64; FEATURES] {
        let mut embedding = [0.0; FEATURES];
        let cterm = peptide.sequence.len().saturating_sub(3);
        let pep_length = peptide.sequence.len() as f64;

        //let default_first_val = 1.0f64;
        for (aa_idx, residue) in peptide.sequence.iter().enumerate() {
            let idx = map[(residue - b'A') as usize];
            embedding[idx] += 1.0;
            // Embed N- and C-terminal AA's
            // 2 on each end
            match aa_idx {
                0 | 1 => embedding[N_TERMINAL + idx] += 1.0,
                x if x > cterm => embedding[C_TERMINAL + idx] += 1.0,
                _ => {}
            }
            let x = idx;

            if BULKY_AA_IDXS.contains(&x) {
                embedding[NUM_BULKY] += 1.0;
            };
            if UNCHARGED_POLAR_AA_IDXS.contains(&x) {
                embedding[NUM_UC_POLAR] += 1.0;
            };
            if POSITIVE_AA_IDXS.contains(&x) {
                embedding[NUM_POSITIVE] += 1.0
            };
            if NEGATIVE_AA_IDXS.contains(&x) {
                embedding[NUM_NEGATIVE] += 1.0
            };
            if TINY_AA_IDXS.contains(&x) {
                embedding[NUM_TINY] += 1.0
            };
            if BRANCHED_AA_IDXS.contains(&x) {
                embedding[NUM_BRANCHED] += 1.0
            };
        }

        // PCT features are just the AA counts divided by the length of the peptide
        for idx in 0..VALID_AA.len() {
            let pct_val = embedding[idx] / pep_length;
            embedding[PCT_FEATURES_START + idx] = pct_val;
        }

        let charge_feature: f64 = *charge as f64;
        embedding[PEPTIDE_CHARGE] = charge_feature;
        embedding[INV_PEPTIDE_CHARGE] = 1. / charge_feature;
        embedding[PEPTIDE_LEN] = peptide.sequence.len() as f64;
        embedding[PEPTIDE_MASS] = (peptide.monoisotopic as f64) / 1000.0;
        embedding[PEPTIDE_MZ] = ((peptide.monoisotopic as f64) / charge_feature) / 1000.0;
        embedding[INTERCEPT] = 1.0;
        embedding
    }

    /// Attempt to fit a linear regression model: peptide sequence + charge ~ retention time
    // Add `decoy_free: bool` to the function signature
    pub fn fit(db: &IndexedDatabase, training_set: &[Feature], decoy_free: bool) -> Option<Self> {

        // Create a mapping from amino acid character to vector embedding
        let mut map = [0; 26];
        for (idx, aa) in VALID_AA.iter().enumerate() {
            map[(aa - b'A') as usize] = idx;
        }

        // Create a filtered iterator of high-quality PSMs ONCE.
        let training_candidates = training_set
            .par_iter()
            .filter(|feat| {
                // This is the new conditional filter logic
                let is_target = if decoy_free {
                    feat.rank == 1
                } else {
                    feat.label == 1
                };
                // In both cases, we only want to train on high-confidence peptides
                is_target && feat.spectrum_q <= 0.01
            });

        // Use the filtered iterator to collect ion mobilities
        let ims = training_candidates.clone().map(|psm| psm.ims as f64).collect::<Vec<f64>>();
		
		// ADD THIS CHECK
		if ims.len() < 10 {
			log::warn!(
				"Not enough high-quality PSMs ({}) to train the ion mobility model.",
				ims.len()
			);
			return None;
		}
		// END ADDITION

        let ims_mean = ims.iter().sum::<f64>() / ims.len() as f64;
        let ims_var = ims.iter().map(|rt| (rt - ims_mean).powi(2)).sum::<f64>();

        let rt = Matrix::col_vector(ims);

        let features = training_set
            .par_iter()
            .filter(|feat| feat.rank == 1 && feat.spectrum_q <= 0.01)
            .flat_map_iter(|psm| Self::embed(&db[psm.peptide_idx], &psm.charge, &map))
            .collect::<Vec<_>>();

        let rows = features.len() / FEATURES;
        let features = Matrix::new(features, rows, FEATURES);

        let f_t = features.transpose();
        let cov = f_t.dot(&features);
        let b = f_t.dot(&rt);

        let beta = Gauss::solve(cov, b)?;

        let predicted_im = features.dot(&beta).take();
        let sum_squared_error = predicted_im
            .iter()
            .zip(rt.take())
            .map(|(pred, act)| (pred - act).powi(2))
            .sum::<f64>();

        let mse: f64 = sum_squared_error / predicted_im.len() as f64;
        let r2 = 1.0 - (sum_squared_error / ims_var);
        log::info!("- fit mobility model, rsq = {}, mse = {}", r2, mse);
        Some(Self {
            beta: beta.take(),
            map,
            r2,
        })
    }

    /// Predict retention times for a collection of PSMs
    pub fn predict_peptide(&self, db: &IndexedDatabase, psm: &Feature) -> f64 {
        let v = Self::embed(&db[psm.peptide_idx], &psm.charge, &self.map);
        v.into_iter()
            .zip(&self.beta)
            .fold(0.0f64, |sum, (x, y)| sum + x * y)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::enzyme::Digest;

    #[test]
    fn test_feature_embed() {
        let peps = vec![
            Peptide::try_from(Digest {
                decoy: false,
                sequence: "LEKSLIEK".into(),
                missed_cleavages: 0,
                ..Default::default()
            })
            .unwrap(),
            Peptide::try_from(Digest {
                decoy: false,
                sequence: "LERSLIEWK".into(),
                missed_cleavages: 0,
                ..Default::default()
            })
            .unwrap(),
            Peptide::try_from(Digest {
                decoy: false,
                sequence: "LWESLIEK".into(),
                missed_cleavages: 0,
                ..Default::default()
            })
            .unwrap(),
            Peptide::try_from(Digest {
                decoy: false,
                sequence: "CHADWICK".into(),
                missed_cleavages: 0,
                ..Default::default()
            })
            .unwrap(),
        ];

        let charge = 2;
        let mut map = [0; 26];
        for (idx, aa) in VALID_AA.iter().enumerate() {
            map[(aa - b'A') as usize] = idx;
        }
        let embeddings: Vec<[f64; FEATURES]> = peps
            .iter()
            .map(|x| MobilityModel::embed(x, &charge, &map))
            .collect();

        let k_idx = map[(b'K' - b'A') as usize];
        let w_idx = map[(b'W' - b'A') as usize];
        let l_idx = map[(b'L' - b'A') as usize];
        let i_idx = map[(b'I' - b'A') as usize];

        let l_nterm_counts = embeddings
            .iter()
            .map(|x| x[N_TERMINAL + l_idx])
            .collect::<Vec<f64>>();
        let k_nterm_counts = embeddings
            .iter()
            .map(|x| x[N_TERMINAL + k_idx])
            .collect::<Vec<f64>>();
        let k_cterm_counts = embeddings
            .iter()
            .map(|x| x[C_TERMINAL + k_idx])
            .collect::<Vec<f64>>();
        let w_nterm_counts = embeddings
            .iter()
            .map(|x| x[N_TERMINAL + w_idx])
            .collect::<Vec<f64>>();
        let w_cterm_counts = embeddings
            .iter()
            .map(|x| x[C_TERMINAL + w_idx])
            .collect::<Vec<f64>>();
        let i_cterm_counts = embeddings
            .iter()
            .map(|x| x[C_TERMINAL + i_idx])
            .collect::<Vec<f64>>();
        assert_eq!(l_nterm_counts, vec![1.0, 1.0, 1.0, 0.0,], "L N-term counts");
        assert_eq!(k_nterm_counts, vec![0.0, 0.0, 0.0, 0.0,], "K N-term counts");
        assert_eq!(w_nterm_counts, vec![0.0, 0.0, 1.0, 0.0,], "W N-term counts");
        assert_eq!(k_cterm_counts, vec![1.0, 1.0, 1.0, 1.0,], "K C-term counts");
        assert_eq!(w_cterm_counts, vec![0.0, 1.0, 0.0, 0.0,], "W C-term counts");
        assert_eq!(i_cterm_counts, vec![0.0, 0.0, 0.0, 0.0,], "I C-term counts");
    }
}
