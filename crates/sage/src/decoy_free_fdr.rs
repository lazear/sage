//! Decoy-free FDR estimation using a Gumbel distribution model.
//!
//! This module implements a method for calculating q-values without a decoy database,
//! primarily for increasing sensitivity in low-input proteomics experiments. The statistical
//! model is built by assuming that low-scoring, high-rank PSMs can be used to
//! estimate the null distribution of all incorrect matches.
//!
//! ## References
//! Modeling Lower-Order Statistics to Enable Decoy-Free FDR Estimation in Proteomics
//! Dominik Madej and Henry Lam
//! Journal of Proteome Research 2023 22 (4), 1159-1171
//! https://doi.org/10.1021/acs.jproteome.2c00604
//! https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00604
//! https://github.com/dommad/pylord
//! 
//! New mixture models for decoy-free false discovery rate estimation in mass spectrometry proteomics
//! Yisu Peng, Shantanu Jain, Yong Fuga Li, Michal Greguš, Alexander R. Ivanov, Olga Vitek, Predrag Radivojac
//! Bioinformatics, Volume 36, Issue Supplement_2, December 2020, Pages i745–i753
//! https://doi.org/10.1093/bioinformatics/btaa807
//! https://academic.oup.com/bioinformatics/article/36/Supplement_2/i745/6055912
//! https://github.com/shawn-peng/DecoyFree-MSFDR

//! A Decoy-Free Approach to the Identification of Peptides
//! Giulia Gonnelli, Michiel Stock, Jan Verwaeren, Davy Maddelein, Bernard De Baets, Lennart Martens, and Sven Degroeve
//! Journal of Proteome Research 2015 14 (4), 1792-1798
//! https://doi.org/10.1021/pr501164r
//! https://pubs.acs.org/doi/10.1021/pr501164r
//! https://bio.tools/nokoi

use crate::database::IndexedDatabase;
use crate::scoring::Feature;
use fnv::FnvHashMap;
use rayon::prelude::*;
use statrs::consts::EULER_MASCHERONI;
use statrs::distribution::{ContinuousCDF, Gumbel};
use crate::lfq::{Peak, PrecursorId};

/// Calculate spectrum-level q-values using a Gumbel-based decoy-free method.
/// This function now takes an immutable slice and returns a new, sorted Vec<Feature>.
pub fn calculate_q_values(psms: &[Feature]) -> Vec<Feature> {
    // Create a mutable copy to work with, guaranteeing the result is fresh
    let mut new_features = psms.to_vec();

    log::info!("Building null distribution from high-rank PSMs for decoy-free FDR...");
    let null_scores: Vec<f64> = new_features
        .par_iter()
        .filter_map(|psm| if psm.rank >= 4 { Some(psm.hyperscore as f64) } else { None })
        .collect();

    if null_scores.len() < 50 {
        log::error!("Not enough high-rank PSMs ({}) to build a stable null distribution. All q-values will be set to 1.0.", null_scores.len());
        new_features.par_iter_mut().for_each(|psm| psm.spectrum_q = 1.0);
        return new_features;
    }
    
    let (mu, beta) = fit_gumbel_moments(&null_scores);
    log::info!("Fitted Gumbel distribution: µ = {:.3}, β = {:.3}", mu, beta);

    let gumbel = Gumbel::new(mu, beta).unwrap();
    new_features.iter_mut().for_each(|psm| {
        if psm.rank == 1 {
            psm.spectrum_q = gumbel.sf(psm.hyperscore as f64) as f32;
        } else {
            psm.spectrum_q = 1.0;
        }
    });

    let mut rank_1_psms: Vec<&mut Feature> = new_features
        .iter_mut()
        .filter(|psm| psm.rank == 1)
        .collect();
        
    rank_1_psms.sort_unstable_by(|a, b| a.spectrum_q.total_cmp(&b.spectrum_q));
    
    let n = rank_1_psms.len() as f32;
    let mut last_q: f32 = 1.0;

    for (i, psm) in rank_1_psms.into_iter().rev().enumerate() {
        let rank = (n - i as f32) as f32;
        let q_val = (psm.spectrum_q * n) / rank;
        
        last_q = last_q.min(q_val);
        psm.spectrum_q = last_q;
    }
	
    // Finally, sort the entire results vector by q-value and return it.
    new_features.sort_unstable_by(|a, b| a.spectrum_q.total_cmp(&b.spectrum_q));
    new_features
}

/// A simplified Gumbel parameter estimator using the method of moments.
fn fit_gumbel_moments(scores: &[f64]) -> (f64, f64) {
    let n = scores.len() as f64;
    let mean = scores.iter().sum::<f64>() / n;
    let variance = scores.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / n;
    
    let beta = (variance * 6.0 / std::f64::consts::PI.powi(2)).sqrt();
    let mu = mean - EULER_MASCHERONI * beta;
    (mu, beta)
}

/// Calculate peptide-level q-values in decoy-free mode.
/// A peptide's q-value is the best q-value of any of its constituent PSMs.
pub fn calculate_peptide_q(features: &mut [Feature], db: &IndexedDatabase) -> usize {
    let mut best_q: FnvHashMap<String, f32> = FnvHashMap::default();

    // Find the best spectrum_q for each peptide sequence
    for feat in features.iter().filter(|f| f.rank == 1) {
        let peptide = db[feat.peptide_idx].to_string();
        best_q.entry(peptide)
              .and_modify(|q| *q = q.min(feat.spectrum_q))
              .or_insert(feat.spectrum_q);
    }

    // Assign the best q-value to all PSMs for that peptide
    for feat in features.iter_mut() {
        let peptide = db[feat.peptide_idx].to_string();
        if let Some(q) = best_q.get(&peptide) {
            feat.peptide_q = *q;
        }
    }
    
    // Return number of peptides passing 1% FDR
    best_q.values().filter(|&&q| q <= 0.01).count()
}

/// Calculate protein-level q-values in decoy-free mode.
/// A protein's q-value is the best q-value of any of its constituent peptides.
pub fn calculate_protein_q(features: &mut [Feature], db: &IndexedDatabase) -> usize {
    let mut best_q: FnvHashMap<String, f32> = FnvHashMap::default();

    // Find the best peptide_q for each protein
    for feat in features.iter() {
        // Use the protein name string as the key
        let protein_key = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
        best_q.entry(protein_key)
              .and_modify(|q| *q = q.min(feat.peptide_q))
              .or_insert(feat.peptide_q);
    }
    
    // Assign the best q-value to all PSMs for that protein
    for feat in features.iter_mut() {
        let protein_key = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
        if let Some(q) = best_q.get(&protein_key) {
            feat.protein_q = *q;
        }
    }
    
    // Return number of proteins passing 1% FDR
    best_q.values().filter(|&&q| q <= 0.01).count()
}

	// The new function will have a similar signature to the original `picked_precursor`
	pub fn decoy_free_precursor(
		peaks: &mut FnvHashMap<(PrecursorId, bool), (Peak, Vec<f64>)>,
	) -> usize {
		// 1. Collect scores from the artificial decoy precursors to build the null model.
		let decoy_scores: Vec<f64> = peaks
			.iter()
			.filter_map(|((_id, is_decoy), (peak, _))| {
				if *is_decoy {
					Some(peak.score)
				} else {
					None
				}
			})
			.collect();

		if decoy_scores.len() < 50 {
			log::error!("Not enough decoy precursors ({}) to build stable LFQ null distribution.", decoy_scores.len());
			// Set all q-values to 1.0
			peaks.par_iter_mut().for_each(|(_, (peak, _))| peak.q_value = 1.0);
			return 0;
		}

		// 2. Fit the Gumbel distribution.
		let (mu, beta) = fit_gumbel_moments(&decoy_scores);
		let gumbel = Gumbel::new(mu, beta).unwrap();

		// 3. Calculate p-values for all target precursors.
		//    Collect mutable references to the target peaks to update them later.
		let mut targets: Vec<&mut Peak> = peaks
			.iter_mut()
			.filter_map(|((_id, is_decoy), (peak, _))| {
				if !*is_decoy {
					// Calculate p-value and store it temporarily in the q_value field
					peak.q_value = gumbel.sf(peak.score) as f32;
					Some(peak)
				} else {
					None
				}
			})
			.collect();
			
		// 4. Apply Benjamini-Hochberg correction to the p-values.
		targets.sort_unstable_by(|a, b| a.q_value.total_cmp(&b.q_value));
		
		let n = targets.len() as f32;
		let mut last_q = 1.0;

		for (i, peak) in targets.into_iter().rev().enumerate() {
			let rank = (n - i as f32) as f32;
			let q_val = (peak.q_value * n) / rank;
			
			if q_val < last_q {
				last_q = q_val;
			}
			peak.q_value = last_q;
		}
		
		// 5. Return the count of significant precursors at 5% FDR.
		peaks
			.values()
			.filter(|(peak, _)| !peak.score.is_nan() && peak.q_value <= 0.05)
			.count()
	}