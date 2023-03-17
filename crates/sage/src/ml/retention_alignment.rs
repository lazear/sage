//! Perform global retention time alignment using a modified algorithm based on
//! Chen AT, et al. DART-ID https://pubmed.ncbi.nlm.nih.gov/31260443/ (Fig 1A/B)
//!
//! 0) Transform all RTs into unit-less percentages (0.0 - 1.0)
//! 1) Assume that the expected RT for a peptide can be estimated from the average
//!    RT across all runs
//! 2) For each run, calculate a linear regression between the observed peptide RTs
//!    and the global average. Transform all PSM retention times by the regression
//!    parameters
//!
//! If LFQ is enabled, MS1 apex times will be used instead of PSM retention times

use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::sync::atomic::AtomicU32;

use super::matrix::Matrix;
use crate::database::PeptideIx;
use crate::scoring::Feature;
use dashmap::DashMap;
use fnv::FnvHasher;
use rayon::prelude::*;

type FnvDashMap<K, V> = DashMap<K, V, BuildHasherDefault<FnvHasher>>;

fn max_rt_by_file(features: &[Feature], n_files: usize) -> Vec<f64> {
    let max_rt = (0..n_files)
        .map(|_| AtomicU32::new(0))
        .map(|_| AtomicU32::new(0))
        .collect::<Vec<_>>();

    features.par_iter().for_each(|feat| {
        // If LFQ is performed, we should use MS1 apex RT, rather than PSM retention times
        let rt = if feat.ms1_apex_rt.is_normal() {
            feat.ms1_apex_rt
        } else {
            feat.rt
        };
        max_rt[feat.file_id].fetch_max(rt.ceil() as u32, std::sync::atomic::Ordering::Relaxed);
    });

    max_rt
        .into_iter()
        .map(|v| v.load(std::sync::atomic::Ordering::Acquire) as f64)
        .collect()
}

/// Return a map from PeptideIx to a map from File ID to average RT of the parent
/// PeptideIX
fn mean_rt_by_file(features: &[Feature]) -> FnvDashMap<PeptideIx, HashMap<usize, f64>> {
    let rts: FnvDashMap<PeptideIx, HashMap<usize, Vec<f64>>> = DashMap::default();
    features.par_iter().for_each(|feat| {
        // If LFQ is performed, we should use MS1 apex RT, rather than PSM retention times
        let rt = if feat.ms1_apex_rt.is_normal() {
            feat.ms1_apex_rt
        } else {
            feat.rt
        };
        rts.entry(feat.peptide_idx)
            .or_default()
            .entry(feat.file_id)
            .or_default()
            .push(rt as f64);
    });

    // Use the average RT for each file
    rts.into_par_iter()
        .map(|(ix, map)| {
            let map = map
                .into_iter()
                .map(|(file_id, rts)| (file_id, super::mean(&rts)))
                .collect();
            (ix, map)
        })
        .collect()
}

fn rt_matrix(features: &[Feature], max_rt: &[f64]) -> (HashMap<PeptideIx, f64>, Matrix) {
    let mean_rt = mean_rt_by_file(features);

    let (means, mat): (HashMap<PeptideIx, f64>, Vec<_>) = mean_rt
        .par_iter()
        .map(|entry| {
            let mut v = vec![f64::NAN; max_rt.len()];

            // While we're here, calculate the mean RT across all runs
            let mut sum = 0.0;
            let mut len = 0.0;
            for (&file_id, &rt) in entry.value() {
                let rt = rt / max_rt[file_id];
                v[file_id] = rt;
                sum += rt;
                len += 1.0;
            }

            ((*entry.key(), sum / len), v)
        })
        .filter(|((_, mean), _)| mean.is_normal())
        .unzip();
    let n = mat.len();
    let mat: Vec<f64> = mat.into_par_iter().flatten().collect();

    (means, Matrix::new(mat, n, max_rt.len()))
}

pub fn global_alignment(features: &mut [Feature], n_files: usize) {
    let max_rt = max_rt_by_file(features, n_files);
    let (peptide_rt, rt) = rt_matrix(features, &max_rt);

    let mean_rts: Vec<f64> = (0..rt.rows)
        .into_par_iter()
        .map(|row| {
            // Don't include NaN values
            let (len, sum) = rt
                .row(row)
                .filter(|rt| rt.is_finite())
                .fold((0, 0.0f64), |(len, sum), x| (len + 1, sum + x));
            sum / len as f64
        })
        .collect();

    // for file_idx in 0..n_files {
    let reg = (0..n_files)
        .into_par_iter()
        .map(|file_idx| {
            // calculate dot product across all ID'ed peptides
            let (len, dot, sum_x, sum_y) = rt
                .col(file_idx)
                .zip(mean_rts.iter())
                .filter(|(x, _)| x.is_finite())
                .fold(
                    (0, 0.0f64, 0.0f64, 0.0f64),
                    |(len, dot, sum_x, sum_y), (x, y)| (len + 1, dot + x * y, sum_x + x, sum_y + y),
                );

            let x_mean = sum_x / len as f64;
            let y_mean = sum_y / len as f64;
            let ssxy = dot - len as f64 * x_mean * y_mean;

            let sx2 = rt
                .col(file_idx)
                .filter(|rt| rt.is_finite())
                .fold(0.0f64, |sum, x| sum + (x - x_mean).powi(2));

            let slope = ssxy / sx2;
            let intercept = y_mean - slope * x_mean;

            log::info!(
                "aligning file #{file}: y = {m:.4}x + {b:.4}",
                file = file_idx,
                m = slope,
                b = intercept
            );

            (slope, intercept)
        })
        .collect::<Vec<(f64, f64)>>();

    log::info!("aligned retention times across {} files", n_files);

    features.par_iter_mut().for_each(|feature| {
        let (slope, intercept) = reg[feature.file_id];

        // Calculate aligned RT
        // - Divide by maximum RT of this run
        // - Multiply by regression parameters
        feature.aligned_rt =
            (feature.rt / max_rt[feature.file_id] as f32) * (slope as f32) + (intercept as f32);

        if let Some(&rt) = peptide_rt.get(&feature.peptide_idx) {
            feature.predicted_rt = rt as f32;
            feature.delta_rt = (rt as f32 - feature.aligned_rt).abs();
        }
    });
}
