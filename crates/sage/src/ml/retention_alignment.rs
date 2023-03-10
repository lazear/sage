use std::collections::HashMap;
use std::hash::BuildHasherDefault;

use super::matrix::Matrix;
use crate::database::PeptideIx;
use crate::scoring::Feature;
use dashmap::DashMap;
use fnv::FnvHasher;
use rayon::prelude::*;

type FnvDashMap<K, V> = DashMap<K, V, BuildHasherDefault<FnvHasher>>;

fn mean_rt_by_file(features: &[Feature]) -> FnvDashMap<PeptideIx, HashMap<usize, f64>> {
    let rts: FnvDashMap<PeptideIx, HashMap<usize, Vec<f64>>> = DashMap::default();
    // let max_rt: FnvDashMap<usize, f64> = DashMap::default();

    features.par_iter().for_each(|feat| {
        // If LFQ is performed, we should use MS1 apex RT, rather than PSM retention times
        if feat.ms1_apex_rt > 0.0 {
            rts.entry(feat.peptide_idx)
                .or_default()
                .entry(feat.file_id)
                .or_default()
                .push(feat.ms1_apex_rt as f64);
        } else {
            rts.entry(feat.peptide_idx)
                .or_default()
                .entry(feat.file_id)
                .or_default()
                .push(feat.rt as f64);
        }

        // let mut curr_max = max_rt.entry(feat.file_id).or_default();
        // *curr_max = curr_max.max(feat.rt as f64);
    });

    // Use the average RT for each file
    rts.into_par_iter()
        .map(|(ix, map)| {
            let map = map
                .into_iter()
                .map(|(file_id, rts)| {
                    // let file_max_rt = max_rt.get(&file_id).map(|e| e.value()).unwrap_or(200.0);

                    (file_id, super::mean(&rts))
                })
                .collect();
            (ix, map)
        })
        .collect()
}

fn rt_matrix(features: &[Feature], n_files: usize) -> Matrix {
    let mean_rt = mean_rt_by_file(features);

    let mat = mean_rt
        .par_iter()
        .filter(|entry| entry.value().len() >= 2)
        .flat_map(|entry| {
            let mut v = vec![f64::NAN; n_files];
            for (file_id, rt) in entry.value() {
                v[*file_id] = *rt;
            }
            v
        })
        .collect::<Vec<_>>();

    Matrix::new(mat, mean_rt.len(), n_files)
}

pub fn global_alignment(features: &mut [Feature], n_files: usize) {
    if n_files < 2 {
        return;
    }

    let rt = rt_matrix(features, n_files);

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
                "running linear regression {file}: y = {m}x+{b}",
                file = file_idx,
                m = slope,
                b = intercept
            );

            (slope, intercept)
        })
        .collect::<Vec<(f64, f64)>>();

    features.par_iter_mut().for_each(|feature| {
        let (slope, intercept) = reg[feature.file_id];
        feature.rt = feature.rt * (slope as f32) + (intercept as f32);
    });
}
