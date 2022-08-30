//! Linear Discriminant Analysis for FDR refinement
//!
//! "What I cannot create, I do not understand" - Richard Feynman
//!
//! One of the major reasons for the creation of Sage is to develop a search
//! engine from first principles - And when I mean first principles, I mean
//! first principles - we are going to implement a basic linear algebra system
//! (complete with Gauss-Jordan elimination and eigenvector calculation) from scratch
//! to enable LDA.

mod gauss;
mod kde;
mod matrix;
use matrix::Matrix;
use rayon::prelude::*;

use crate::scoring::Percolator;

#[allow(dead_code)]
fn all_close(lhs: &[f64], rhs: &[f64], eps: f64) -> bool {
    lhs.iter()
        .zip(rhs.iter())
        .all(|(l, r)| (l - r).abs() <= eps)
}

pub fn norm(slice: &[f64]) -> f64 {
    slice.iter().fold(0.0, |acc, x| acc + x.powi(2)).sqrt()
}

pub fn mean(slice: &[f64]) -> f64 {
    slice.iter().sum::<f64>() / slice.len() as f64
}

pub fn std(slice: &[f64]) -> f64 {
    let mean = mean(slice);
    let x = slice.iter().fold(0.0, |acc, x| acc + (x - mean).powi(2));
    (x / slice.len() as f64).sqrt()
}

pub struct LinearDiscriminantAnalysis {
    eigenvector: Vec<f64>,
}

impl LinearDiscriminantAnalysis {
    pub fn train(features: &Matrix, decoy: &[bool]) -> Option<LinearDiscriminantAnalysis> {
        assert_eq!(features.rows, decoy.len());

        // Calculate class means, and overall mean
        let x_bar = features.mean();
        let mut scatter_within = Matrix::zeros(features.cols, features.cols);
        let mut scatter_between = Matrix::zeros(features.cols, features.cols);

        let mut class_means = Vec::new();

        for class in [true, false] {
            let count = decoy.iter().filter(|&label| *label == class).count();

            let class_data = (0..features.rows)
                .into_iter()
                .zip(decoy)
                .filter(|&(_, label)| *label == class)
                .flat_map(|(row, _)| features.row(row))
                .collect::<Vec<_>>();

            let mut class_data = Matrix::new(class_data, count, features.cols);
            let class_mean = class_data.mean();

            for row in 0..class_data.rows {
                for col in 0..class_data.cols {
                    class_data[(row, col)] -= class_mean[col];
                }
            }

            let cov = class_data.transpose().dot(&class_data);
            scatter_within += cov;

            let diff = Matrix::col_vector(
                class_mean
                    .iter()
                    .zip(x_bar.iter())
                    .map(|(x, y)| x - y)
                    .collect::<Vec<_>>(),
            );

            scatter_between += diff.dot(&diff.transpose());
            class_means.extend(class_mean);
        }

        let class_means = Matrix::new(class_means, 2, features.cols);

        // Use overall mean as the initial vector for power method... seems
        // unlikely to be the actual best eigenvector!
        let mut evec = gauss::Gauss::solve(scatter_within, scatter_between)
            .map(|mat| mat.power_method(&x_bar))?;

        // In some cases, power method can return eigenvector with signs flipped
        // Make it so that Target class scores are higher than Decoy, so that
        // we can make assumptions about this for ranking
        let coef = class_means.dotv(&evec);
        if coef[1] < coef[0] {
            evec.iter_mut().for_each(|c| *c *= -1.0);
        }

        log::trace!("linear model fit with eigenvector: {:?}", evec);

        Some(LinearDiscriminantAnalysis { eigenvector: evec })
    }

    pub fn score(&self, features: &Matrix) -> Vec<f64> {
        features.dotv(&self.eigenvector)
    }
}

pub fn score_psms(scores: &mut [Percolator]) -> Option<()> {
    log::trace!("fitting linear discriminant model");

    // Declare, so that we have compile time checking of matrix dimensions
    const FEATURES: usize = 13;
    let features = scores
        .into_par_iter()
        .flat_map(|perc| {
            let poisson = match -perc.poisson.log10() {
                x if x.is_finite() => x.ln_1p(),
                _ => 3.5,
            };

            let x: [f64; FEATURES] = [
                (perc.hyperscore.min(255.) as f64).ln_1p(),
                (perc.delta_hyperscore.min(255.) as f64).ln_1p(),
                (perc.delta_mass as f64).ln_1p(),
                perc.average_ppm as f64,
                poisson as f64,
                (perc.matched_intensity_pct as f64 * 100.0).ln_1p(),
                (perc.matched_peaks as f64).ln_1p(),
                // (perc.matched_neutral_loss as f64).ln_1p(),
                (perc.longest_b as f64).ln_1p(),
                (perc.longest_y as f64).ln_1p(),
                (perc.peptide_len as f64).ln_1p(),
                (perc.scored_candidates as f64).ln_1p(),
                (perc.rt as f64).ln_1p(),
                (perc.charge as f64).ln_1p(),
            ];
            x
        })
        .collect::<Vec<_>>();

    let decoys = scores.iter().map(|sc| sc.label == -1).collect::<Vec<_>>();
    let features = Matrix::new(features, scores.len(), FEATURES);

    let lda = LinearDiscriminantAnalysis::train(&features, &decoys)?;
    let discriminants = lda.score(&features);

    log::trace!("fitting non-parametric model for posterior error probabilities");
    let kde = kde::Estimator::fit(&discriminants, &decoys);

    scores
        .iter_mut()
        .zip(&discriminants)
        .for_each(|(perc, score)| {
            perc.discriminant_score = *score as f32;
            perc.posterior_error = kde.posterior_error(*score) as f32;
        });

    Some(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn linear_discriminant() {
        let a = Matrix::new([1., 2., 3., 4.], 2, 2);
        let eigenvector = [0.4159736, 0.90937671];
        assert!(all_close(
            &a.power_method(&[0.54, 0.34]),
            &eigenvector,
            1E-5
        ));

        #[rustfmt::skip]
        let feats = Matrix::new(
            [
                5., 4., 3., 2., 
                4., 5., 4., 3., 
                6., 3., 4., 5., 
                1., 0., 2., 9., 
                5., 4., 4., 3., 
                2., 1., 1., 9.5, 
                1., 0., 2., 8., 
                3., 2., -2., 10.,
            ],
            8,
            4,
        );

        let lda = LinearDiscriminantAnalysis::train(
            &feats,
            &[false, false, false, true, false, true, true, true],
        )
        .expect("error training LDA");

        let mut scores = lda.score(&feats);
        let norm = norm(&scores);
        scores = scores.into_iter().map(|s| s / norm).collect();

        let expected = [
            0.49706043,
            0.48920177,
            0.48920177,
            -0.07209359,
            0.51204672,
            -0.02849527,
            -0.04924864,
            -0.06055943,
        ];

        assert!(
            all_close(&scores, &expected, 1E-8),
            "{:?} {:?}",
            scores,
            expected
        );
    }
}
