//! Linear Algebra, Machine Learning & FDR refinement

pub mod gauss;
pub mod kde;
pub mod linear_discriminant;
pub mod matrix;
pub mod mobility_model;
pub mod qvalue;
pub mod retention_alignment;
pub mod retention_model;

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
