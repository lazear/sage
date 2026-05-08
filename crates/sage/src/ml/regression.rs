//! Streaming OLS linear regression.
//!
//! Fits `beta = (X^T X)^-1 X^T y` without materializing the `n x D` design
//! matrix: one parallel fold/reduce pass accumulates `X^T X`, `X^T y`,
//! `sum(y)`, `sum(y^2)`, `n` from per-row outer products. A second pass
//! evaluates `sum((X beta - y)^2)` to report r^2 that matches the
//! materialized formula numerically.
//!
//! Per-worker scratch is `O(D^2)` (the cov accumulator), independent of `n`.

use super::{gauss::Gauss, matrix::Matrix};
use rayon::prelude::*;

pub struct LinearRegression {
    pub beta: Vec<f64>,
    pub r2: f64,
}

struct Acc {
    cov: Vec<f64>, // D*D row-major
    b: Vec<f64>,   // D
    sum_y: f64,
    sum_y2: f64,
    n: usize,
}

impl Acc {
    fn zero(d: usize) -> Self {
        Self {
            cov: vec![0.0; d * d],
            b: vec![0.0; d],
            sum_y: 0.0,
            sum_y2: 0.0,
            n: 0,
        }
    }

    fn add_row(&mut self, row: &[f64], y: f64) {
        let d = row.len();
        for j in 0..d {
            let rj = row[j];
            self.b[j] += rj * y;
            let off = j * d;
            for k in 0..d {
                self.cov[off + k] += rj * row[k];
            }
        }
        self.sum_y += y;
        self.sum_y2 += y * y;
        self.n += 1;
    }

    fn merge(mut self, other: Acc) -> Acc {
        for i in 0..self.cov.len() {
            self.cov[i] += other.cov[i];
        }
        for i in 0..self.b.len() {
            self.b[i] += other.b[i];
        }
        self.sum_y += other.sum_y;
        self.sum_y2 += other.sum_y2;
        self.n += other.n;
        self
    }
}

impl LinearRegression {
    /// Fit OLS over `items` with predicate `filter`. `embed(item)` produces a
    /// design row of length `D`; `target(item)` produces the response.
    ///
    /// Returns `None` if no items pass the filter or `X^T X` is singular.
    pub fn fit<T: Sync, const D: usize>(
        items: &[T],
        filter: impl Fn(&T) -> bool + Sync,
        embed: impl Fn(&T) -> [f64; D] + Sync,
        target: impl Fn(&T) -> f64 + Sync,
    ) -> Option<Self> {
        let acc = items
            .par_iter()
            .filter(|x| filter(x))
            .fold(
                || Acc::zero(D),
                |mut acc, x| {
                    let row = embed(x);
                    acc.add_row(&row, target(x));
                    acc
                },
            )
            .reduce(|| Acc::zero(D), Acc::merge);

        if acc.n == 0 {
            return None;
        }

        let nf = acc.n as f64;
        let y_mean = acc.sum_y / nf;
        let y_var = acc.sum_y2 - nf * y_mean * y_mean;

        let cov = Matrix::new(acc.cov, D, D);
        let b_mat = Matrix::col_vector(acc.b);
        let beta = Gauss::solve(cov, b_mat)?.take();

        // Streaming pass for SSE = sum((X beta - y)^2). O(N*D), small vs fit.
        let sse: f64 = items
            .par_iter()
            .filter(|x| filter(x))
            .map(|x| {
                let row = embed(x);
                let pred: f64 = row.iter().zip(&beta).map(|(v, w)| v * w).sum();
                let act = target(x);
                (pred - act).powi(2)
            })
            .sum();

        let r2 = 1.0 - sse / y_var;
        Some(Self { beta, r2 })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn fit_perfect_line() {
        // y = 2 x + 1, with intercept embedded as the last column.
        let items: Vec<(f64, f64)> = (0..50).map(|i| (i as f64, 2.0 * i as f64 + 1.0)).collect();
        let lr = LinearRegression::fit::<_, 2>(&items, |_| true, |&(x, _)| [x, 1.0], |&(_, y)| y)
            .unwrap();
        assert!((lr.beta[0] - 2.0).abs() < 1e-9, "slope: {}", lr.beta[0]);
        assert!((lr.beta[1] - 1.0).abs() < 1e-9, "intercept: {}", lr.beta[1]);
        assert!((lr.r2 - 1.0).abs() < 1e-9, "r2: {}", lr.r2);
    }

    #[test]
    fn fit_with_noise() {
        // y ~= 3 x + 2 with a deterministic perturbation; r^2 should be high.
        let items: Vec<(f64, f64)> = (0..200)
            .map(|i| {
                let x = i as f64 / 10.0;
                let noise = ((i as f64) * 0.7).sin() * 0.1;
                (x, 3.0 * x + 2.0 + noise)
            })
            .collect();
        let lr = LinearRegression::fit::<_, 2>(&items, |_| true, |&(x, _)| [x, 1.0], |&(_, y)| y)
            .unwrap();
        assert!((lr.beta[0] - 3.0).abs() < 0.05, "slope: {}", lr.beta[0]);
        assert!((lr.beta[1] - 2.0).abs() < 0.1, "intercept: {}", lr.beta[1]);
        assert!(lr.r2 > 0.99, "r2: {}", lr.r2);
    }

    #[test]
    fn empty_filter_returns_none() {
        let items: Vec<f64> = vec![1.0, 2.0, 3.0];
        let lr = LinearRegression::fit::<_, 1>(&items, |_| false, |_| [1.0], |&y| y);
        assert!(lr.is_none());
    }
}
