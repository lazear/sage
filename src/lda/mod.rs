//! Linear Discriminant Analysis for FDR refinement
//!
//! "What I cannot create, I do not understand" - Richard Feynman
//!
//! One of the major reasons for the creation of Sage is to develop a search
//! engine from first principles - And when I mean first principles, I mean
//! first principles - we are going to implement a basic linear algebra system
//! (complete with Gauss-Jordan elimination and eigenvector calculation) from scratch
//! to enable LDA.

mod impls;

use rayon::prelude::*;
use std::fmt::{self, Debug};
use std::marker::PhantomData;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Row;
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Col;

#[derive(Clone, PartialEq, PartialOrd)]
pub struct Matrix {
    data: Vec<f64>,
    rows: usize,
    cols: usize,
}

pub struct Iter<'a, Axes> {
    data: &'a Matrix,
    row: usize,
    col: usize,
    axes: PhantomData<Axes>,
}

pub struct IterMut<'a, Axes> {
    data: &'a mut Matrix,
    row: usize,
    col: usize,
    axes: PhantomData<Axes>,
}

impl<'a> Iterator for Iter<'a, Row> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        let data = self.data.get(self.row, self.col);
        self.col += 1;
        data
    }
}

impl<'a> Iterator for Iter<'a, Col> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        let data = self.data.get(self.row, self.col);
        self.row += 1;
        data
    }
}

pub fn norm(slice: &[f64]) -> f64 {
    slice.iter().fold(0.0, |acc, x| acc + x.powi(2)).sqrt()
}

impl Matrix {
    pub fn new<T: Into<Vec<f64>>>(t: T, rows: usize, cols: usize) -> Matrix {
        Matrix {
            data: t.into(),
            rows,
            cols,
        }
    }

    pub fn zeros(rows: usize, cols: usize) -> Matrix {
        Matrix {
            data: vec![0.0; rows * cols],
            rows,
            cols,
        }
    }

    pub fn row_vector(data: Vec<f64>) -> Matrix {
        let rows = data.len();
        Matrix {
            data,
            rows,
            cols: 1,
        }
    }

    pub fn col_vector(data: Vec<f64>) -> Matrix {
        let cols = data.len();
        Matrix {
            data,
            rows: 1,
            cols,
        }
    }

    pub const fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row >= self.rows || col >= self.cols {
            None
        } else {
            self.data.get(self.cols * row + col).copied()
        }
    }

    pub fn get_mut(&mut self, row: usize, col: usize) -> Option<&mut f64> {
        self.data.get_mut(self.cols * row + col)
    }

    pub fn row(&self, row: usize) -> Iter<'_, Row> {
        Iter {
            data: &self,
            row,
            col: 0,
            axes: PhantomData,
        }
    }

    pub fn col(&self, col: usize) -> Iter<'_, Col> {
        Iter {
            data: &self,
            row: 0,
            col,
            axes: PhantomData,
        }
    }

    // Use power method to find the eigenvector with the largest
    // corresponding eigenvalue
    pub fn power_method(&self) -> Vec<f64> {
        // Pick a starting vector
        // let mut v = Vector::col(vec![0.3; self.cols]);
        let mut v = self.row(0).collect::<Vec<_>>();
        let mut last_eig = 0.0;
        for _ in 0..10 {
            let mut v1 = self.dotv(&v);
            let norm = norm(&v1);
            if (norm - last_eig).abs() < f32::EPSILON as f64 {
                break;
            }
            last_eig = norm;
            v1.iter_mut().for_each(|x| *x /= norm);
            v = v1;
        }
        v
    }

    pub fn transpose(&self) -> Matrix {
        let mut mat = Matrix::zeros(self.cols, self.rows);
        for row in 0..self.rows {
            for col in 0..self.cols {
                mat[(col, row)] = self[(row, col)]
            }
        }
        mat
    }

    pub fn dotv(&self, rhs: &[f64]) -> Vec<f64> {
        assert_eq!(
            self.rows,
            rhs.len(),
            "lhs has shape ({},{}), rhs has shape (,{})",
            self.rows,
            self.cols,
            rhs.len()
        );
        (0..self.rows)
            .into_par_iter()
            .map(|row| self.row(row).zip(rhs).fold(0.0, |acc, (x, y)| acc + x * y))
            .collect::<Vec<_>>()
    }

    pub fn dot(&self, rhs: &Matrix) -> Matrix {
        assert_eq!(
            self.cols, rhs.rows,
            "lhs has shape ({},{}), rhs has shape ({},{})",
            self.rows, self.cols, rhs.rows, rhs.cols
        );
        let data = (0..self.rows)
            .into_par_iter()
            .flat_map(|row| {
                (0..rhs.cols).into_par_iter().map(move |col| {
                    self.row(row)
                        .zip(rhs.col(col))
                        .fold(0.0, |acc, (x, y)| acc + x * y)
                })
            })
            .collect::<Vec<_>>();
        Matrix {
            data,
            rows: self.rows,
            cols: rhs.cols,
        }
    }

    /// Calculate mean of each column
    pub fn mean(&self) -> Vec<f64> {
        (0..self.cols)
            .into_par_iter()
            .map(|col| {
                let sum = self.col(col).sum::<f64>();
                sum / self.rows as f64
            })
            .collect()
    }
}

pub fn lda(features: Matrix, decoy: &[bool]) {
    assert_eq!(features.rows, decoy.len());
    // Calculate class means, and overall mean
    // let x_bar = features.mean(std::iter::repeat(true));
    let x_bar = features.mean();
    let mut scatter_within = Matrix::zeros(features.cols, features.cols);
    let mut scatter_between = Matrix::zeros(features.cols, features.cols);

    for class in [true, false] {
        let count = decoy.iter().filter(|&label| *label == class).count();

        // NB: Using par_bridge() removes ordering guarantees - should be fine
        // for this use though
        let class_data = (0..features.rows)
            .into_par_iter()
            .zip(decoy)
            .filter(|&(_, label)| *label == class)
            .flat_map(|(row, _)| features.row(row).par_bridge())
            .collect::<Vec<_>>();

        let mut class_data = Matrix::new(class_data, count, features.cols);
        let class_mean = class_data.mean();
        dbg!(&class_mean);

        for row in 0..class_data.rows {
            for col in 0..class_data.cols {
                class_data[(row, col)] -= class_mean[col];
            }
        }

        let cov = class_data.transpose().dot(&class_data);
        dbg!(&cov);
        scatter_within += cov;

        let diff = Matrix::new(
            class_mean
                .iter()
                .zip(x_bar.iter())
                .map(|(x, y)| x - y)
                .collect::<Vec<_>>(),
            features.cols,
            1,
        );
        let mut diff2 = diff.clone();
        std::mem::swap(&mut diff2.cols, &mut diff2.rows);
        scatter_between += diff.dot(&diff2);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn eigenvector() {
        let a = Matrix::new([1., 2., 3., 4.], 2, 2);
        let eigenvector = [0.4159736, 0.90937671];
        assert!(a
            .power_method()
            .iter()
            .zip(eigenvector)
            .all(|(x, y)| (x - y).abs() <= 0.0001));

        let feats = Matrix::new(
            [
                5., 4., 3., 2., 4., 5., 4., 3., 6., 3., 4., 5., 1., 0., 2., 9., 5., 4., 4., 3., 2.,
                1., 1., 9.5, 1., 0., 2., 8., 3., 2., -2., 10.,
            ],
            8,
            4,
        );
        lda(feats, &[false, false, false, true, false, true, true, true]);
        panic!("oops");
    }

    #[test]
    fn dot2() {
        let a = Matrix::new([1., 2., 3., 4.], 2, 2);
        // let c = &a * b;
        let v0 = a.dotv(&[0.5, 0.5]);
        assert_eq!(v0, vec![1.5, 3.5]);
        let n = norm(&v0);

        let c = v0.iter().map(|v| v / n).collect::<Vec<_>>();
        assert!(c
            .iter()
            .zip(&[0.3939193, 0.91914503])
            .all(|(x, y)| (x - y).abs() <= 0.0001));
    }
}
