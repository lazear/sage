use super::norm;
use rayon::prelude::*;
use std::fmt::{self, Debug};
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Index, IndexMut};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Row;
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Col;

#[derive(Clone, PartialEq, PartialOrd)]
pub struct Matrix {
    data: Vec<f64>,
    pub rows: usize,
    pub cols: usize,
}

pub struct Iter<'a, Axes> {
    data: &'a Matrix,
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

impl Matrix {
    /// Create a new `Matrix`
    ///
    /// # Panics
    ///
    /// * Panics if `data` does not have len == rows * cols
    pub fn new<T: Into<Vec<f64>>>(t: T, rows: usize, cols: usize) -> Matrix {
        let data = t.into();
        assert_eq!(
            data.len(),
            rows * cols,
            "data passed to Matrix::new() does not have shape ({}, {})",
            rows,
            cols
        );
        Matrix { data, rows, cols }
    }

    pub fn zeros(rows: usize, cols: usize) -> Matrix {
        Matrix {
            data: vec![0.0; rows * cols],
            rows,
            cols,
        }
    }

    pub fn identity(size: usize) -> Matrix {
        let mut matrix = Matrix::zeros(size, size);
        for i in 0..size {
            matrix[(i, i)] = 1.0
        }
        matrix
    }

    pub fn col_vector(data: Vec<f64>) -> Matrix {
        let rows = data.len();
        Matrix {
            data,
            rows,
            cols: 1,
        }
    }

    pub fn row_vector(data: Vec<f64>) -> Matrix {
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
            data: self,
            row,
            col: 0,
            axes: PhantomData,
        }
    }

    pub fn col(&self, col: usize) -> Iter<'_, Col> {
        Iter {
            data: self,
            row: 0,
            col,
            axes: PhantomData,
        }
    }

    pub fn is_close(&self, rhs: &Self, eps: f64) -> bool {
        if self.cols != self.rows {
            return false;
        }

        for i in 0..self.rows {
            for j in 0..self.rows {
                if (self[(i, j)] - rhs[(i, j)]).abs() > eps {
                    return false;
                }
            }
        }
        true
    }

    // Use power method to find the eigenvector with the largest
    // corresponding eigenvalue
    pub fn power_method(&self, initial: &[f64]) -> Vec<f64> {
        // Pick a starting vector

        let n = norm(initial);
        let mut v = initial.iter().map(|i| i / n).collect::<Vec<_>>();

        let mut last_eig = 0.0;
        for _ in 0..50 {
            let mut v1 = self.dotv(&v);
            let norm = norm(&v1);
            if (norm - last_eig).abs() < 1E-8 {
                break;
            }
            last_eig = norm;
            v1.iter_mut().for_each(|x| *x /= norm);
            v = v1;
        }
        v
    }

    pub fn transpose(&self) -> Matrix {
        if self.cols == 1 || self.rows == 1 {
            let mut mat = self.clone();
            std::mem::swap(&mut mat.cols, &mut mat.rows);
            return mat;
        }
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
            self.cols,
            rhs.len(),
            "lhs has shape ({},{}), rhs has shape (1,{})",
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

impl Index<(usize, usize)> for Matrix {
    type Output = f64;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.data[self.cols * row + col]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.data[self.cols * row + col]
    }
}

impl Add<Matrix> for Matrix {
    type Output = Matrix;

    fn add(mut self, rhs: Matrix) -> Self::Output {
        assert_eq!(
            self.shape(),
            rhs.shape(),
            "matrices must have equal shape to add"
        );
        for i in 0..self.data.len() {
            self.data[i] += rhs.data[i];
        }
        self
    }
}

impl AddAssign<Matrix> for Matrix {
    fn add_assign(&mut self, rhs: Matrix) {
        assert_eq!(
            self.shape(),
            rhs.shape(),
            "matrices must have equal shape to add"
        );
        for i in 0..self.data.len() {
            self.data[i] += rhs.data[i];
        }
    }
}

impl Debug for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "[")?;
        for row in 0..self.rows {
            writeln!(f, "{:?}", self.row(row).collect::<Vec<_>>())?;
        }
        writeln!(f, "]")
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn dotv() {
        let a = Matrix::new([1., 2., 3., 4.], 2, 2);

        let v0 = a.dotv(&[0.5, 0.5]);
        assert_eq!(v0, vec![1.5, 3.5]);
        let n = norm(&v0);

        let c = v0.iter().map(|v| v / n).collect::<Vec<_>>();
        assert!(c
            .iter()
            .zip(&[0.3939193, 0.91914503])
            .all(|(x, y)| (x - y).abs() <= 0.0001));
    }

    #[test]
    fn tranpose() {
        let mut mat = Matrix {
            data: vec![1., 2., 3., 4., 5., 6.],
            rows: 3,
            cols: 2,
        };

        assert_eq!(mat[(0, 0)], 1., "{:?}", mat);
        assert_eq!(mat[(0, 1)], 2., "{:?}", mat);
        assert_eq!(mat[(1, 0)], 3., "{:?}", mat);
        assert_eq!(mat[(1, 1)], 4., "{:?}", mat);
        assert_eq!(mat[(2, 0)], 5., "{:?}", mat);
        assert_eq!(mat[(2, 1)], 6., "{:?}", mat);

        mat = mat.transpose();

        assert_eq!(mat[(0, 0)], 1., "{:?}", mat);
        assert_eq!(mat[(0, 1)], 3., "{:?}", mat);
        assert_eq!(mat[(0, 2)], 5., "{:?}", mat);
        assert_eq!(mat[(1, 0)], 2., "{:?}", mat);
        assert_eq!(mat[(1, 1)], 4., "{:?}", mat);
        assert_eq!(mat[(1, 2)], 6., "{:?}", mat);
    }

    #[test]
    fn dot() {
        #[rustfmt::skip]
        let a = vec![
            1., 0., 1., 
            2., 1., 1., 
            0., 1., 1., 
            1., 1., 2.
        ];
        let a = Matrix::new(a, 4, 3);

        #[rustfmt::skip]
        let b = vec![
            1., 2., 1., 
            2., 3., 1., 
            4., 2., 2.
        ];
        let b = Matrix::new(b, 3, 3);

        let c = a.dot(&b);
        assert_eq!(c.rows, 4);
        assert_eq!(c.cols, 3);
        #[rustfmt::skip]
        assert_eq!(
            c.data,
            vec![
                5., 4., 3., 
                8., 9., 5., 
                6., 5., 3., 
                11., 9., 6.
            ]
        );

        let d = vec![1., 2., 3., 4., 5., 6.];
        let d = Matrix::new(d, 2, 3);
        let e = Matrix::col_vector(vec![7., 9., 11.]);

        assert_eq!(
            d.dot(&e),
            Matrix {
                data: vec![58., 139.],
                cols: 1,
                rows: 2
            }
        );
    }
}
