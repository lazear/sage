use super::Matrix;

use rayon::prelude::*;
use std::fmt::{self, Debug};
use std::ops::{Add, AddAssign, Deref, DerefMut, Div, Index, IndexMut, Mul, Sub};

impl Mul<&[f64]> for &Matrix {
    type Output = Vec<f64>;

    fn mul(self, rhs: &[f64]) -> Self::Output {
        self.dotv(rhs)
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

    fn add(self, rhs: Matrix) -> Self::Output {
        assert_eq!(
            self.shape(),
            rhs.shape(),
            "matrices must have equal shape to add"
        );
        let mut output = self.clone();
        for i in 0..output.data.len() {
            output.data[i] += rhs.data[i];
        }
        output
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

/// Parallel matrix multiplication using Rayon
// macro_rules! impl_mul {
//     (($lhs:ty, $rhs:ty)) => {
//         impl Mul<$rhs> for $lhs {
//             type Output = Matrix;
//             fn mul(self, rhs: $rhs) -> Self::Output {
//                 assert_eq!(self.cols, rhs.rows, "matrices must match M,N * N,P shape");
//                 let lhs = &self;
//                 let rhs = &rhs;
//                 let output = (0..self.rows)
//                     .into_par_iter()
//                     .flat_map(|row| {
//                         (0..rhs.cols).into_par_iter().map(move |col| {
//                             lhs.row(row)
//                                 .zip(rhs.col(col))
//                                 .fold(0.0, |acc, (x, y)| acc + x * y)
//                         })
//                     })
//                     .collect::<Vec<_>>();

//                 Matrix {
//                     data: output,
//                     rows: self.rows,
//                     cols: rhs.cols,
//                 }
//             }
//         }
//     };
// }

macro_rules! impl_scalar {
    ($lhs:ty) => {
        impl Add<f64> for $lhs {
            type Output = Matrix;
            fn add(self, rhs: f64) -> Self::Output {
                let mut output = self.clone();
                output.data.iter_mut().for_each(|x| *x += rhs);
                output
            }
        }

        impl Sub<f64> for $lhs {
            type Output = Matrix;
            fn sub(self, rhs: f64) -> Self::Output {
                let mut output = self.clone();
                output.data.iter_mut().for_each(|x| *x -= rhs);
                output
            }
        }

        impl Mul<f64> for $lhs {
            type Output = Matrix;
            fn mul(self, rhs: f64) -> Self::Output {
                let mut output = self.clone();
                output.data.iter_mut().for_each(|x| *x *= rhs);
                output
            }
        }

        impl Div<f64> for $lhs {
            type Output = Matrix;
            fn div(self, rhs: f64) -> Self::Output {
                let mut output = self.clone();
                output.data.iter_mut().for_each(|x| *x /= rhs);
                output
            }
        }
    };
}

// impl_mul!((Matrix, Matrix));
// impl_mul!((Matrix, &Matrix));
// impl_mul!((&Matrix, Matrix));
// impl_mul!((&Matrix, &Matrix));
impl_scalar!(Matrix);
impl_scalar!(&Matrix);

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
    fn multiply() {
        let a = vec![1., 0., 1., 2., 1., 1., 0., 1., 1., 1., 1., 2.];
        let b = vec![1., 2., 1., 2., 3., 1., 4., 2., 2.];
        let a = Matrix {
            data: a,
            rows: 4,
            cols: 3,
        };
        let b = Matrix {
            data: b,
            rows: 3,
            cols: 3,
        };

        let c = a.dot(&b);
        assert_eq!(c.rows, 4);
        assert_eq!(c.cols, 3);
        assert_eq!(
            c.data,
            vec![5., 4., 3., 8., 9., 5., 6., 5., 3., 11., 9., 6.]
        );

        let d = vec![1., 2., 3., 4., 5., 6.];
        let d = Matrix {
            data: d,
            rows: 2,
            cols: 3,
        };
        let e = Matrix::row_vector(vec![7., 9., 11.]);

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
