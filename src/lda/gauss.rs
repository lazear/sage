//! Gauss-Jordan elimination for solution of systems of linear equations
//!
//! LDA requires solving the generalized eigenvalue problem for scatter matrices
//! Sb and Sw. We can actually solve this as the standard eigenvalue problem for
//! the matrix inv(Sw).dot(Sb) - or we solve the linear sysem Sw.dot(x) = Sb,
//! then calculate the eigenvalue for x. This is the approach we take.

use super::*;

#[derive(Debug)]
pub struct Gauss {
    pub left: Matrix,
    pub right: Matrix,
}

impl Matrix {
    fn swap_rows(&mut self, i: usize, j: usize) {
        for k in 0..self.cols {
            let tmp = self[(i, k)];
            self[(i, k)] = self[(j, k)];
            self[(j, k)] = tmp;
        }
    }
}

impl Gauss {
    pub fn solve(left: Matrix, right: Matrix) -> Option<Matrix> {
        let mut g = Gauss { left, right };
        g.echelon();
        g.reduce();
        g.backfill();

        // let eye = Matrix::identity(g.left.rows);

        // If `left` is the identity matrix, then `right` contains
        // the solution to the system of equations
        // match g.left.is_close(&eye, 0.00001) {
        match g.left_solved() {
            true => Some(g.right),
            false => None,
        }
    }

    // Is `left` an identity matrix, or else contains rows of all zeros?
    fn left_solved(&self) -> bool {
        let n = self.left.cols;
        for i in 0..n {
            for j in 0..n {
                let x = self.left[(i, j)];
                if i == j {
                    if x != 1.0 && x != 0.0 {
                        return false;
                    }
                } else if x != 0.0 {
                    return false;
                }
            }
        }
        true
    }

    fn echelon(&mut self) {
        let (m, n) = self.left.shape();
        let mut h = 0;
        let mut k = 0;

        while h < m && k < n {
            // find the row with the largest value in the current pivot column (k)
            let mut max = (0, f64::MIN);
            for i in h..m {
                if self.left[(i, k)] >= max.1 {
                    max = (i, self.left[(i, k)])
                }
            }
            let i = max.0;
            if self.left[(i, k)] == 0.0 {
                k += 1;
                continue;
            }

            // Swap rows (partial pivoting)
            if h != max.0 {
                self.left.swap_rows(h, i);
                self.right.swap_rows(h, i);
            }

            // Clear rows below pivot row
            for i in h + 1..m {
                let factor = self.left[(i, k)] / self.left[(h, k)];
                self.left[(i, k)] = 0.0;
                for j in k + 1..n {
                    self.left[(i, j)] -= self.left[(h, j)] * factor;
                }
                for j in 0..self.right.cols {
                    self.right[(i, j)] -= self.right[(h, j)] * factor;
                }
            }
            h += 1;
            k += 1;
        }
    }

    // Reduce left matrix to reduced echelon form - diagonal is all ones
    fn reduce(&mut self) {
        for i in (0..self.left.rows).rev() {
            for j in 0..self.left.cols {
                let x = self.left[(i, j)];
                if x == 0.0 {
                    continue;
                }
                for k in j..self.left.cols {
                    self.left[(i, k)] /= x;
                }
                for k in 0..self.right.cols {
                    self.right[(i, k)] /= x;
                }
                break;
            }
        }
    }

    // Solve the upper triangular matrix
    fn backfill(&mut self) {
        for i in (0..self.left.rows).rev() {
            for j in 0..self.left.cols {
                if self.left[(i, j)] == 0.0 {
                    continue;
                }
                for k in 0..i {
                    let factor = self.left[(k, j)] / self.left[(i, j)];
                    for h in 0..self.left.cols {
                        self.left[(k, h)] -= self.left[(i, h)] * factor;
                    }
                    for h in 0..self.right.cols {
                        self.right[(k, h)] -= self.right[(i, h)] * factor;
                    }
                }
                break;
            }
        }
    }
}
