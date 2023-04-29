use std::{fmt::Write, iter::Sum, ops::Mul};

use serde::{Deserialize, Serialize};

pub const H2O: f32 = 18.010565;
pub const PROTON: f32 = 1.0072764;
pub const NEUTRON: f32 = 1.00335;
pub const NH3: f32 = 17.026548;

#[derive(Copy, Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd)]
#[serde(rename_all = "lowercase")]
pub enum Tolerance {
    Ppm(f32, f32),
    Da(f32, f32),
}

impl Tolerance {
    /// Compute the (`lower`, `upper`) window (in Da) for for a monoisotopic
    /// mass and a given tolerance
    pub fn bounds(&self, center: f32) -> (f32, f32) {
        match self {
            Tolerance::Ppm(lo, hi) => {
                let delta_lo = center * lo / 1_000_000.0;
                let delta_hi = center * hi / 1_000_000.0;
                (center + delta_lo, center + delta_hi)
            }
            Tolerance::Da(lo, hi) => (center + lo, center + hi),
        }
    }

    pub fn contains(&self, center: f32, rhs: f32) -> bool {
        let (lo, hi) = self.bounds(center);
        rhs >= lo && rhs <= hi
    }

    pub fn ppm_to_delta_mass(center: f32, ppm: f32) -> f32 {
        ppm * center / 1_000_000.0
    }
}

impl Mul<f32> for Tolerance {
    type Output = Tolerance;

    fn mul(self, rhs: f32) -> Self::Output {
        match self {
            Tolerance::Ppm(lo, hi) => Tolerance::Ppm(lo * rhs, hi * rhs),
            Tolerance::Da(lo, hi) => Tolerance::Da(lo * rhs, hi * rhs),
        }
    }
}

pub trait Mass {
    fn monoisotopic(&self) -> f32;
    fn composition(&self) -> Composition;
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize)]
pub enum Residue {
    // Standard amino acid residue
    Just(u8),
    // Amino acid residue with a mass modification
    Mod(u8, f32),
}

impl Mass for Residue {
    fn monoisotopic(&self) -> f32 {
        match self {
            Residue::Just(c) => c.monoisotopic(),
            Residue::Mod(c, m) => c.monoisotopic() + m,
        }
    }

    fn composition(&self) -> Composition {
        match self {
            Residue::Just(c) => c.composition(),
            Residue::Mod(c, _) => c.composition(),
        }
    }
}

pub const VALID_AA: [u8; 22] = [
    b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S',
    b'T', b'V', b'W', b'Y', b'U', b'O',
];

impl Mass for u8 {
    fn monoisotopic(&self) -> f32 {
        match self {
            b'A' => 71.03711,
            b'R' => 156.1011,
            b'N' => 114.04293,
            b'D' => 115.02694,
            b'C' => 103.00919,
            b'E' => 129.04259,
            b'Q' => 128.05858,
            b'G' => 57.02146,
            b'H' => 137.05891,
            b'I' => 113.08406,
            b'L' => 113.08406,
            b'K' => 128.09496,
            b'M' => 131.0405,
            b'F' => 147.0684,
            b'P' => 97.05276,
            b'S' => 87.03203,
            b'T' => 101.04768,
            b'W' => 186.07931,
            b'Y' => 163.06333,
            b'V' => 99.06841,
            b'U' => 150.95363,
            b'O' => 237.14773,
            _ => unreachable!("BUG: invalid amino acid {}", *self as char),
        }
    }

    fn composition(&self) -> Composition {
        match self {
            b'A' => Composition::new(3, 2, 0),
            b'R' => Composition::new(6, 2, 0),
            b'N' => Composition::new(4, 3, 0),
            b'D' => Composition::new(4, 4, 0),
            b'C' => Composition::new(3, 2, 1),
            b'E' => Composition::new(5, 4, 0),
            b'Q' => Composition::new(5, 3, 0),
            b'G' => Composition::new(2, 2, 0),
            b'H' => Composition::new(6, 2, 0),
            b'I' => Composition::new(6, 2, 0),
            b'L' => Composition::new(6, 2, 0),
            b'K' => Composition::new(6, 2, 0),
            b'M' => Composition::new(5, 2, 1),
            b'F' => Composition::new(9, 2, 0),
            b'P' => Composition::new(5, 2, 0),
            b'S' => Composition::new(3, 3, 0),
            b'T' => Composition::new(4, 3, 0),
            b'W' => Composition::new(11, 2, 0),
            b'Y' => Composition::new(9, 3, 0),
            b'V' => Composition::new(5, 2, 0),
            b'U' => Composition::new(3, 2, 0),
            b'O' => Composition::new(12, 3, 0),
            _ => unreachable!("BUG: invalid amino acid {}", *self as char),
        }
    }
}

pub struct Composition {
    pub carbon: u16,
    pub sulfur: u16,
}

impl Composition {
    pub const fn new(carbon: u16, _oxygen: u8, sulfur: u16) -> Self {
        Self { carbon, sulfur }
    }
}

impl Sum for Composition {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut comp = Composition::new(0, 0, 0);
        for i in iter {
            comp.carbon += i.carbon;
            // comp.oxygen += i.oxygen;
            comp.sulfur += i.sulfur;
        }
        comp
    }
}

impl std::fmt::Display for Residue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Residue::Just(c) => f.write_char(*c as char),
            Residue::Mod(c, m) => {
                if m.is_sign_positive() {
                    write!(f, "{}[+{}]", *c as char, m)
                } else {
                    write!(f, "{}[{}]", *c as char, m)
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::{Mass, Tolerance, VALID_AA};

    #[test]
    fn smoke() {
        for ch in VALID_AA {
            assert!(ch.monoisotopic() > 0.0);
        }
    }

    #[test]
    fn tolerances() {
        assert_eq!(
            Tolerance::Ppm(-10.0, 20.0).bounds(1000.0),
            (999.99, 1000.02)
        );
        assert_eq!(
            Tolerance::Ppm(-10.0, 10.0).bounds(487.0),
            (486.99513, 487.00487)
        );
        assert_eq!(
            Tolerance::Ppm(-50.0, 50.0).bounds(1000.0),
            (999.95, 1000.05)
        );
    }
}
