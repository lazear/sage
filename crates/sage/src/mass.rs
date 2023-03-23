use std::{fmt::Write, iter::Sum};

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

    pub fn ppm_to_delta_mass(center: f32, ppm: f32) -> f32 {
        ppm * center / 1_000_000.0
    }
}

pub trait Mass {
    fn monoisotopic(&self) -> f32;
    fn composition(&self) -> Composition;
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize)]
pub enum Residue {
    // Standard amino acid residue
    Just(char),
    // Amino acid residue with a mass modification
    Mod(char, f32),
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

pub const VALID_AA: [char; 22] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y', 'U', 'O',
];

impl Mass for char {
    fn monoisotopic(&self) -> f32 {
        match self {
            'A' => 71.03711,
            'R' => 156.1011,
            'N' => 114.04293,
            'D' => 115.02694,
            'C' => 103.00919,
            'E' => 129.04259,
            'Q' => 128.05858,
            'G' => 57.02146,
            'H' => 137.05891,
            'I' => 113.08406,
            'L' => 113.08406,
            'K' => 128.09496,
            'M' => 131.0405,
            'F' => 147.0684,
            'P' => 97.05276,
            'S' => 87.03203,
            'T' => 101.04768,
            'W' => 186.07931,
            'Y' => 163.06333,
            'V' => 99.06841,
            'U' => 150.95363,
            'O' => 237.14773,
            _ => unreachable!("BUG: invalid amino acid {}", self),
        }
    }

    fn composition(&self) -> Composition {
        match self {
            'A' => Composition::new(3, 2, 0),
            'R' => Composition::new(6, 2, 0),
            'N' => Composition::new(4, 3, 0),
            'D' => Composition::new(4, 4, 0),
            'C' => Composition::new(3, 2, 1),
            'E' => Composition::new(5, 4, 0),
            'Q' => Composition::new(5, 3, 0),
            'G' => Composition::new(2, 2, 0),
            'H' => Composition::new(6, 2, 0),
            'I' => Composition::new(6, 2, 0),
            'L' => Composition::new(6, 2, 0),
            'K' => Composition::new(6, 2, 0),
            'M' => Composition::new(5, 2, 1),
            'F' => Composition::new(9, 2, 0),
            'P' => Composition::new(5, 2, 0),
            'S' => Composition::new(3, 3, 0),
            'T' => Composition::new(4, 3, 0),
            'W' => Composition::new(11, 2, 0),
            'Y' => Composition::new(9, 3, 0),
            'V' => Composition::new(5, 2, 0),
            'U' => Composition::new(3, 2, 0),
            'O' => Composition::new(12, 3, 0),
            _ => unreachable!("BUG: invalid amino acid {}", self),
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
            Residue::Just(c) => f.write_char(*c),
            Residue::Mod(c, m) => {
                if m.is_sign_positive() {
                    write!(f, "{}[+{}]", c, m)
                } else {
                    write!(f, "{}[{}]", c, m)
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
