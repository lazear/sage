use crate::{
    mass::{Mass, H2O},
    peptide::Peptide,
};

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub enum Kind {
    B,
    Y,
}

/// Theoretical B/Y ion
#[derive(Copy, Clone, Debug)]
pub struct Ion {
    /// B or Y ion
    pub kind: Kind,
    /// Neutral fragment mass (no charge)
    pub monoisotopic_mass: f32,
}

/// Generate B/Y ions for a candidate peptide under a given charge state
pub struct IonSeries<'p> {
    pub kind: Kind,
    peptide: &'p Peptide,
    idx: usize,
}

impl<'p> IonSeries<'p> {
    /// Create a new [`IonSeries`] iterator for a specified peptide
    pub fn new(peptide: &'p Peptide, kind: Kind) -> Self {
        Self {
            kind,
            peptide,
            idx: 1,
        }
    }
}

impl<'p> Iterator for IonSeries<'p> {
    type Item = Ion;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == self.peptide.sequence.len() {
            return None;
        }
        let seq = match self.kind {
            Kind::B => self.peptide.sequence.get(..self.idx),
            Kind::Y => self.peptide.sequence.get(self.idx..),
        }?;
        self.idx += 1;

        let mut monoisotopic_mass: f32 = seq.iter().map(|r| r.monoisotopic()).sum();

        match (self.kind, self.peptide.nterm) {
            (Kind::B, None) => {}
            (Kind::B, Some(m)) => {
                monoisotopic_mass += m;
            }
            (Kind::Y, _) => {
                monoisotopic_mass += H2O;
            }
        }

        Some(Ion {
            kind: self.kind,
            monoisotopic_mass,
        })
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use super::*;
    use crate::{
        mass::{Residue, PROTON},
        peptide::Peptide,
    };

    impl FromStr for Peptide {
        type Err = char;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let mut sequence = Vec::with_capacity(s.len());
            let mut monoisotopic = H2O;

            for c in s.chars() {
                monoisotopic += c.monoisotopic();
                sequence.push(Residue::Just(c));
            }
            Ok(Peptide {
                protein: String::default(),
                sequence,
                nterm: None,
                monoisotopic,
            })
        }
    }

    fn check_within<I: Iterator<Item = Ion>>(iter: I, expected_mz: &[f32]) {
        let observed = iter.map(|ion| ion.monoisotopic_mass).collect::<Vec<f32>>();
        assert_eq!(expected_mz.len(), observed.len());
        assert!(expected_mz
            .iter()
            .zip(observed.iter())
            .all(|(a, b)| (a - b).abs() < 0.001));
    }

    macro_rules! ions {
        ($peptide:expr, $kind:expr, $charge:expr) => {{
            IonSeries::new($peptide, $kind).map(|mut ion| {
                ion.monoisotopic_mass = (ion.monoisotopic_mass + $charge * PROTON) / $charge;
                ion
            })
        }};
    }

    #[test]
    fn iterate_b_ions() {
        let peptide = "PEPTIDE".parse::<Peptide>().unwrap();

        // Charge state 2
        let expected_mz = vec![
            98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
        ];

        check_within(ions!(&peptide, Kind::B, 1.0), &expected_mz);
    }

    #[test]
    fn iterate_y_ions() {
        let peptide = "PEPTIDE".parse::<Peptide>().unwrap();

        // Charge state 1
        let expected_mz = vec![
            703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
        ];

        check_within(ions!(&peptide, Kind::Y, 1.0), &expected_mz);
    }

    #[test]
    fn decoy() {
        let peptide = "PEPTIDE".parse::<Peptide>().unwrap();

        // Charge state 2
        let expected_mz = vec![
            352.16087, 287.639_6, 239.11319, 188.58935, 132.04732, 74.53385,
        ];

        check_within(ions!(&peptide, Kind::Y, 2.0), &expected_mz);

        let peptide = "EDITPEP".parse::<Peptide>().unwrap();

        // Charge state 2
        let expected_mz = vec![
            336.16596, 278.652_5, 222.110_46, 171.586_62, 123.060237, 58.538_94,
        ];

        check_within(ions!(&peptide, Kind::Y, 2.0), &expected_mz);
    }

    #[test]
    fn nterm_mod() {
        let mut peptide = "PEPTIDE".parse::<Peptide>().unwrap();
        peptide.set_nterm_mod(229.01);

        // Charge state 1, TMT tagged
        let expected_b = [
            98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
        ].into_iter().map(|x| x + 229.01).collect::<Vec<_>>();

        let expected_y = vec![
            703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
        ];

        check_within(ions!(&peptide, Kind::B, 1.0), &expected_b);
        check_within(ions!(&peptide, Kind::Y, 1.0), &expected_y);

    }
}
