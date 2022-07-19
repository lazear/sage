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
    pub kind: Kind,
    pub mz: f32,
}

/// Generate B/Y ions for a candidate peptide under a given charge state
pub struct IonSeries<'p> {
    pub kind: Kind,
    peptide: &'p Peptide,
    idx: usize,
}

impl<'p> IonSeries<'p> {
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

        let mut neutral: f32 = seq.iter().map(|r| r.monoisotopic()).sum();
        if Kind::Y == self.kind {
            neutral += H2O;
        }
        // neutral = (neutral + self.charge as f32 * PROTON) / self.charge as f32;

        Some(Ion {
            kind: self.kind,
            mz: neutral,
        })
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use super::*;
    use crate::{mass::Residue, peptide::Peptide};

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
                monoisotopic,
            })
        }
    }

    fn check_within(iter: IonSeries, expected_mz: &[f32]) {
        let observed = iter.map(|ion| ion.mz).collect::<Vec<f32>>();
        assert_eq!(expected_mz.len(), observed.len());
        assert!(expected_mz
            .iter()
            .zip(observed.iter())
            .all(|(a, b)| (a - b).abs() < 0.001));
    }

    #[test]
    fn iterate_b_ions() {
        let peptide = "PEPTIDE".parse::<Peptide>().unwrap();
        let iter = IonSeries::new(&peptide, Kind::B);

        let expected_mz = vec![
            98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
        ];

        check_within(iter, &expected_mz);
    }

    #[test]
    fn iterate_y_ions() {
        let peptide = "PEPTIDE".parse::<Peptide>().unwrap();
        let iter = IonSeries::new(&peptide, Kind::Y);

        let expected_mz = vec![
            703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
        ];

        check_within(iter, &expected_mz);
    }

    #[test]
    fn iterate_y_ions_2() {
        let peptide = "PEPTIDE".parse::<Peptide>().unwrap();
        let iter = IonSeries::new(&peptide, Kind::Y);

        let expected_mz = vec![
            352.16087, 287.639_6, 239.11319, 188.58935, 132.04732, 74.53385,
        ];

        check_within(iter, &expected_mz);
    }

    #[test]
    fn decoy() {
        let peptide = "PEPTIDE".parse::<Peptide>().unwrap();
        let iter = IonSeries::new(&peptide, Kind::Y);

        let expected_mz = vec![
            352.16087, 287.639_6, 239.11319, 188.58935, 132.04732, 74.53385,
        ];

        check_within(iter, &expected_mz);

        let peptide = "EDITPEP".parse::<Peptide>().unwrap();
        let iter = IonSeries::new(&peptide, Kind::Y);
        let expected_mz = vec![
            336.16596, 278.652_5, 222.110_46, 171.586_62, 123.060237, 58.538_94,
        ];
        check_within(iter, &expected_mz);
    }
}
