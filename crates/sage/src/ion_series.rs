use serde::Serialize;

use crate::mass::Mass;
use crate::peptide::Peptide;

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize)]
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
    cumulative_mass: f32,
    peptide: &'p Peptide,
    idx: usize,
}

impl<'p> IonSeries<'p> {
    /// Create a new [`IonSeries`] iterator for a specified peptide
    pub fn new(peptide: &'p Peptide, kind: Kind) -> Self {
        let cumulative_mass = match kind {
            Kind::B => peptide.nterm.unwrap_or_default(),
            Kind::Y => peptide.monoisotopic - peptide.nterm.unwrap_or_default(),
        };
        Self {
            kind,
            cumulative_mass,
            peptide,
            idx: 0,
        }
    }
}

impl<'p> Iterator for IonSeries<'p> {
    type Item = Ion;

    // Dynamic programming solution - memoize cumulative mass of
    // peptide fragment for fast fragment ion generation
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == self.peptide.sequence.len() - 1 {
            return None;
        }
        let r = self.peptide.sequence.get(self.idx)?;

        self.cumulative_mass += match self.kind {
            Kind::B => r.monoisotopic(),
            Kind::Y => -r.monoisotopic(),
        };
        self.idx += 1;

        Some(Ion {
            kind: self.kind,
            monoisotopic_mass: self.cumulative_mass,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{fasta::Digest, mass::PROTON, peptide::Peptide};

    fn peptide(s: &str) -> Peptide {
        Peptide::try_from(&Digest {
            sequence: s.into(),
            missed_cleavages: 0,
            decoy: false,
            start: 0,
            end: 0,
        })
        .unwrap()
    }

    fn check_within<I: Iterator<Item = Ion>>(iter: I, expected_mz: &[f32]) {
        let observed = iter.map(|ion| ion.monoisotopic_mass).collect::<Vec<f32>>();
        assert_eq!(expected_mz.len(), observed.len());
        assert!(
            expected_mz
                .iter()
                .zip(observed.iter())
                .all(|(a, b)| (a - b).abs() < 0.01),
            "{:?}",
            expected_mz
                .iter()
                .zip(observed.iter())
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>()
        );
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
        let peptide = peptide("PEPTIDE");

        // Charge state 2
        let expected_mz = vec![
            98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
        ];

        check_within(ions!(&peptide, Kind::B, 1.0), &expected_mz);
    }

    #[test]
    fn iterate_y_ions() {
        let peptide = peptide("PEPTIDE");

        // Charge state 1
        let expected_mz = vec![
            703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
        ];

        check_within(ions!(&peptide, Kind::Y, 1.0), &expected_mz);
    }

    #[test]
    fn y_index() {
        let peptide = peptide("PEPTIDE");
        // Charge state 1
        let expected_ion: Vec<(usize, f32)> = vec![
            (6, 703.31447),
            (5, 574.27188),
            (4, 477.21912),
            (3, 376.17144),
            (2, 263.08737),
            (1, 148.06043),
        ];
        assert!(IonSeries::new(&peptide, Kind::Y)
            .enumerate()
            .map(|(idx, ion)| (peptide.sequence.len().saturating_sub(1) - idx, ion))
            .zip(expected_ion.into_iter())
            .all(|((idx, ion), (idx_, mz))| {
                idx == idx_ && (ion.monoisotopic_mass + PROTON - mz).abs() <= 0.01
            }),)
    }

    #[test]
    fn index_filtering() {
        let peptide = &peptide("PEPTIDE");
        let ions = IonSeries::new(peptide, Kind::B)
            .enumerate()
            .chain(IonSeries::new(peptide, Kind::Y).enumerate())
            .filter(|(ion_idx, ion)| {
                // Don't store b1, b2, y1, y2 ions for preliminary scoring
                let ion_idx_filter = match ion.kind {
                    Kind::B => (ion_idx + 1) > 2,
                    Kind::Y => peptide.sequence.len().saturating_sub(1) - ion_idx > 2,
                };
                ion_idx_filter
            })
            .map(|(_, mut ion)| {
                ion.monoisotopic_mass += PROTON;
                ion
            })
            .collect::<Vec<_>>();

        #[rustfmt::skip]
        let expected = vec![
            Ion { kind: Kind::B, monoisotopic_mass: 324.155397 },
            Ion { kind: Kind::B, monoisotopic_mass: 425.203076 },
            Ion { kind: Kind::B, monoisotopic_mass: 538.287140 },
            Ion { kind: Kind::B, monoisotopic_mass: 653.314083 },
            Ion { kind: Kind::Y, monoisotopic_mass: 703.314477 },
            Ion { kind: Kind::Y, monoisotopic_mass: 574.271884 },
            Ion { kind: Kind::Y, monoisotopic_mass: 477.219120 },
            Ion { kind: Kind::Y, monoisotopic_mass: 376.171441 },
        ];

        assert_eq!(expected.len(), ions.len(), "{:?}\n{:?}", ions, expected);
        assert!(
            ions.iter().zip(expected.iter()).all(|(left, right)| {
                left.kind == right.kind && (left.monoisotopic_mass - right.monoisotopic_mass) <= 0.1
            }),
            "{:?}",
            ions
        );
    }

    #[test]
    fn decoy() {
        let peptide_ = peptide("PEPTIDE");

        // Charge state 2
        let expected_mz = vec![
            352.16087, 287.639_6, 239.11319, 188.58935, 132.04732, 74.53385,
        ];

        check_within(ions!(&peptide_, Kind::Y, 2.0), &expected_mz);

        let peptide = peptide("EDITPEP");

        // Charge state 2
        let expected_mz = vec![
            336.16596, 278.652_5, 222.110_46, 171.586_62, 123.060237, 58.538_94,
        ];

        check_within(ions!(&peptide, Kind::Y, 2.0), &expected_mz);
    }

    #[test]
    fn nterm_mod() {
        let mut peptide = peptide("PEPTIDE");
        peptide.static_mod('^', 229.01);

        // Charge state 1, b-ions should be TMT tagged
        let expected_b = [
            98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
        ]
        .into_iter()
        .map(|x| x + 229.01)
        .collect::<Vec<_>>();

        // y-ions shouldn't have TMT tag
        let expected_y = vec![
            703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
        ];

        check_within(ions!(&peptide, Kind::B, 1.0), &expected_b);
        check_within(ions!(&peptide, Kind::Y, 1.0), &expected_y);
    }
}
