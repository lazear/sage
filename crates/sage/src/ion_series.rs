use serde::{Deserialize, Serialize};

use crate::mass::monoisotopic;
use crate::modification::Modifications;
use crate::peptide::Peptide;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    A,
    B,
    C,
    X,
    Y,
    Z,
    /// Radical z ion (y - NH2).
    #[serde(rename = "z_radical")]
    ZRadical,
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
    mods: &'p Modifications,
    idx: usize,
}

impl<'p> IonSeries<'p> {
    /// Create a new [`IonSeries`] iterator for a specified peptide. `mods` must
    /// be the registry of the database the peptide belongs to.
    pub fn new(peptide: &'p Peptide, kind: Kind, mods: &'p Modifications) -> Self {
        const C: f32 = 12.0;
        const O: f32 = 15.994914;
        const H: f32 = 1.007825;
        const N: f32 = 14.003074;
        const NH2: f32 = N + H * 2.0;
        const NH3: f32 = N + H * 3.0;

        let nterm = mods.mass(peptide.nterm);
        let cumulative_mass = match kind {
            Kind::A => nterm - (C + O),
            Kind::B => nterm,
            Kind::C => nterm + NH3,
            Kind::X => peptide.monoisotopic - nterm + (C + O - NH3 + N + H),
            Kind::Y => peptide.monoisotopic - nterm,
            Kind::Z => peptide.monoisotopic - nterm - NH3,
            Kind::ZRadical => peptide.monoisotopic - nterm - NH2,
        };
        Self {
            kind,
            cumulative_mass,
            peptide,
            mods,
            idx: 0,
        }
    }
}

impl<'p> Iterator for IonSeries<'p> {
    type Item = Ion;

    // Dynamic programming solution - memoize cumulative mass of
    // peptide fragment for fast fragment ion generation
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.peptide.sequence.len() - 1 {
            return None;
        }
        let r = self.peptide.sequence[self.idx];
        let m = self.mods.mass(self.peptide.modifications[self.idx]);

        self.cumulative_mass += match self.kind {
            Kind::A | Kind::B | Kind::C => monoisotopic(r) + m,
            Kind::X | Kind::Y | Kind::Z | Kind::ZRadical => -(monoisotopic(r) + m),
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
    use crate::modification::ModificationSpecificity;
    use crate::peptide::Peptide;
    use crate::{enzyme::Digest, mass::PROTON};

    fn peptide(s: &str) -> Peptide {
        Peptide::try_from(Digest {
            sequence: s.into(),
            ..Default::default()
        })
        .unwrap()
    }

    /// Apply a single static modification to `peptide`
    fn with_static_mod(
        peptide: Peptide,
        target: ModificationSpecificity,
        mass: f32,
    ) -> (Peptide, Modifications) {
        let mut mods = Modifications::default();
        let id = mods.insert(mass, None).unwrap();
        let peptide = peptide
            .apply(&[], &[(target, id)], &mods, 1, None)
            .remove(0);
        (peptide, mods)
    }

    fn check_within<I: Iterator<Item = Ion>>(iter: I, expected_mz: &[f32]) {
        let observed = iter.map(|ion| ion.monoisotopic_mass).collect::<Vec<f32>>();
        assert_eq!(expected_mz.len(), observed.len());
        assert!(
            expected_mz
                .iter()
                .zip(observed.iter())
                .all(|(a, b)| (a - b).abs() < 0.005),
            "{:?}",
            expected_mz
                .iter()
                .zip(observed.iter())
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>()
        );
    }

    macro_rules! ions {
        ($peptide:expr, $kind:expr, $charge:expr) => {
            ions!($peptide, $kind, $charge, &Modifications::default())
        };
        ($peptide:expr, $kind:expr, $charge:expr, $mods:expr) => {{
            IonSeries::new($peptide, $kind, $mods).map(|mut ion| {
                ion.monoisotopic_mass = (ion.monoisotopic_mass + $charge * PROTON) / $charge;
                ion
            })
        }};
    }

    #[test]
    fn abc_xyz() {
        let peptide = peptide("PEPTIDE");
        let expected_a = vec![70.065, 199.108, 296.160, 397.208, 510.292, 625.32];
        let expected_b = vec![98.0600, 227.1026, 324.155, 425.2030, 538.287, 653.314];
        let expected_c = vec![115.086, 244.129, 341.182, 442.229, 555.314, 670.341];
        let expected_x = vec![729.294, 600.251, 503.198, 402.151, 289.066, 174.039];
        let expected_y = vec![703.314, 574.2719, 477.219, 376.171, 263.0874, 148.0604];
        let expected_z = vec![686.288, 557.245, 460.193, 359.145, 246.061, 131.034];

        check_within(ions!(&peptide, Kind::A, 1.0), &expected_a);
        check_within(ions!(&peptide, Kind::B, 1.0), &expected_b);
        check_within(ions!(&peptide, Kind::C, 1.0), &expected_c);
        check_within(ions!(&peptide, Kind::X, 1.0), &expected_x);
        check_within(ions!(&peptide, Kind::Y, 1.0), &expected_y);
        check_within(ions!(&peptide, Kind::Z, 1.0), &expected_z);
    }

    #[test]
    fn c_ion_mass() {
        let peptide = peptide("PEPTIDE");
        let observed = ions!(&peptide, Kind::C, 1.0)
            .next()
            .unwrap()
            .monoisotopic_mass;
        // Detect the 0.00055 Da ammonia error.
        assert!((observed - 115.08659).abs() < 0.00003, "{observed}");
    }

    #[test]
    fn radical_z_config() {
        let json = r#"["c","z_radical","z"]"#;
        let kinds: Vec<Kind> = serde_json::from_str(json).unwrap();
        assert_eq!(kinds, vec![Kind::C, Kind::ZRadical, Kind::Z]);
        assert_eq!(serde_json::to_string(&kinds).unwrap(), json);
    }

    #[test]
    fn iterate_z_radical_ions() {
        let peptide = peptide("PEPTIDE");
        let expected_mz = vec![
            687.29575, 558.25316, 461.20040, 360.15272, 247.06865, 132.04171,
        ];
        check_within(ions!(&peptide, Kind::ZRadical, 1.0), &expected_mz);

        let expected_mz = vec![
            344.15152, 279.63022, 231.10384, 180.58000, 124.03797, 66.52449,
        ];
        check_within(ions!(&peptide, Kind::ZRadical, 2.0), &expected_mz);
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
        assert!(IonSeries::new(&peptide, Kind::Y, &Modifications::default())
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
        let mods = Modifications::default();
        let ions = IonSeries::new(peptide, Kind::B, &mods)
            .enumerate()
            .chain(IonSeries::new(peptide, Kind::Y, &mods).enumerate())
            .filter(|(ion_idx, ion)| {
                // Don't store b1, b2, y1, y2 ions for preliminary scoring
                let ion_idx_filter = match ion.kind {
                    Kind::A | Kind::B | Kind::C => (ion_idx + 1) > 2,
                    Kind::X | Kind::Y | Kind::Z | Kind::ZRadical => {
                        peptide.sequence.len().saturating_sub(1) - ion_idx > 2
                    }
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
        let (peptide, mods) = with_static_mod(
            peptide("PEPTIDE"),
            ModificationSpecificity::PeptideN(None),
            229.01,
        );

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

        let expected_z_radical = expected_y
            .iter()
            .map(|mz| mz - 16.018724)
            .collect::<Vec<_>>();

        check_within(ions!(&peptide, Kind::B, 1.0, &mods), &expected_b);
        check_within(ions!(&peptide, Kind::Y, 1.0, &mods), &expected_y);
        check_within(
            ions!(&peptide, Kind::ZRadical, 1.0, &mods),
            &expected_z_radical,
        );
    }

    #[test]
    fn cterm_mod() {
        let (peptide, mods) = with_static_mod(
            peptide("PEPTIDE"),
            ModificationSpecificity::PeptideC(None),
            229.01,
        );
        assert!((peptide.monoisotopic - 1028.37).abs() < 0.001);

        // b-ions should not be tagged
        let expected_b = [
            98.06004, 227.10263, 324.155_4, 425.203_06, 538.287_2, 653.314_1,
        ];

        // y-ions should be tagged
        let expected_y = vec![
            703.31447, 574.27188, 477.21912, 376.17144, 263.08737, 148.06043,
        ]
        .into_iter()
        .map(|x| x + 229.01)
        .collect::<Vec<_>>();

        let expected_z_radical = expected_y
            .iter()
            .map(|mz| mz - 16.018724)
            .collect::<Vec<_>>();

        check_within(ions!(&peptide, Kind::B, 1.0, &mods), &expected_b);
        check_within(ions!(&peptide, Kind::Y, 1.0, &mods), &expected_y);
        check_within(
            ions!(&peptide, Kind::ZRadical, 1.0, &mods),
            &expected_z_radical,
        );
    }

    #[test]
    fn internal_mod() {
        let (peptide, mods) = with_static_mod(
            peptide("PEPTIDE"),
            ModificationSpecificity::Residue(b'I'),
            29.0,
        );

        let expected_b = [
            98.06004,
            227.10263,
            324.155_4,
            425.203_06,
            538.287_2 + 29.0,
            653.314_1 + 29.0,
        ];

        let expected_y = vec![
            703.31447 + 29.0,
            574.27188 + 29.0,
            477.21912 + 29.0,
            376.17144 + 29.0,
            263.08737,
            148.06043,
        ];

        let expected_z_radical = expected_y
            .iter()
            .map(|mz| mz - 16.018724)
            .collect::<Vec<_>>();

        check_within(ions!(&peptide, Kind::B, 1.0, &mods), &expected_b);
        check_within(ions!(&peptide, Kind::Y, 1.0, &mods), &expected_y);
        check_within(
            ions!(&peptide, Kind::ZRadical, 1.0, &mods),
            &expected_z_radical,
        );
    }
}
