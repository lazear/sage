use std::collections::HashMap;

use itertools::Itertools;

use crate::{
    enzyme::Digest,
    mass::{Mass, Residue, H2O, VALID_AA},
};

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Vec<Residue>,
    /// Modification on peptide C-terminus
    pub nterm: Option<f32>,
    /// Modification on peptide C-terminus
    pub cterm: Option<f32>,
    /// Monoisotopic mass, inclusive of N/C-terminal mods
    pub monoisotopic: f32,
    /// Number of missed cleavages for this sequence
    pub missed_cleavages: u8,
}

impl Peptide {
    fn set_nterm_mod(&mut self, m: f32) {
        if self.nterm.is_none() {
            self.nterm = Some(m);
            self.monoisotopic += m;
        }
    }

    fn set_cterm_mod(&mut self, m: f32) {
        if self.cterm.is_none() {
            self.cterm = Some(m);
            self.monoisotopic += m;
        }
    }

    pub fn label(&self) -> i32 {
        match self.decoy {
            true => -1,
            false => 1,
        }
    }

    /// Apply a static modification to a peptide in-place
    pub fn static_mod(&mut self, residue: char, mass: f32) {
        if residue == '^' {
            return self.set_nterm_mod(mass);
        } else if residue == '$' {
            return self.set_cterm_mod(mass);
        } else {
            for resi in self.sequence.iter_mut() {
                // Don't overwrite an already modified amino acid!
                match resi {
                    Residue::Just(c) if *c == residue => {
                        self.monoisotopic += mass;
                        *resi = Residue::Mod(residue, mass);
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply all variable mods in `sites` to self
    fn apply_variable_mods(&mut self, sites: &[&(Site, f32)]) {
        for (site, mass) in sites {
            match site {
                Site::N => self.set_nterm_mod(*mass),
                Site::C => self.set_cterm_mod(*mass),
                Site::Sequence(index) => {
                    if let Residue::Just(c) = self.sequence[*index as usize] {
                        self.sequence[*index as usize] = Residue::Mod(c, *mass);
                        self.monoisotopic += mass;
                    }
                }
            }
        }
    }

    fn modification_sites(&self, residue: char, mass: f32) -> ModificationSites {
        ModificationSites {
            peptide: self,
            index: 0,
            residue,
            mass,
        }
    }

    /// Apply variable modifications, then static modifications to a peptide
    pub fn apply(
        mut self,
        variable_mods: &[(char, f32)],
        static_mods: &HashMap<char, f32>,
        combinations: usize,
    ) -> Vec<Peptide> {
        if variable_mods.is_empty() {
            for (resi, mass) in static_mods {
                self.static_mod(*resi, *mass);
            }
            vec![self]
        } else {
            // Create list of all possible variable modifications
            let mods = variable_mods
                .iter()
                .fold(vec![], |mut acc, (residue, mass)| {
                    acc.extend(self.modification_sites(*residue, *mass));
                    acc
                });

            let mut modified = Vec::new();
            modified.push(self.clone());

            for n in 1..=combinations {
                for combination in mods.iter().combinations(n) {
                    let mut peptide = self.clone();
                    peptide.apply_variable_mods(&combination);
                    modified.push(peptide);
                }
            }

            // Apply static mods to all peptides
            for peptide in modified.iter_mut() {
                for (&residue, &mass) in static_mods {
                    peptide.static_mod(residue, mass);
                }
            }

            modified
        }
    }

    /// If `self` is a decoy peptide, un-reverse it
    pub fn pseudo_forward(&self) -> Option<Peptide> {
        if self.decoy {
            let mut fwd = self.clone();
            if fwd.sequence.len() > 2 {
                let n = fwd.sequence.len().saturating_sub(1);
                fwd.sequence[1..n].reverse();
            }
            return Some(fwd);
        }
        None
    }

    pub fn ambiguous(&self, other: &Peptide) -> bool {
        self.sequence.len() == other.sequence.len()
            && self.monoisotopic == other.monoisotopic
            && self
                .sequence
                .iter()
                .zip(other.sequence.iter())
                .all(|(l, r)| l.monoisotopic() == r.monoisotopic())
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Site {
    N,
    C,
    Sequence(u32),
}

struct ModificationSites<'a> {
    peptide: &'a Peptide,
    index: usize,
    residue: char,
    mass: f32,
}

impl<'a> Iterator for ModificationSites<'a> {
    type Item = (Site, f32);

    fn next(&mut self) -> Option<Self::Item> {
        match (self.residue, self.index) {
            ('^', 0) => {
                self.index = self.peptide.sequence.len();
                return Some((Site::N, self.mass));
            }
            ('$', 0) => {
                self.index = self.peptide.sequence.len();
                return Some((Site::C, self.mass));
            }
            _ => {}
        };
        while self.index < self.peptide.sequence.len() {
            let idx = self.index;
            self.index += 1;
            match self.peptide.sequence[idx] {
                Residue::Just(r) if r == self.residue => {
                    return Some((Site::Sequence(idx as u32), self.mass))
                }
                _ => continue,
            }
        }
        None
    }
}

impl TryFrom<&Digest> for Peptide {
    type Error = char;

    fn try_from(value: &Digest) -> Result<Self, Self::Error> {
        let mut sequence = Vec::with_capacity(value.sequence.len());
        let mut monoisotopic = H2O;

        for c in value.sequence.chars() {
            if !VALID_AA.contains(&c) {
                return Err(c);
            }
            monoisotopic += c.monoisotopic();
            sequence.push(Residue::Just(c));
        }

        Ok(Peptide {
            decoy: value.decoy,
            sequence,
            monoisotopic,
            nterm: None,
            cterm: None,
            missed_cleavages: value.missed_cleavages,
        })
    }
}

impl std::fmt::Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = self.nterm {
            if m.is_sign_positive() {
                write!(f, "[+{}]-", m)?;
            } else {
                write!(f, "[{}]-", m)?;
            }
        }
        f.write_str(
            &self
                .sequence
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(""),
        )?;
        if let Some(m) = self.cterm {
            if m.is_sign_positive() {
                write!(f, "-[+{}]", m)?;
            } else {
                write!(f, "-[{}]", m)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::enzyme::Enzyme;

    use super::*;

    #[test]
    fn test_variable_mods() {
        let variable_mods = [('M', 16.0f32), ('C', 57.)];
        let peptide = Peptide::try_from(&Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GCMGCMG",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "GC[+57]MGCMG",
            "GCMGC[+57]MG",
            "GCM[+16]GCM[+16]G",
            "GC[+57]M[+16]GCMG",
            "GCM[+16]GC[+57]MG",
            "GC[+57]MGCM[+16]G",
            "GCMGC[+57]M[+16]G",
            "GC[+57]MGC[+57]MG",
        ];

        let peptides = peptide
            .apply(&variable_mods, &HashMap::default(), 2)
            .into_iter()
            .map(|peptide| peptide.to_string())
            .collect::<Vec<_>>();
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_no_effeect() {
        let variable_mods = [('M', 16.0f32), ('C', 57.)];
        let peptide = Peptide::try_from(&Digest {
            sequence: "AAAAAAAA".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec!["AAAAAAAA"];

        let peptides = peptide
            .apply(&variable_mods, &HashMap::default(), 2)
            .into_iter()
            .map(|peptide| peptide.to_string())
            .collect::<Vec<_>>();
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_nterm() {
        let variable_mods = [('^', 42.), ('M', 16.)];
        let peptide = Peptide::try_from(&Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GCMGCMG",
            "[+42]-GCMGCMG",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "[+42]-GCM[+16]GCMG",
            "[+42]-GCMGCM[+16]G",
            "GCM[+16]GCM[+16]G",
            "[+42]-GCM[+16]GCM[+16]G",
        ];

        let peptides = peptide
            .apply(&variable_mods, &HashMap::default(), 3)
            .into_iter()
            .map(|peptide| peptide.to_string())
            .collect::<Vec<_>>();
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_cterm() {
        let variable_mods = [('$', 42.), ('M', 16.)];
        let peptide = Peptide::try_from(&Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GCMGCMG",
            "GCMGCMG-[+42]",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "GCM[+16]GCMG-[+42]",
            "GCMGCM[+16]G-[+42]",
            "GCM[+16]GCM[+16]G",
            "GCM[+16]GCM[+16]G-[+42]",
        ];

        let peptides = peptide
            .apply(&variable_mods, &HashMap::default(), 3)
            .into_iter()
            .map(|peptide| peptide.to_string())
            .collect::<Vec<_>>();
        assert_eq!(peptides, expected);
    }

    /// Check that picked-peptide approach will match forward and reverse peptides
    #[test]
    fn test_psuedo_forward() {
        let trypsin = crate::enzyme::EnzymeParameters {
            missed_cleavages: 0,
            min_len: 3,
            max_len: 30,
            enyzme: Enzyme::new("KR", Some('P')),
        };

        let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        for digest in trypsin.digest(fwd) {
            let mut fwd = Peptide::try_from(&digest).unwrap();
            let mut rev = Peptide::try_from(&digest.reverse()).unwrap();

            assert_eq!(fwd.decoy, false);
            assert_eq!(rev.decoy, true);
            assert!(
                fwd.sequence.len() < 4 || fwd.sequence != rev.sequence,
                "{} {}",
                fwd,
                rev
            );
            assert_eq!(fwd.pseudo_forward(), None);
            assert_eq!(rev.pseudo_forward().unwrap().to_string(), fwd.to_string());

            fwd.static_mod('E', 15.0);
            rev.static_mod('E', 15.0);
            assert_eq!(rev.pseudo_forward().unwrap().to_string(), fwd.to_string());
        }
    }

    #[test]
    fn apply_mods() {
        let peptide = Peptide::try_from(&Digest {
            sequence: "AACAACAA".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "AAC[+57]AAC[+57]AA",
            "AAC[+30]AAC[+57]AA",
            "AAC[+57]AAC[+30]AA",
            "AAC[+30]AAC[+30]AA",
        ];

        let mut static_mods = HashMap::new();
        static_mods.insert('C', 57.0);

        let variable_mods = [('C', 30.0)];

        let actual = peptide
            .apply(&variable_mods, &static_mods, 2)
            .into_iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
    }

    fn mod_sites(peptide: &Peptide, residue: char, mass: f32) -> ModificationSites {
        ModificationSites {
            peptide,
            index: 0,
            residue,
            mass,
        }
    }

    #[test]
    fn modification_sites() {
        use Site::*;
        let peptide = Peptide::try_from(&Digest {
            sequence: "AACAACAA".into(),
            ..Default::default()
        })
        .unwrap();

        let mods = mod_sites(&peptide, 'C', 16.0);
        assert_eq!(
            mods.collect::<Vec<_>>(),
            vec![(Sequence(2), 16.0), (Sequence(5), 16.0)]
        );

        let mods = mod_sites(&peptide, '$', 16.0);
        assert_eq!(mods.collect::<Vec<_>>(), vec![(C, 16.0)]);

        let mods = mod_sites(&peptide, '^', 16.0);
        assert_eq!(mods.collect::<Vec<_>>(), vec![(N, 16.0)]);

        let mods = [('^', 12.0), ('$', 200.0), ('C', 57.0), ('A', 43.0)];
        let mods = mods.iter().fold(vec![], |mut acc, m| {
            acc.extend(mod_sites(&peptide, m.0, m.1));
            acc
        });

        assert_eq!(
            mods,
            vec![
                (N, 12.0),
                (C, 200.0),
                (Sequence(2), 57.0),
                (Sequence(5), 57.0),
                (Sequence(0), 43.0),
                (Sequence(1), 43.0),
                (Sequence(3), 43.0),
                (Sequence(4), 43.0),
                (Sequence(6), 43.0),
                (Sequence(7), 43.0),
            ]
        );
    }
}
