use std::{collections::HashMap, fmt::Debug, sync::Arc};

use fnv::{FnvHashSet, FnvHasher};
use itertools::Itertools;

use crate::{
    enzyme::{Digest, Position},
    mass::{Mass, H2O, VALID_AA},
};

#[derive(Clone, PartialEq, PartialOrd)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Arc<String>,
    pub modifications: Vec<f32>,
    /// Modification on peptide C-terminus
    pub nterm: Option<f32>,
    /// Modification on peptide C-terminus
    pub cterm: Option<f32>,
    /// Monoisotopic mass, inclusive of N/C-terminal mods
    pub monoisotopic: f32,
    /// Number of missed cleavages for this sequence
    pub missed_cleavages: u8,
    /// Where is this peptide located in the protein?
    pub position: Position,

    pub modified_sequence: String,
}

impl Debug for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Peptide")
            .field("decoy", &self.decoy)
            .field("sequence", &self.modified_sequence)
            .field("nterm", &self.nterm)
            .field("cterm", &self.cterm)
            .field("monoisotopic", &self.monoisotopic)
            .field("missed_cleavages", &self.missed_cleavages)
            .field("position", &self.position)
            .finish()
    }
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

    fn set_nterm_protein_mod(&mut self, m: f32) {
        if self.position == Position::Full || self.position == Position::Nterm {
            self.set_nterm_mod(m)
        }
    }

    fn set_cterm_protein_mod(&mut self, m: f32) {
        if self.position == Position::Full || self.position == Position::Cterm {
            self.set_cterm_mod(m)
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
        match residue {
            '^' => self.set_nterm_mod(mass),
            '$' => self.set_cterm_mod(mass),
            '[' => self.set_nterm_protein_mod(mass),
            ']' => self.set_cterm_protein_mod(mass),
            _ => {
                for (idx, resi) in self.sequence.as_bytes().iter().enumerate() {
                    // Don't overwrite an already modified amino acid!
                    if residue as u8 == *resi && self.modifications[idx] == 0.0 {
                        self.modifications[idx] = mass;
                        self.monoisotopic += mass;
                    }
                }
            }
        }
    }

    /// Apply all variable mods in `sites` to self
    fn apply_variable_mods(&mut self, sites: &[&(Site, f32)]) {
        for (site, mass) in sites {
            match site {
                Site::PeptideN => self.set_nterm_mod(*mass),
                Site::PeptideC => self.set_cterm_mod(*mass),
                Site::ProteinN => self.set_nterm_protein_mod(*mass),
                Site::ProteinC => self.set_cterm_protein_mod(*mass),
                Site::Sequence(index) => {
                    if self.modifications[*index as usize] == 0.0 {
                        self.modifications[*index as usize] = *mass;
                        self.monoisotopic += mass;
                    }
                }
            }
        }
    }

    fn push_resi(&self, acc: &mut Vec<(Site, f32)>, residue: char, mass: f32) {
        match (residue, self.position) {
            ('^', _) => acc.push((Site::PeptideN, mass)),
            ('$', _) => acc.push((Site::PeptideC, mass)),
            ('[', Position::Nterm | Position::Full) => acc.push((Site::ProteinN, mass)),
            (']', Position::Cterm | Position::Full) => acc.push((Site::ProteinC, mass)),
            _ => {
                acc.extend(self.sequence.as_bytes().iter().enumerate().filter_map(
                    |(idx, resi)| {
                        if *resi == residue as u8 {
                            Some((Site::Sequence(idx as u32), mass))
                        } else {
                            None
                        }
                    },
                ));
            }
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
            self.modified_sequence = self.to_string();
            vec![self]
        } else {
            let mut mods = Vec::new();
            for (residue, mass) in variable_mods.iter() {
                self.push_resi(&mut mods, *residue, *mass);
            }

            let mut modified = Vec::new();
            modified.push(self.clone());

            for n in 1..=combinations {
                'next: for combination in mods.iter().combinations(n).filter(no_duplicates) {
                    let mut set = FnvHashSet::default();
                    for (site, _) in &combination {
                        if !set.insert(*site) {
                            continue 'next;
                        }
                    }
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

                peptide.modified_sequence = peptide.to_string();
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
                let mut s = fwd.sequence.to_string();
                unsafe { s.as_bytes_mut()[1..n].reverse() }
                fwd.sequence = Arc::new(s);
                fwd.modifications[1..n].reverse();
                fwd.modified_sequence = fwd.to_string();
            }
            return Some(fwd);
        }
        None
    }

    // pub fn ambiguous(&self, other: &Peptide) -> bool {
    //     self.sequence.len() == other.sequence.len()
    //         && self.monoisotopic == other.monoisotopic
    //         && self
    //             .sequence
    //             .iter()
    //             .zip(other.sequence.iter())
    //             .all(|(l, r)| l.monoisotopic() == r.monoisotopic())
    // }
}

fn no_duplicates(combination: &Vec<&(Site, f32)>) -> bool {
    let mut n = 0;
    let mut c = 0;
    for (site, _) in combination {
        match site {
            Site::PeptideN => n += 1,
            Site::PeptideC => c += 1,
            Site::ProteinN => n += 1,
            Site::ProteinC => c += 1,
            _ => {}
        }
    }

    n <= 1 && c <= 1
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Site {
    PeptideN,
    PeptideC,
    ProteinN,
    ProteinC,
    Sequence(u32),
}

impl TryFrom<Digest> for Peptide {
    type Error = char;

    fn try_from(value: Digest) -> Result<Self, Self::Error> {
        let mut monoisotopic = H2O;

        for c in value.sequence.as_bytes() {
            let mono = c.monoisotopic();
            if mono == 0.0 {
                return Err(*c as char);
            }
            monoisotopic += mono;
        }

        Ok(Peptide {
            decoy: value.decoy,
            position: value.position,
            modified_sequence: value.sequence.clone(),
            modifications: vec![0.0; value.sequence.len()],
            sequence: Arc::new(value.sequence),
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
            write!(f, "[{:+}]-", m)?;
        }
        for (c, m) in self.sequence.chars().zip(self.modifications.iter()) {
            if *m != 0.0 {
                write!(f, "{}[{:+}]", c, m)?;
            } else {
                write!(f, "{}", c)?;
            }
        }
        if let Some(m) = self.cterm {
            write!(f, "-[{:+}]", m)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::enzyme::{Enzyme, EnzymeParameters};

    use super::*;

    fn var_mod_sequence(peptide: &Peptide, mods: &[(char, f32)], combo: usize) -> Vec<String> {
        let static_mods = HashMap::default();
        peptide
            .clone()
            .apply(&mods, &static_mods, combo)
            .into_iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
    }

    #[test]
    fn full() {
        let sequence = "MPEPTIDEKMSAGEKEND";
        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", Some('P'), true),
        };

        let peptides = tryp
            .digest(sequence)
            .into_iter()
            .map(Peptide::try_from)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(peptides.len(), 3);
        assert_eq!(peptides[0].to_string(), "MPEPTIDEK");
        assert_eq!(peptides[0].position, Position::Nterm);
        assert_eq!(peptides[1].to_string(), "MSAGEK");
        assert_eq!(peptides[1].position, Position::Internal);
        assert_eq!(peptides[2].to_string(), "END");
        assert_eq!(peptides[2].position, Position::Cterm);

        let mods = [('[', 42.0), (']', 11.0), ('^', 12.0), ('$', 19.0)];
        let a = var_mod_sequence(&peptides[0], &mods, 2);
        let b = var_mod_sequence(&peptides[1], &mods, 2);
        let c = var_mod_sequence(&peptides[2], &mods, 2);

        // Make sure no duplicates exist
        assert_eq!(
            a,
            vec![
                "MPEPTIDEK",
                "[+42]-MPEPTIDEK",
                "[+12]-MPEPTIDEK",
                "MPEPTIDEK-[+19]",
                "[+42]-MPEPTIDEK-[+19]",
                "[+12]-MPEPTIDEK-[+19]",
            ]
        );

        assert_eq!(
            b,
            vec![
                "MSAGEK",
                "[+12]-MSAGEK",
                "MSAGEK-[+19]",
                "[+12]-MSAGEK-[+19]",
            ]
        );

        assert_eq!(
            c,
            vec![
                "END",
                "END-[+11]",
                "[+12]-END",
                "END-[+19]",
                "[+12]-END-[+11]",
                "[+12]-END-[+19]",
            ]
        );
    }

    #[test]
    fn test_variable_mods() {
        let variable_mods = [('M', 16.0f32), ('C', 57.)];
        let peptide = Peptide::try_from(Digest {
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

        let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_no_effeect() {
        let variable_mods = [('M', 16.0f32), ('C', 57.)];
        let peptide = Peptide::try_from(Digest {
            sequence: "AAAAAAAA".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec!["AAAAAAAA"];
        let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_nterm() {
        let variable_mods = [('^', 42.), ('M', 16.)];
        let peptide = Peptide::try_from(Digest {
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

        let peptides = var_mod_sequence(&peptide, &variable_mods, 3);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_cterm() {
        let variable_mods = [('$', 42.), ('M', 16.)];
        let peptide = Peptide::try_from(Digest {
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

        let peptides = var_mod_sequence(&peptide, &variable_mods, 3);
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_multi() {
        let variable_mods = [('S', 79.), ('S', 541.)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GGGSGGGS".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec![
            "GGGSGGGS",
            "GGGS[+79]GGGS",
            "GGGSGGGS[+79]",
            "GGGS[+541]GGGS",
            "GGGSGGGS[+541]",
            "GGGS[+79]GGGS[+79]",
            "GGGS[+79]GGGS[+541]",
            "GGGS[+541]GGGS[+79]",
            "GGGS[+541]GGGS[+541]",
        ];

        let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
        assert_eq!(peptides, expected);
    }

    /// Check that picked-peptide approach will match forward and reverse peptides
    #[test]
    fn test_psuedo_forward() {
        let trypsin = crate::enzyme::EnzymeParameters {
            missed_cleavages: 0,
            min_len: 3,
            max_len: 30,
            enyzme: Enzyme::new("KR", Some('P'), true),
        };

        let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        for digest in trypsin.digest(fwd) {
            let mut fwd = Peptide::try_from(digest.clone()).unwrap();
            let mut rev = Peptide::try_from(digest.reverse()).unwrap();

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
        let peptide = Peptide::try_from(Digest {
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

        let peptides = peptide
            .apply(&variable_mods, &static_mods, 2)
            .into_iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>();

        assert_eq!(peptides, expected);
    }

    #[test]
    fn modification_sites() {
        use Site::*;
        let peptide = Peptide::try_from(Digest {
            sequence: "AACAACAA".into(),
            ..Default::default()
        })
        .unwrap();

        let mut mods = vec![];
        peptide.push_resi(&mut mods, 'C', 16.0);
        assert_eq!(mods, vec![(Sequence(2), 16.0), (Sequence(5), 16.0)]);
        mods.clear();

        peptide.push_resi(&mut mods, '$', 16.0);
        assert_eq!(mods, vec![(PeptideC, 16.0)]);
        mods.clear();

        peptide.push_resi(&mut mods, '^', 16.0);
        assert_eq!(mods, vec![(PeptideN, 16.0)]);
        mods.clear();

        let mut mods = vec![];
        for (residue, mass) in [('^', 12.0), ('$', 200.0), ('C', 57.0), ('A', 43.0)] {
            peptide.push_resi(&mut mods, residue, mass);
        }

        assert_eq!(
            mods,
            vec![
                (PeptideN, 12.0),
                (PeptideC, 200.0),
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
