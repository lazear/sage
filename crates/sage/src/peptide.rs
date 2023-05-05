use std::{collections::HashMap, fmt::Debug, sync::Arc};

use crate::modification::ModificationSpecificity;
use crate::{
    enzyme::{Digest, Position},
    mass::{Mass, H2O},
};
use fnv::FnvHashSet;
use itertools::Itertools;

#[derive(Clone, PartialEq, PartialOrd)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Arc<Box<[u8]>>,
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

    // pub modified_sequence: String,
    pub proteins: Vec<Arc<String>>,
}

impl Debug for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Peptide")
            .field("proteins", &self.proteins)
            .field("decoy", &self.decoy)
            .field("sequence", &self.sequence)
            .field("nterm", &self.nterm)
            .field("cterm", &self.cterm)
            .field("monoisotopic", &self.monoisotopic)
            .field("missed_cleavages", &self.missed_cleavages)
            .field("position", &self.position)
            .finish()
    }
}

impl Peptide {
    pub fn label(&self) -> i32 {
        match self.decoy {
            true => -1,
            false => 1,
        }
    }

    pub fn proteins(&self, decoy_tag: &str) -> String {
        if self.decoy {
            self.proteins
                .iter()
                .map(|s| format!("{}{}", decoy_tag, s))
                .join(";")
        } else {
            self.proteins.iter().join(";")
        }
    }

    /// Apply all variable mods in `sites` to self
    fn apply_site(&mut self, site: Site, mass: f32) {
        match site {
            Site::Nterm => {
                if self.nterm.is_none() {
                    self.nterm = Some(mass);
                    self.monoisotopic += mass;
                }
            }
            Site::Cterm => {
                if self.cterm.is_none() {
                    self.cterm = Some(mass);
                    self.monoisotopic += mass;
                }
            }
            Site::Sequence(index) => {
                if self.modifications[index as usize] == 0.0 {
                    self.modifications[index as usize] = mass;
                    self.monoisotopic += mass;
                }
            }
        }
    }

    fn push_resi(&self, acc: &mut Vec<(Site, f32)>, target: ModificationSpecificity, mass: f32) {
        match (target, self.position) {
            (ModificationSpecificity::PeptideN(None), _) => acc.push((Site::Nterm, mass)),
            (ModificationSpecificity::PeptideN(Some(resi)), _)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                acc.push((Site::Sequence(0), mass))
            }
            (ModificationSpecificity::PeptideC(None), _) => acc.push((Site::Cterm, mass)),
            (ModificationSpecificity::PeptideC(Some(resi)), _)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                acc.push((
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    mass,
                ))
            }
            (ModificationSpecificity::ProteinN(None), Position::Nterm | Position::Full) => {
                acc.push((Site::Nterm, mass))
            }
            (ModificationSpecificity::ProteinN(Some(resi)), Position::Nterm | Position::Full)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                acc.push((Site::Sequence(0), mass))
            }
            (ModificationSpecificity::ProteinC(None), Position::Cterm | Position::Full) => {
                acc.push((Site::Cterm, mass))
            }
            (ModificationSpecificity::ProteinC(Some(resi)), Position::Cterm | Position::Full)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                acc.push((
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    mass,
                ))
            }
            (ModificationSpecificity::Residue(resi), _) => {
                acc.extend(
                    self.sequence
                        .iter()
                        .enumerate()
                        .filter_map(|(idx, residue)| {
                            if resi == *residue {
                                Some((Site::Sequence(idx as u32), mass))
                            } else {
                                None
                            }
                        }),
                );
            }
            _ => {}
        }
    }

    fn static_mods(&mut self, target: ModificationSpecificity, mass: f32) {
        match (target, self.position) {
            (ModificationSpecificity::PeptideN(resi), _)
                if resi
                    .map(|r| r == *self.sequence.first().unwrap_or(&0))
                    .unwrap_or(true) =>
            {
                self.apply_site(Site::Nterm, mass)
            }
            (ModificationSpecificity::PeptideC(resi), _)
                if resi
                    .map(|r| r == *self.sequence.last().unwrap_or(&0))
                    .unwrap_or(true) =>
            {
                self.apply_site(Site::Cterm, mass)
            }
            (ModificationSpecificity::ProteinN(resi), Position::Nterm | Position::Full)
                if resi
                    .map(|r| r == *self.sequence.first().unwrap_or(&0))
                    .unwrap_or(true) =>
            {
                self.apply_site(Site::Nterm, mass)
            }
            (ModificationSpecificity::ProteinC(resi), Position::Cterm | Position::Full)
                if resi
                    .map(|r| r == *self.sequence.last().unwrap_or(&0))
                    .unwrap_or(true) =>
            {
                self.apply_site(Site::Cterm, mass)
            }
            (ModificationSpecificity::Residue(resi), _) => {
                for (idx, residue) in self.sequence.iter().enumerate() {
                    if resi == *residue && self.modifications[idx] == 0.0 {
                        self.modifications[idx] = mass;
                        self.monoisotopic += mass;
                    }
                }
            }
            _ => {}
        }
    }

    /// Apply variable modifications, then static modifications to a peptide
    pub fn apply(
        mut self,
        variable_mods: &[(ModificationSpecificity, f32)],
        static_mods: &HashMap<ModificationSpecificity, f32>,
        combinations: usize,
    ) -> Vec<Peptide> {
        if variable_mods.is_empty() {
            for (target, mass) in static_mods {
                self.static_mods(*target, *mass);
            }
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
                    for (site, mass) in combination {
                        peptide.apply_site(*site, *mass);
                    }
                    modified.push(peptide);
                }
            }

            // Apply static mods to all peptides
            for peptide in modified.iter_mut() {
                for (target, mass) in static_mods {
                    peptide.static_mods(*target, *mass);
                }
            }

            modified
        }
    }

    /// If `self` is a decoy peptide, un-reverse it
    // pub fn pseudo_forward(&self) -> Option<Peptide> {
    //     if self.decoy {
    //         let mut fwd = self.clone();
    //         if fwd.sequence.len() > 2 {
    //             let n = fwd.sequence.len().saturating_sub(1);
    //             let mut s = (*fwd.sequence).clone();
    //             s[1..n].reverse();
    //             fwd.sequence = Arc::new(s);
    //             fwd.modifications[1..n].reverse();
    //         }
    //         return Some(fwd);
    //     }
    //     None
    // }

    pub fn reverse(&self) -> Peptide {
        let mut pep = self.clone();
        pep.decoy = !self.decoy;
        let n = pep.sequence.len().saturating_sub(1);
        if n > 1 {
            let mut s = (*pep.sequence).clone();
            s[1..n].reverse();
            pep.sequence = Arc::new(s);
            pep.modifications[1..n].reverse();
        }
        pep
    }
}

fn no_duplicates(combination: &Vec<&(Site, f32)>) -> bool {
    let mut n = 0;
    let mut c = 0;
    for (site, _) in combination {
        match site {
            Site::Nterm => n += 1,
            Site::Cterm => c += 1,
            _ => {}
        }
    }

    n <= 1 && c <= 1
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Site {
    Nterm,
    Cterm,
    Sequence(u32),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PeptideError {
    InvalidSequence(String),
}

impl TryFrom<Digest> for Peptide {
    type Error = PeptideError;

    fn try_from(value: Digest) -> Result<Self, Self::Error> {
        let mut monoisotopic = H2O;
        // This is an important invariant to enforce, that ensures safety
        // while reversing peptide sequences
        if !value.sequence.is_ascii() {
            return Err(PeptideError::InvalidSequence(value.sequence));
        }
        for c in value.sequence.as_bytes() {
            let mono = c.monoisotopic();
            if mono == 0.0 {
                return Err(PeptideError::InvalidSequence(value.sequence));
            }
            monoisotopic += mono;
        }

        Ok(Peptide {
            decoy: value.decoy,
            position: value.position,
            // modified_sequence: value.sequence.clone(),
            modifications: vec![0.0; value.sequence.len()],
            sequence: Arc::new(value.sequence.into_bytes().into_boxed_slice()),
            monoisotopic,
            nterm: None,
            cterm: None,
            missed_cleavages: value.missed_cleavages,
            proteins: vec![value.protein],
        })
    }
}

impl std::fmt::Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = self.nterm {
            write!(f, "[{:+}]-", m)?;
        }
        for (c, m) in self.sequence.iter().zip(self.modifications.iter()) {
            if *m != 0.0 {
                write!(f, "{}[{:+}]", *c as char, m)?;
            } else {
                write!(f, "{}", *c as char)?;
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

    fn var_mod_sequence(
        peptide: &Peptide,
        mods: &[(ModificationSpecificity, f32)],
        combo: usize,
    ) -> Vec<String> {
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
            .digest(sequence, Default::default())
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

        use ModificationSpecificity::*;

        let mods = [
            (ProteinN(None), 42.0),
            (ProteinC(None), 11.0),
            (PeptideN(None), 12.0),
            (PeptideC(None), 19.0),
        ];
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
        use ModificationSpecificity::*;
        let variable_mods = [(Residue(b'M'), 16.0f32), (Residue(b'C'), 57.)];
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
        use ModificationSpecificity::*;
        let variable_mods = [(Residue(b'M'), 16.0f32), (Residue(b'C'), 57.)];
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
        use ModificationSpecificity::*;
        let variable_mods = [(PeptideN(None), 42.), (Residue(b'M'), 16.)];
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
        use ModificationSpecificity::*;
        let variable_mods = [(PeptideC(None), 42.), (Residue(b'M'), 16.)];
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
        use ModificationSpecificity::*;
        let variable_mods = [(Residue(b'S'), 79.), (Residue(b'S'), 541.)];
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
        for digest in trypsin.digest(fwd, Default::default()) {
            let fwd = Peptide::try_from(digest.clone()).unwrap();
            let rev = Peptide::try_from(digest.reverse()).unwrap();

            assert_eq!(fwd.decoy, false);
            assert_eq!(rev.decoy, true);
            assert!(
                fwd.sequence.len() < 4 || fwd.sequence != rev.sequence,
                "{} {}",
                fwd,
                rev
            );
            assert_eq!(rev.reverse().to_string(), fwd.to_string());
        }
    }

    #[test]
    fn apply_mods() {
        use ModificationSpecificity::*;
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
        static_mods.insert(Residue(b'C'), 57.0);

        let variable_mods = [(Residue(b'C'), 30.0)];

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
        peptide.push_resi(&mut mods, ModificationSpecificity::Residue(b'C'), 16.0);
        assert_eq!(mods, vec![(Sequence(2), 16.0), (Sequence(5), 16.0)]);
        mods.clear();

        peptide.push_resi(&mut mods, ModificationSpecificity::PeptideC(None), 16.0);
        assert_eq!(mods, vec![(Cterm, 16.0)]);
        mods.clear();

        peptide.push_resi(&mut mods, ModificationSpecificity::PeptideN(None), 16.0);
        assert_eq!(mods, vec![(Nterm, 16.0)]);
        mods.clear();

        let mut mods = vec![];
        for (residue, mass) in [("^", 12.0), ("$", 200.0), ("C", 57.0), ("A", 43.0)] {
            peptide.push_resi(&mut mods, residue.parse().unwrap(), mass);
        }

        assert_eq!(
            mods,
            vec![
                (Nterm, 12.0),
                (Cterm, 200.0),
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
