use std::{fmt::Debug, sync::Arc};

use crate::modification::{ModId, ModificationSpecificity, Modifications};
use crate::{
    enzyme::{Digest, DigestGroup, Position},
    mass::{monoisotopic, H2O},
};
use fnv::FnvHashSet;
use itertools::Itertools;

#[derive(Clone, PartialEq, Default)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Arc<[u8]>,
    /// Modification at each residue, as ids into the database's [`Modifications`]
    pub modifications: Box<[Option<ModId>]>,
    /// Modification on peptide N-terminus
    pub nterm: Option<ModId>,
    /// Modification on peptide C-terminus
    pub cterm: Option<ModId>,
    /// Monoisotopic mass, inclusive of N/C-terminal mods
    pub monoisotopic: f32,
    /// Number of missed cleavages for this sequence
    pub missed_cleavages: u8,
    /// Is this a semi-enzymatic peptide?
    pub semi_enzymatic: bool,
    /// Where is this peptide located in the protein?
    pub position: Position,

    pub proteins: Vec<Arc<str>>,
}

impl Peptide {
    pub fn initial_sort(&self, other: &Self) -> std::cmp::Ordering {
        self.sequence
            .cmp(&other.sequence)
            .then_with(|| self.modifications.cmp(&other.modifications))
            .then_with(|| self.nterm.cmp(&other.nterm))
            .then_with(|| self.cterm.cmp(&other.cterm))
    }
}

impl Debug for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Peptide")
            .field("proteins", &self.proteins)
            .field("decoy", &self.decoy)
            .field(
                "sequence",
                &std::str::from_utf8(&self.sequence).unwrap_or("error"),
            )
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

    pub fn proteins(&self, decoy_tag: &str, generate_decoys: bool) -> String {
        if self.decoy {
            self.proteins
                .iter()
                .map(|s| {
                    if generate_decoys {
                        format!("{}{}", decoy_tag, s)
                    } else {
                        s.to_string()
                    }
                })
                .join(";")
        } else {
            self.proteins.iter().join(";")
        }
    }

    fn modification_mass(&self, mods: &Modifications) -> f32 {
        self.modifications
            .iter()
            .map(|&m| mods.mass(m))
            .sum::<f32>()
            + mods.mass(self.nterm)
            + mods.mass(self.cterm)
    }

    /// Apply a modification to `site`, unless it is already modified
    fn apply_site(&mut self, site: Site, id: ModId) {
        let slot = match site {
            Site::Nterm => &mut self.nterm,
            Site::Cterm => &mut self.cterm,
            Site::Sequence(index) => &mut self.modifications[index as usize],
        };
        slot.get_or_insert(id);
    }

    /// Enumerate all possible variable modification sites
    fn push_resi(
        &self,
        acc: &mut Vec<(Site, ModId, usize)>,
        target: ModificationSpecificity,
        id: ModId,
        mod_idx: usize,
    ) {
        match (target, self.position) {
            (ModificationSpecificity::PeptideN(None), _) => acc.push((Site::Nterm, id, mod_idx)),
            (ModificationSpecificity::PeptideN(Some(resi)), _)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                acc.push((Site::Sequence(0), id, mod_idx))
            }
            (ModificationSpecificity::PeptideC(None), _) => acc.push((Site::Cterm, id, mod_idx)),
            (ModificationSpecificity::PeptideC(Some(resi)), _)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                acc.push((
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    id,
                    mod_idx,
                ))
            }
            (ModificationSpecificity::ProteinN(None), Position::Nterm | Position::Full) => {
                acc.push((Site::Nterm, id, mod_idx))
            }
            (ModificationSpecificity::ProteinN(Some(resi)), Position::Nterm | Position::Full)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                acc.push((Site::Sequence(0), id, mod_idx))
            }
            (ModificationSpecificity::ProteinC(None), Position::Cterm | Position::Full) => {
                acc.push((Site::Cterm, id, mod_idx))
            }
            (ModificationSpecificity::ProteinC(Some(resi)), Position::Cterm | Position::Full)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                acc.push((
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    id,
                    mod_idx,
                ))
            }
            (ModificationSpecificity::Residue(resi), _) => {
                acc.extend(
                    self.sequence
                        .iter()
                        .enumerate()
                        .filter_map(|(idx, residue)| {
                            if resi == *residue {
                                Some((Site::Sequence(idx as u32), id, mod_idx))
                            } else {
                                None
                            }
                        }),
                );
            }
            _ => {}
        }
    }

    fn static_mods(&mut self, target: ModificationSpecificity, id: ModId) {
        match (target, self.position) {
            (ModificationSpecificity::PeptideN(None), _) => self.apply_site(Site::Nterm, id),
            (ModificationSpecificity::PeptideN(Some(resi)), _)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                self.apply_site(Site::Sequence(0), id)
            }
            (ModificationSpecificity::PeptideC(None), _) => self.apply_site(Site::Cterm, id),
            (ModificationSpecificity::PeptideC(Some(resi)), _)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                self.apply_site(
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    id,
                )
            }
            (ModificationSpecificity::ProteinN(None), Position::Nterm | Position::Full) => {
                self.apply_site(Site::Nterm, id)
            }
            (ModificationSpecificity::ProteinN(Some(resi)), Position::Nterm | Position::Full)
                if resi == *self.sequence.first().unwrap_or(&0) =>
            {
                self.apply_site(Site::Sequence(0), id)
            }
            (ModificationSpecificity::ProteinC(None), Position::Cterm | Position::Full) => {
                self.apply_site(Site::Cterm, id)
            }
            (ModificationSpecificity::ProteinC(Some(resi)), Position::Cterm | Position::Full)
                if resi == *self.sequence.last().unwrap_or(&0) =>
            {
                self.apply_site(
                    Site::Sequence(self.sequence.len().saturating_sub(1) as u32),
                    id,
                )
            }
            (ModificationSpecificity::Residue(resi), _) => {
                for (idx, residue) in self.sequence.iter().enumerate() {
                    if resi == *residue {
                        self.modifications[idx].get_or_insert(id);
                    }
                }
            }
            _ => {}
        }
    }

    /// Apply variable modifications, then static modifications to a peptide.
    /// `variable_mods` entries are `(specificity, modification, per_mod_limit)`,
    /// and `mods` is the registry those ids (and `static_mods`) refer to.
    /// `max_combinations` caps total variants (unmodified + modified); fewer-PTM
    /// variants are always generated first so they are preferred when truncating.
    pub fn apply(
        mut self,
        variable_mods: &[(ModificationSpecificity, ModId, Option<usize>)],
        static_mods: &[(ModificationSpecificity, ModId)],
        mods: &Modifications,
        combinations: usize,
        max_combinations: Option<usize>,
    ) -> Vec<Peptide> {
        if variable_mods.is_empty() {
            for (target, id) in static_mods {
                self.static_mods(*target, *id);
            }
            self.monoisotopic += self.modification_mass(mods);
            vec![self]
        } else {
            let mut sites: Vec<(Site, ModId, usize)> = Vec::new();
            for (mod_idx, (residue, id, _limit)) in variable_mods.iter().enumerate() {
                self.push_resi(&mut sites, *residue, *id, mod_idx);
            }

            let mut modified = Vec::new();
            modified.push(self.clone());
            let mut mod_counts = vec![0usize; variable_mods.len()];
            let mut touched_mods = Vec::new();

            'outer: for n in 1..=combinations {
                'next: for combination in sites.iter().combinations(n).filter(no_duplicates) {
                    let mut set = FnvHashSet::default();
                    for (site, _, _) in &combination {
                        if !set.insert(*site) {
                            continue 'next;
                        }
                    }

                    // Enforce per-mod count limits
                    for (_, _, mod_idx) in &combination {
                        if mod_counts[*mod_idx] == 0 {
                            touched_mods.push(*mod_idx);
                        }
                        mod_counts[*mod_idx] += 1;
                    }
                    let exceeds_limit = touched_mods.iter().any(|&mod_idx| {
                        variable_mods[mod_idx]
                            .2
                            .map_or(false, |limit| mod_counts[mod_idx] > limit)
                    });
                    for mod_idx in touched_mods.drain(..) {
                        mod_counts[mod_idx] = 0;
                    }
                    if exceeds_limit {
                        continue 'next;
                    }

                    if max_combinations.map_or(false, |cap| modified.len() >= cap) {
                        break 'outer;
                    }

                    let mut peptide = self.clone();
                    for (site, id, _) in combination {
                        peptide.apply_site(*site, *id);
                    }
                    modified.push(peptide);
                }
            }

            // Apply static mods to all peptides
            for peptide in modified.iter_mut() {
                for (target, id) in static_mods {
                    peptide.static_mods(*target, *id);
                }
                peptide.monoisotopic += peptide.modification_mass(mods);
            }

            modified
        }
    }

    pub fn reverse(&self) -> Peptide {
        let mut pep = self.clone();
        pep.decoy = !self.decoy;
        let n = pep.sequence.len().saturating_sub(1);
        if n > 1 {
            let mut s = Vec::from(pep.sequence.as_ref());
            s[1..n].reverse();
            pep.sequence = Arc::from(s.into_boxed_slice());
            pep.modifications[1..n].reverse();
        }
        pep
    }
}

fn no_duplicates(combination: &Vec<&(Site, ModId, usize)>) -> bool {
    let mut n = 0;
    let mut c = 0;
    for (site, _, _) in combination {
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

impl TryFrom<DigestGroup> for Peptide {
    type Error = PeptideError;

    fn try_from(value: DigestGroup) -> Result<Self, Self::Error> {
        let mut pep = Peptide::try_from(value.reference)?;
        pep.proteins = value.proteins;
        Ok(pep)
    }
}

impl TryFrom<Digest> for Peptide {
    type Error = PeptideError;

    fn try_from(value: Digest) -> Result<Self, Self::Error> {
        let mut mass = H2O;
        // This is an important invariant to enforce, that ensures safety
        // while reversing peptide sequences
        if !value.sequence.is_ascii() {
            return Err(PeptideError::InvalidSequence(value.sequence));
        }
        for c in value.sequence.as_bytes() {
            let mono = monoisotopic(*c);
            if mono == 0.0 {
                return Err(PeptideError::InvalidSequence(value.sequence));
            }
            mass += mono;
        }

        Ok(Peptide {
            decoy: value.decoy,
            position: value.position,
            modifications: vec![None; value.sequence.len()].into_boxed_slice(),
            sequence: Arc::from(value.sequence.into_bytes().into_boxed_slice()),
            monoisotopic: mass,
            nterm: None,
            cterm: None,
            missed_cleavages: value.missed_cleavages,
            semi_enzymatic: value.semi_enzymatic,
            proteins: vec![value.protein],
        })
    }
}

impl Peptide {
    /// ProForma representation of this peptide. `mods` must be the registry of
    /// the database this peptide belongs to.
    pub fn display<'a>(&'a self, mods: &'a Modifications) -> Display<'a> {
        Display {
            peptide: self,
            mods,
        }
    }
}

/// See [`Peptide::display`]
pub struct Display<'a> {
    peptide: &'a Peptide,
    mods: &'a Modifications,
}

impl std::fmt::Display for Display<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let peptide = self.peptide;
        if let Some(id) = peptide.nterm {
            write!(f, "[{}]-", self.mods.get(id))?;
        }
        for (c, m) in peptide.sequence.iter().zip(peptide.modifications.iter()) {
            match m {
                Some(id) => write!(f, "{}[{}]", *c as char, self.mods.get(*id))?,
                None => write!(f, "{}", *c as char)?,
            }
        }
        if let Some(id) = peptide.cterm {
            write!(f, "-[{}]", self.mods.get(id))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::enzyme::{Enzyme, EnzymeParameters};

    use super::*;

    /// Register unnamed, mass-only modifications and apply them to `peptide`
    fn apply_masses(
        peptide: &Peptide,
        variable_mods: &[(ModificationSpecificity, f32, Option<usize>)],
        static_mods: &[(ModificationSpecificity, f32)],
        combinations: usize,
        max_combinations: Option<usize>,
    ) -> Vec<String> {
        let mut mods = Modifications::default();
        let variable_mods = variable_mods
            .iter()
            .map(|&(s, m, limit)| (s, mods.insert(m, None).unwrap(), limit))
            .collect::<Vec<_>>();
        let static_mods = static_mods
            .iter()
            .map(|&(s, m)| (s, mods.insert(m, None).unwrap()))
            .collect::<Vec<_>>();
        peptide
            .clone()
            .apply(
                &variable_mods,
                &static_mods,
                &mods,
                combinations,
                max_combinations,
            )
            .iter()
            .map(|p| p.display(&mods).to_string())
            .collect()
    }

    fn var_mod_sequence(
        peptide: &Peptide,
        mods: &[(ModificationSpecificity, f32)],
        combo: usize,
    ) -> Vec<String> {
        let mods_with_limits: Vec<(ModificationSpecificity, f32, Option<usize>)> =
            mods.iter().map(|&(s, m)| (s, m, None)).collect();
        apply_masses(peptide, &mods_with_limits, &[], combo, None)
    }

    #[test]
    fn full() {
        let sequence = "MPEPTIDEKMSAGEKEND";
        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: 50,
            missed_cleavages: 0,
            enzyme: Enzyme::new("KR", "P", true, false),
        };

        let peptides = tryp
            .digest(sequence, Default::default())
            .into_iter()
            .map(Peptide::try_from)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(peptides.len(), 3);
        assert_eq!(
            peptides[0].display(&Modifications::default()).to_string(),
            "MPEPTIDEK"
        );
        assert_eq!(peptides[0].position, Position::Nterm);
        assert_eq!(
            peptides[1].display(&Modifications::default()).to_string(),
            "MSAGEK"
        );
        assert_eq!(peptides[1].position, Position::Internal);
        assert_eq!(
            peptides[2].display(&Modifications::default()).to_string(),
            "END"
        );
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
            enzyme: Enzyme::new("KR", "P", true, false),
        };

        let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        for digest in trypsin.digest(fwd, Default::default()) {
            let fwd = Peptide::try_from(digest.clone()).unwrap();
            let rev = Peptide::try_from(digest.reverse()).unwrap();

            assert_eq!(fwd.decoy, false);
            assert_eq!(rev.decoy, true);
            let mods = Modifications::default();
            assert!(
                fwd.sequence.len() < 4 || fwd.sequence != rev.sequence,
                "{} {}",
                fwd.display(&mods),
                rev.display(&mods)
            );
            assert_eq!(
                rev.reverse().display(&mods).to_string(),
                fwd.display(&mods).to_string()
            );
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

        let static_mods = [(Residue(b'C'), 57.0)];
        let variable_mods = [(Residue(b'C'), 30.0, None)];

        let peptides = apply_masses(&peptide, &variable_mods, &static_mods, 2, None);

        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_per_mod_limit() {
        use ModificationSpecificity::*;
        // GCMGCMG has two M residues; limit oxidation to max 1 per peptide
        let variable_mods = [(Residue(b'M'), 16.0f32, Some(1))];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, None);

        // Should get unmodified + each single-M variant, but NOT the double-M variant
        let expected = vec!["GCMGCMG", "GCM[+16]GCMG", "GCMGCM[+16]G"];
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_max_combinations() {
        use ModificationSpecificity::*;
        // GCMGCMG with oxidation and carbamidomethylation would normally yield many variants;
        // cap at 4 total (unmodified + 3 modified)
        let variable_mods = [(Residue(b'M'), 16.0f32, None), (Residue(b'C'), 57.0, None)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, Some(4));

        // Cap at 4: unmodified + the first 3 single-mod variants (fewest PTMs first)
        assert_eq!(peptides.len(), 4);
        assert_eq!(peptides[0], "GCMGCMG");
    }

    #[test]
    fn modification_sites() {
        use Site::*;
        let peptide = Peptide::try_from(Digest {
            sequence: "AACAACAA".into(),
            ..Default::default()
        })
        .unwrap();

        let mut registry = Modifications::default();
        let ox = registry.insert(16.0, None).unwrap();

        let mut mods = vec![];
        peptide.push_resi(&mut mods, ModificationSpecificity::Residue(b'C'), ox, 0);
        assert_eq!(mods, vec![(Sequence(2), ox, 0), (Sequence(5), ox, 0)]);
        mods.clear();

        peptide.push_resi(&mut mods, ModificationSpecificity::PeptideC(None), ox, 0);
        assert_eq!(mods, vec![(Cterm, ox, 0)]);
        mods.clear();

        peptide.push_resi(&mut mods, ModificationSpecificity::PeptideN(None), ox, 0);
        assert_eq!(mods, vec![(Nterm, ox, 0)]);
        mods.clear();

        let ids = [12.0, 200.0, 57.0, 43.0].map(|mass| registry.insert(mass, None).unwrap());
        let mut mods = vec![];
        for (idx, residue) in ["^", "$", "C", "A"].iter().enumerate() {
            peptide.push_resi(&mut mods, residue.parse().unwrap(), ids[idx], idx);
        }

        let [n, c, cam, a] = ids;
        assert_eq!(
            mods,
            vec![
                (Nterm, n, 0),
                (Cterm, c, 1),
                (Sequence(2), cam, 2),
                (Sequence(5), cam, 2),
                (Sequence(0), a, 3),
                (Sequence(1), a, 3),
                (Sequence(3), a, 3),
                (Sequence(4), a, 3),
                (Sequence(6), a, 3),
                (Sequence(7), a, 3),
            ]
        );
    }

    #[test]
    fn test_per_mod_limit_exactly_met() {
        use ModificationSpecificity::*;
        // Limit of 2 on a peptide with exactly 2 M residues — all combos should be allowed
        let variable_mods = [(Residue(b'M'), 16.0f32, Some(2))];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, None);

        // No restriction: unmodified + each single + double
        let expected = vec![
            "GCMGCMG",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "GCM[+16]GCM[+16]G",
        ];
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_per_mod_limit_zero() {
        use ModificationSpecificity::*;
        // Limit of 0 means this mod is entirely suppressed
        let variable_mods = [(Residue(b'M'), 16.0f32, Some(0))];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, None);

        assert_eq!(peptides, vec!["GCMGCMG"]);
    }

    #[test]
    fn test_mixed_limited_and_unlimited() {
        use ModificationSpecificity::*;
        // M oxidation limited to 1; C carbamidomethylation unlimited
        // GCMGCMG has 2 M and 2 C
        let variable_mods = [
            (Residue(b'M'), 16.0f32, Some(1)),
            (Residue(b'C'), 57.0f32, None),
        ];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, None);

        // Should include all combos with ≤1 oxidized M,
        // but never both M residues oxidized simultaneously
        for p in &peptides {
            let oxid_count = p.matches("[+16]").count();
            assert!(oxid_count <= 1, "too many oxidations in: {}", p);
        }
        // Both C residues carbamidomethylated simultaneously should be present
        assert!(
            peptides.contains(&"GC[+57]MGC[+57]MG".to_string()),
            "expected double-C mod"
        );
        // Double oxidation should be absent
        assert!(
            !peptides.contains(&"GCM[+16]GCM[+16]G".to_string()),
            "double oxidation should be suppressed"
        );
    }

    #[test]
    fn test_limits_are_per_mod_not_per_residue() {
        use ModificationSpecificity::*;
        // Both modifications target M, but only oxidation is limited to one.
        let variable_mods = [
            (Residue(b'M'), 16.0f32, Some(1)),
            (Residue(b'M'), 32.0f32, None),
        ];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, None);

        assert!(!peptides.contains(&"GCM[+16]GCM[+16]G".to_string()));
        assert!(peptides.contains(&"GCM[+32]GCM[+32]G".to_string()));
    }

    #[test]
    fn test_limits_support_more_than_64_mod_entries() {
        use ModificationSpecificity::*;
        let mut variable_mods = (1..=65)
            .map(|mass| (Residue(b'M'), mass as f32, None))
            .collect::<Vec<_>>();
        variable_mods[64].2 = Some(0);
        let peptide = Peptide::try_from(Digest {
            sequence: "GMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 1, None);

        assert_eq!(peptides.len(), 65); // unmodified + 64 allowed entries
        assert!(!peptides.iter().any(|peptide| peptide.contains("[+65]")));
    }

    #[test]
    fn test_max_combinations_only_unmodified() {
        use ModificationSpecificity::*;
        // cap of 1 means only the unmodified peptide is returned
        let variable_mods = [(Residue(b'M'), 16.0f32, None)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, Some(1));

        assert_eq!(peptides, vec!["GCMGCMG"]);
    }

    #[test]
    fn test_max_combinations_prefers_fewer_ptms() {
        use ModificationSpecificity::*;
        // GCMGCMG with oxidation (2 sites) — normally 3 variants (unmod + 2 single + 1 double)
        // cap at 3 means we get unmod + both singles but not the double
        let variable_mods = [(Residue(b'M'), 16.0f32, None)];
        let peptide = Peptide::try_from(Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();

        let peptides = apply_masses(&peptide, &variable_mods, &[], 2, Some(3));

        assert_eq!(peptides, vec!["GCMGCMG", "GCM[+16]GCMG", "GCMGCM[+16]G"]);
        // Double-mod must not appear — it would require cap > 3
        assert!(!peptides.contains(&"GCM[+16]GCM[+16]G".to_string()));
    }

    #[test]
    fn named_mods_display() {
        use ModificationSpecificity::*;
        let peptide = Peptide::try_from(Digest {
            sequence: "PEMAT".into(),
            ..Default::default()
        })
        .unwrap();

        let mut mods = Modifications::default();
        let ox = mods
            .insert(15.9949, Some(&"UNIMOD:35".parse().unwrap()))
            .unwrap();
        let acetyl = mods.insert(42.0106, None).unwrap();
        let peptides = peptide.apply(
            &[(Residue(b'M'), ox, None)],
            &[(PeptideN(None), acetyl)],
            &mods,
            1,
            None,
        );
        let sequences = peptides
            .iter()
            .map(|p| p.display(&mods).to_string())
            .collect::<Vec<_>>();
        assert_eq!(
            sequences,
            vec!["[+42.0106]-PEMAT", "[+42.0106]-PEM[UNIMOD:35]AT"]
        );

        // Names don't change masses
        assert!((peptides[1].monoisotopic - peptides[0].monoisotopic - 15.9949).abs() < 1e-3);
    }
}
