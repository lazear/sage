use crate::fasta::Trypsin;
use crate::mass::{Modification, Residue};
use crate::peptide::{Peptide, TargetDecoy};
use crate::spectrum::ProcessedSpectrum;
use crate::{
    fasta::Fasta,
    ion_series::{IonSeries, Kind},
    mass::{Tolerance, PROTON},
};

use rayon::prelude::*;
use std::hash::Hash;

use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

pub const FRAGMENT_BUCKET_SIZE: usize = 256;

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[repr(transparent)]
pub struct PeptideIx(u32);

#[derive(Copy, Clone)]
#[repr(C)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub precursor_mz: f32,
    pub fragment_mz: f32,
    pub kind: Kind,
    pub charge: u8,
}

#[derive(Copy, Clone)]
#[repr(C)]
pub(crate) struct FragmentPage {
    pub(crate) fragments: [Theoretical; FRAGMENT_BUCKET_SIZE],
}

pub struct IndexedDatabase {
    pub(crate) peptides: Vec<TargetDecoy>,
    pub(crate) pages: Vec<FragmentPage>,
    pub(crate) min_value: Vec<f32>,
    pub fragment_min_mz: f32,
    pub fragment_max_mz: f32,
}

impl IndexedDatabase {
    pub fn new<P: AsRef<Path>>(
        p: P,
        static_mods: HashMap<Residue, Modification>,
        fragment_min_mz: f32,
        fragment_max_mz: f32,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let fasta = Fasta::open(p)?;

        let trypsin = Trypsin::new(true, true);
        let peptides = fasta
            .proteins
            .par_iter()
            .flat_map(|(protein, sequence)| trypsin.digest(protein, sequence))
            .filter(|dig| dig.sequence.len() >= 7 && dig.sequence.len() <= 50)
            .collect::<HashSet<_>>();

        let mut target_decoys = peptides
            .par_iter()
            .filter_map(|f| Peptide::try_from(f).ok().map(|pep| (f, pep)))
            .map(|(digest, mut peptide)| {
                for (resi, modi) in &static_mods {
                    peptide.static_mod(resi, *modi);
                }
                match digest.reversed {
                    true => TargetDecoy::Decoy(peptide),
                    false => TargetDecoy::Target(peptide),
                }
                // [TargetDecoy::Target(peptide), TargetDecoy::Decoy(reversed)]
            })
            .collect::<Vec<TargetDecoy>>();

        // Sort by precursor neutral mass
        target_decoys.sort_by(|a, b| a.neutral().total_cmp(&b.neutral()));

        let mut fragments = Vec::new();

        for (idx, peptide) in target_decoys.iter().enumerate() {
            for charge in 1..4 {
                for kind in [Kind::B, Kind::Y] {
                    fragments.extend(
                        IonSeries::new(peptide.peptide(), kind, charge)
                            .map(|ion| Theoretical {
                                peptide_index: PeptideIx(idx as u32),
                                precursor_mz: peptide.neutral(),
                                fragment_mz: ion.mz,
                                kind: ion.kind,
                                charge: ion.charge,
                            })
                            .filter(|frag| {
                                frag.fragment_mz >= fragment_min_mz
                                    && frag.fragment_mz <= fragment_max_mz
                            }),
                    );
                }
            }
        }

        fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

        let (pages, min_value): (Vec<FragmentPage>, Vec<f32>) = fragments
            .par_chunks(FRAGMENT_BUCKET_SIZE)
            .map(|chunk| {
                let default = Theoretical {
                    peptide_index: PeptideIx(0),
                    precursor_mz: f32::MAX,
                    fragment_mz: f32::MAX,
                    kind: Kind::B,
                    charge: u8::MAX,
                };

                let mut array = [default; FRAGMENT_BUCKET_SIZE];
                array[..chunk.len().min(FRAGMENT_BUCKET_SIZE)].copy_from_slice(chunk);
                array.sort_by(|a, b| a.precursor_mz.total_cmp(&b.precursor_mz));
                let m = array
                    .iter()
                    .fold(f32::MAX, |acc, frag| acc.min(frag.fragment_mz));
                (FragmentPage { fragments: array }, m)
            })
            .unzip();

        Ok(Self {
            peptides: target_decoys,
            pages,
            min_value,
            fragment_max_mz,
            fragment_min_mz,
        })
    }

    pub fn query<'d, 'q>(
        &'d self,
        query: &'q ProcessedSpectrum,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
    ) -> IndexedQuery<'d, 'q> {
        IndexedQuery {
            db: self,
            query,
            precursor_tol,
            fragment_tol,
        }
    }

    pub fn size(&self) -> usize {
        self.pages.len() * FRAGMENT_BUCKET_SIZE
    }
}

impl std::ops::Index<PeptideIx> for IndexedDatabase {
    type Output = TargetDecoy;

    fn index(&self, index: PeptideIx) -> &Self::Output {
        &self.peptides[index.0 as usize]
    }
}

pub struct IndexedQuery<'d, 'q> {
    db: &'d IndexedDatabase,
    query: &'q ProcessedSpectrum,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
}

impl<'d, 'q> IndexedQuery<'d, 'q> {
    pub fn page_search(&self, fragment_mz: f32) -> impl Iterator<Item = &Theoretical> {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(fragment_mz);

        let left_idx = match self
            .db
            .min_value
            .binary_search_by(|a| a.total_cmp(&fragment_lo))
        {
            Ok(idx) | Err(idx) => {
                let mut idx = idx.saturating_sub(1);
                while idx > 0 && self.db.min_value[idx] >= fragment_lo {
                    idx -= 1;
                }
                idx
            }
        };

        let right_idx = match self
            .db
            .min_value
            .binary_search_by(|a| a.total_cmp(&fragment_hi))
        {
            Ok(idx) | Err(idx) => idx.min(self.db.min_value.len().saturating_sub(1)),
        };

        debug_assert!(self.db.min_value[left_idx] <= fragment_lo || left_idx == 0);
        debug_assert!(
            self.db.min_value[right_idx] >= fragment_hi || right_idx == self.db.min_value.len() - 1
        );

        self.db.pages[left_idx..right_idx]
            .iter()
            .flat_map(move |page| {
                let (left, right) = self.precursor_tol.bounds(self.query.precursor_mz - PROTON);
                binary_search_slice(&page.fragments, |frag| frag.precursor_mz, left, right)
                    .iter()
                    .filter(move |frag| {
                        frag.precursor_mz >= left
                            && frag.precursor_mz <= right
                            && frag.fragment_mz >= fragment_lo
                            && frag.fragment_mz <= fragment_hi
                    })
            })
    }
}

#[inline]
pub fn binary_search_slice<T, F>(slice: &[T], key: F, left: f32, right: f32) -> &[T]
where
    F: Fn(&T) -> f32,
{
    let left_idx = match slice.binary_search_by(|a| key(a).total_cmp(&left)) {
        Ok(idx) | Err(idx) => {
            let mut idx = idx.saturating_sub(1);
            while idx > 0 && key(&slice[idx]) >= left {
                idx -= 1;
            }
            idx
        }
    };

    let right_idx = match slice.binary_search_by(|a| key(a).total_cmp(&right)) {
        Ok(idx) | Err(idx) => idx.min(slice.len().saturating_sub(1)),
    };

    &slice[left_idx..right_idx]
}
