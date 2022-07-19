use crate::fasta::Trypsin;
use crate::peptide::{Peptide, TargetDecoy};
use crate::spectrum::ProcessedSpectrum;
use crate::{
    fasta::Fasta,
    ion_series::{IonSeries, Kind},
    mass::{Tolerance, PROTON},
};

use log::{error, info, warn};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::hash::Hash;

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

#[derive(Deserialize)]
pub struct Builder {
    // This parameter allows tuning of the internal search structure
    bucket_size: Option<usize>,
    // Minimum fragment m/z that will be stored in the database
    fragment_min_mz: Option<f32>,
    // Maximum fragment m/z that will be stored in the database
    fragment_max_mz: Option<f32>,
    // Minimum peptide length that will be fragmented
    peptide_min_len: Option<usize>,
    // Maximum peptide length that will be fragmented
    peptide_max_len: Option<usize>,
    // Use target-decoy
    decoy: Option<bool>,
    // How many missed cleavages to use
    missed_cleavages: Option<u8>,
    // Static modification to add to the N-terminus of a peptide
    n_term_mod: Option<f32>,
    // Static modifications to add to matching amino acids
    static_mods: Option<HashMap<char, f32>>,
    // Fasta database
    fasta: PathBuf,
}

impl Builder {
    pub fn make_parameters(self) -> Parameters {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        let mut static_mods = HashMap::new();
        if let Some(map) = self.static_mods {
            for (ch, mass) in map {
                if crate::mass::VALID_AA.contains(&ch) {
                    static_mods.insert(ch, mass);
                } else {
                    error!(
                        "invalid residue: {}, proceeding without this static mod",
                        ch
                    );
                }
            }
        }
        Parameters {
            bucket_size,
            fragment_min_mz: self.fragment_min_mz.unwrap_or(75.0),
            fragment_max_mz: self.fragment_max_mz.unwrap_or(4000.0),
            peptide_min_len: self.peptide_min_len.unwrap_or(5),
            peptide_max_len: self.peptide_max_len.unwrap_or(65),
            decoy: self.decoy.unwrap_or(true),
            missed_cleavages: self.missed_cleavages.unwrap_or(0),
            n_term_mod: self.n_term_mod,
            static_mods,
            fasta: self.fasta,
        }
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    bucket_size: usize,
    pub fragment_min_mz: f32,
    pub fragment_max_mz: f32,
    peptide_min_len: usize,
    peptide_max_len: usize,
    decoy: bool,
    missed_cleavages: u8,
    n_term_mod: Option<f32>,
    static_mods: HashMap<char, f32>,
    pub fasta: PathBuf,
}

impl Parameters {
    pub fn build(self) -> Result<IndexedDatabase, Box<dyn std::error::Error>> {
        let fasta = Fasta::open(self.fasta)?;

        let trypsin = Trypsin::new(
            self.decoy,
            self.missed_cleavages,
            self.peptide_min_len,
            self.peptide_max_len,
        );
        let peptides = fasta
            .proteins
            .par_iter()
            .flat_map(|(protein, sequence)| trypsin.digest(protein, sequence))
            .collect::<HashSet<_>>();

        let mut target_decoys = peptides
            .par_iter()
            .filter_map(|f| Peptide::try_from(f).ok().map(|pep| (f, pep)))
            .map(|(digest, mut peptide)| {
                // First modification we apply takes priority
                if let Some(m) = self.n_term_mod {
                    peptide.set_nterm_mod(m);
                }

                for (resi, mass) in &self.static_mods {
                    peptide.static_mod(*resi, *mass);
                }

                match digest.reversed {
                    true => TargetDecoy::Decoy(peptide),
                    false => TargetDecoy::Target(peptide),
                }
            })
            .collect::<Vec<TargetDecoy>>();

        // Sort by precursor neutral mass
        target_decoys.sort_by(|a, b| a.neutral().total_cmp(&b.neutral()));

        let mut fragments = Vec::new();

        for (idx, peptide) in target_decoys.iter().enumerate() {
            // for charge in 1..4 {
            for kind in [Kind::B, Kind::Y] {
                fragments.extend(
                    IonSeries::new(peptide.peptide(), kind)
                        .map(|ion| Theoretical {
                            peptide_index: PeptideIx(idx as u32),
                            precursor_mz: peptide.neutral(),
                            fragment_mz: ion.mz,
                            kind: ion.kind,
                        })
                        .filter(|frag| {
                            frag.fragment_mz >= self.fragment_min_mz
                                && frag.fragment_mz <= self.fragment_max_mz
                        }),
                );
            }
            // }
        }

        fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

        let min_value = fragments
            .par_chunks_mut(self.bucket_size)
            .map(|chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].fragment_mz;
                chunk.sort_by(|a, b| a.precursor_mz.total_cmp(&b.precursor_mz));
                min
            })
            .collect::<Vec<_>>();

        Ok(IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
        })
    }
}

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[repr(transparent)]
pub struct PeptideIx(u32);

#[derive(Copy, Clone)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub precursor_mz: f32,
    pub fragment_mz: f32,
    pub kind: Kind,
}

pub struct IndexedDatabase {
    pub(crate) peptides: Vec<TargetDecoy>,
    pub fragments: Vec<Theoretical>,
    pub(crate) min_value: Vec<f32>,
    bucket_size: usize,
}

impl IndexedDatabase {
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
        self.fragments.len()
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

        let (left_idx, right_idx) =
            binary_search_slice(&self.db.min_value, |m| *m, fragment_lo, fragment_hi);

        let left_idx = left_idx * self.db.bucket_size;
        // last chunk not guaranted to be modulo bucket size
        let right_idx = (right_idx * self.db.bucket_size).min(self.db.fragments.len());

        let (left, right) = self.precursor_tol.bounds(self.query.precursor_mz - PROTON);

        let slice = &&self.db.fragments[left_idx..right_idx];

        let (inner_left, inner_right) =
            binary_search_slice(slice, |frag| frag.precursor_mz, left, right);
        slice[inner_left..inner_right].iter().filter(move |frag| {
            frag.precursor_mz >= left
                && frag.precursor_mz <= right
                && frag.fragment_mz >= fragment_lo
                && frag.fragment_mz <= fragment_hi
        })
    }
}

#[inline]
pub fn binary_search_slice<T, F>(slice: &[T], key: F, left: f32, right: f32) -> (usize, usize)
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

    let right_idx = match slice[left_idx..].binary_search_by(|a| key(a).total_cmp(&right)) {
        Ok(idx) | Err(idx) => (left_idx + idx).min(slice.len().saturating_sub(1)),
    };
    (left_idx, right_idx)
}
