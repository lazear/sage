use crate::fasta::{Fasta, Trypsin};
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::peptide::{Peptide, TargetDecoy};
use crate::spectrum::ProcessedSpectrum;
use log::error;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::path::PathBuf;

#[derive(Deserialize)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// This parameter allows tuning of the internal search structure
    bucket_size: Option<usize>,
    /// Minimum fragment m/z that will be stored in the database
    fragment_min_mz: Option<f32>,
    /// Maximum fragment m/z that will be stored in the database
    fragment_max_mz: Option<f32>,
    /// Minimum peptide length that will be fragmented
    peptide_min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    peptide_max_len: Option<usize>,
    /// Minimum peptide monoisotopic mass that will be fragmented
    peptide_min_mass: Option<f32>,
    /// Maximum peptide monoisotopic mass that will be fragmented
    peptide_max_mass: Option<f32>,
    /// Use target-decoy
    decoy: Option<bool>,
    /// How many missed cleavages to use
    missed_cleavages: Option<u8>,
    /// Static modification to add to the N-terminus of a peptide
    n_term_mod: Option<f32>,
    /// Static modifications to add to matching amino acids
    static_mods: Option<HashMap<char, f32>>,
    /// Path to fasta database
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
            fragment_max_mz: self.fragment_max_mz.unwrap_or(2500.0),
            peptide_min_len: self.peptide_min_len.unwrap_or(5),
            peptide_max_len: self.peptide_max_len.unwrap_or(65),
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
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
    peptide_min_mass: f32,
    peptide_max_mass: f32,
    decoy: bool,
    missed_cleavages: u8,
    n_term_mod: Option<f32>,
    static_mods: HashMap<char, f32>,
    pub fasta: PathBuf,
}

impl Parameters {
    fn digest(&self, fasta: &Fasta) -> Vec<TargetDecoy> {
        let trypsin = Trypsin::new(
            self.missed_cleavages,
            self.peptide_min_len,
            self.peptide_max_len,
        );

        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        //
        // Then, collect in a HashSet so that we only keep unique tryptic peptides
        let targets = fasta
            .proteins
            .par_iter()
            .flat_map(|(protein, sequence)| trypsin.digest(protein, sequence))
            .collect::<HashSet<_>>();

        // Handle reversed/decoy sequences separately
        // We want to make sure that there is no decoy peptide with the
        // exact sequence as a target peptide - so this leads to some
        // code duplication due to the dependency
        let reversed = fasta
            .proteins
            .iter()
            .map(|(acc, sequence)| {
                let sequence = sequence.chars().rev().collect::<String>();
                (format!("DECOY_{}", acc), sequence)
            })
            .collect::<Vec<_>>();

        let decoys = reversed
            .par_iter()
            .flat_map(|(protein, sequence)| trypsin.digest(protein, &sequence))
            .filter(|digest| !targets.contains(digest))
            .collect::<HashSet<_>>();

        log::trace!(
            "generated {} target and {} decoy peptides",
            targets.len(),
            decoys.len()
        );

        // From our set of unique peptide sequence, apply any modifications
        // and convert to [`TargetDecoy`] enum
        let mut target_decoys = targets
            .par_iter()
            .filter_map(|f| Peptide::try_from(f).ok())
            .filter(|p| {
                p.monoisotopic >= self.peptide_min_mass && p.monoisotopic <= self.peptide_max_mass
            })
            .map(|mut peptide| {
                // First modification we apply takes priority
                if let Some(m) = self.n_term_mod {
                    peptide.set_nterm_mod(m);
                }

                // Apply any relevant static modifications
                for (resi, mass) in &self.static_mods {
                    peptide.static_mod(*resi, *mass);
                }

                // If this is a reversed digest, annotate it as a Decoy
                TargetDecoy::Target(peptide)
            })
            .collect::<Vec<TargetDecoy>>();

        target_decoys.extend(
            decoys
                .par_iter()
                .filter_map(|f| Peptide::try_from(f).ok())
                .filter(|p| {
                    p.monoisotopic >= self.peptide_min_mass
                        && p.monoisotopic <= self.peptide_max_mass
                })
                .map(|mut peptide| {
                    // First modification we apply takes priority
                    if let Some(m) = self.n_term_mod {
                        peptide.set_nterm_mod(m);
                    }

                    // Apply any relevant static modifications
                    for (resi, mass) in &self.static_mods {
                        peptide.static_mod(*resi, *mass);
                    }

                    // If this is a reversed digest, annotate it as a Decoy
                    TargetDecoy::Decoy(peptide)
                })
                .collect::<Vec<TargetDecoy>>(),
        );

        (&mut target_decoys).par_sort_unstable_by(|a, b| a.neutral().total_cmp(&b.neutral()));
        target_decoys
    }

    pub fn build(self) -> Result<IndexedDatabase, Box<dyn std::error::Error>> {
        let fasta = Fasta::open(&self.fasta)?;
        let target_decoys = self.digest(&fasta);
        let mut fragments = Vec::new();

        // Finally, perform in silico digest for our target sequences
        // Note that multiple charge states are actually handled by the
        // [`SpectrumProcessor`], so we don't annotate charge states anywhere
        // else - we used to though, and could always revert to that - but it
        // saves a ton of memory to not generate 3-4x as many theoretical fragments
        for (idx, peptide) in target_decoys.iter().enumerate() {
            // Generate both B and Y ions, then filter down to make sure that
            // theoretical fragments are within the search space
            for kind in [Kind::B, Kind::Y] {
                fragments.extend(
                    IonSeries::new(peptide.peptide(), kind)
                        .map(|ion| Theoretical {
                            peptide_index: PeptideIx(idx as u32),
                            precursor_mz: peptide.neutral(),
                            fragment_mz: ion.monoisotopic_mass,
                            kind: ion.kind,
                        })
                        .filter(|frag| {
                            frag.fragment_mz >= self.fragment_min_mz
                                && frag.fragment_mz <= self.fragment_max_mz
                        }),
                );
            }
        }

        // Sort all of our theoretical fragments by m/z, from low to high
        (&mut fragments).par_sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

        // Now, we bucket all of our theoretical fragments, and within each bucket
        // sort by precursor m/z - and save the minimum *fragment* m/z in a separate
        // vector so that we can perform an efficient binary search to reduce
        // the number of in silico fragments we evaluate
        //
        // Imagine our theoretical fragments look like this
        //
        // Fragment        A      B       C       D       E       F       G       H
        // Fragment m/z [ 1.0    1.2     1.3     2.5     2.5     2.6     3.5     4.0 ]
        // Parent m/z   [ 500    439     291     800     142     515     517     232 ]
        //
        // If we apply a bucket size of 4 we will end up with the following:
        //
        // Fragment        C      B       A       D       E       H       F       G
        // Fragment m/z [ 1.3    1.2     1.0     2.5     2.5     4.0     3.5     2.6 ]
        // Parent m/z   [ 291    439     500     800     142     232     515     517 ]
        //              |___________________________|   |____________________________|
        //               Bucket 1: min m/z 1.0          Bucket 2: min m/z 2.5
        //
        // * Example query: Fragment m/z 1.3 - 1.9 & Precursor m/z: 450 - 900
        // 1) Perform a binary search to narrow down our window to Bucket 1 only
        //      * Bucket 2 has a min m/z outside of our query range - nothing here can match
        //
        // Fragment        C      B       A       D
        // Fragment m/z [ 1.3    1.2     1.0     2.5
        // Parent m/z   [ 291    439     500     800
        //                            |_____________|
        //                                    ^
        //                                    |
        // Window with matching precursors ___|

        // and within Bucket 1, we can perform another binary search to find fragments
        // matching our desired precursor m/z tolerance

        let min_value = fragments
            .par_chunks_mut(self.bucket_size)
            .map(|chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].fragment_mz;
                chunk.sort_unstable_by(|a, b| a.precursor_mz.total_cmp(&b.precursor_mz));
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

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
#[repr(transparent)]
pub struct PeptideIx(pub u32);

#[derive(Copy, Clone, Debug)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    // technically not m/z currently - precursor mono. mass
    pub precursor_mz: f32,
    pub fragment_mz: f32,
    pub kind: Kind,
}

pub struct IndexedDatabase {
    pub peptides: Vec<TargetDecoy>,
    pub fragments: Vec<Theoretical>,
    pub(crate) min_value: Vec<f32>,
    bucket_size: usize,
}

impl IndexedDatabase {
    /// Create a new [`IndexedQuery`] for a specific [`ProcessedSpectrum`]
    ///
    /// All matches returned by the query will be within the specified tolerance
    /// parameters
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

    pub fn buckets(&self) -> &[f32] {
        &self.min_value
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
    /// Search for a specified `fragment_mz` within the database
    pub fn page_search(&self, fragment_mz: f32) -> impl Iterator<Item = &Theoretical> {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(fragment_mz);
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.query.monoisotopic_mass);

        // Locate the left and right page indices that contain matching fragments
        // Note that we need to multiply by `bucket_size` to transform these into
        // indices that can be used with `self.db.fragments`
        let (left_idx, right_idx) =
            binary_search_slice(&self.db.min_value, |m| *m, fragment_lo, fragment_hi);

        // It is absolutely critical that we do not cross page boundaries!
        // If we do, we can no longer rely on total ordering of precursor m/z
        (left_idx..right_idx).flat_map(move |page| {
            let left_idx = page * self.db.bucket_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = ((page + 1) * self.db.bucket_size).min(self.db.fragments.len());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &&self.db.fragments[left_idx..right_idx];
            let (inner_left, inner_right) =
                binary_search_slice(slice, |frag| frag.precursor_mz, precursor_lo, precursor_hi);

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
                frag.precursor_mz >= precursor_lo
                    && frag.precursor_mz <= precursor_hi
                    && frag.fragment_mz >= fragment_lo
                    && frag.fragment_mz <= fragment_hi
            })
        })
    }
}

/// Return the widest `left` and `right` indices into a `slice` (sorted by the
/// function `key`) such that all values between `low` and `high` are
/// contained in `slice[left..right]`
///
/// # Invariants
///
/// * `slice[left] <= low || left == 0`
/// * `slice[right] <= high && (slice[right+1] > high || right == slice.len())`
/// * `0 <= left <= right <= slice.len()`
#[inline]
pub fn binary_search_slice<T, F>(slice: &[T], key: F, low: f32, high: f32) -> (usize, usize)
where
    F: Fn(&T) -> f32,
{
    let left_idx = match slice.binary_search_by(|a| key(a).total_cmp(&low)) {
        Ok(idx) | Err(idx) => {
            let mut idx = idx.saturating_sub(1);
            while idx > 0 && key(&slice[idx]) >= low {
                idx -= 1;
            }
            idx
        }
    };

    let right_idx = match slice[left_idx..].binary_search_by(|a| key(a).total_cmp(&high)) {
        Ok(idx) | Err(idx) => {
            let mut idx = idx + left_idx;
            while idx < slice.len() && key(&slice[idx]) <= high {
                idx = idx.saturating_add(1);
            }
            idx.min(slice.len())
        }
    };
    (left_idx, right_idx)
}

impl PeptideIx {
    pub fn for_testing_only_seriously_though(idx: usize) -> Self {
        Self(idx as u32)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn binary_search_slice_smoke() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
        let bounds = binary_search_slice(&data, |x| *x, 1.75, 3.5);
        assert_eq!(bounds, (1, 6));
        assert!(data[bounds.0] <= 1.75);
        assert_eq!(&data[bounds.0..bounds.1], &[1.5, 2.0, 2.5, 3.0, 3.5]);

        let bounds = binary_search_slice(&data, |x| *x, 0.0, 5.0);
        assert_eq!(bounds, (0, data.len()));
    }

    #[test]
    fn binary_search_slice_run() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0];
        let (left, right) = binary_search_slice(&data, |x| *x, 1.5, 3.25);
        assert!(data[left] <= 1.5);
        assert!(data[right] > 3.25);
        assert_eq!(
            &data[left..right],
            &[1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0]
        );
    }
}
