use crate::fasta::{Fasta, Trypsin};
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON};
use crate::peptide::{Peptide, TargetDecoy};
use log::error;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
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
    /// How many missed cleavages to use
    missed_cleavages: Option<u8>,
    /// Static modifications to add to matching amino acids
    static_mods: Option<HashMap<char, f32>>,
    /// Variable modifications to add to matching amino acids
    variable_mods: Option<HashMap<char, f32>>,
    /// Use this prefix for decoy proteins
    decoy_prefix: Option<String>,
    /// Path to fasta database
    fasta: PathBuf,
}

impl Builder {
    fn validate_mods(input: Option<HashMap<char, f32>>) -> HashMap<char, f32> {
        let mut output = HashMap::new();
        if let Some(input) = input {
            for (ch, mass) in input {
                if crate::mass::VALID_AA.contains(&ch) || ch == '^' {
                    output.insert(ch, mass);
                } else {
                    error!(
                        "invalid residue: {}, proceeding without this static mod",
                        ch
                    );
                }
            }
        }
        output
    }

    pub fn make_parameters(self) -> Parameters {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        Parameters {
            bucket_size,
            fragment_min_mz: self.fragment_min_mz.unwrap_or(150.0),
            fragment_max_mz: self.fragment_max_mz.unwrap_or(2000.0),
            peptide_min_len: self.peptide_min_len.unwrap_or(5),
            peptide_max_len: self.peptide_max_len.unwrap_or(50),
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            decoy_prefix: self.decoy_prefix.unwrap_or_else(|| "rev_".into()),
            missed_cleavages: self.missed_cleavages.unwrap_or(0),
            static_mods: Self::validate_mods(self.static_mods),
            variable_mods: Self::validate_mods(self.variable_mods),
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
    missed_cleavages: u8,
    static_mods: HashMap<char, f32>,
    variable_mods: HashMap<char, f32>,
    decoy_prefix: String,
    pub fasta: PathBuf,
}

impl Parameters {
    fn digest(&self, fasta: &Fasta, trypsin: &Trypsin) -> Vec<TargetDecoy> {
        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        //
        // Then, collect in a HashSet so that we only keep unique tryptic peptides
        let targets = fasta
            .targets
            .par_iter()
            .flat_map(|(protein, sequence)| trypsin.digest(protein, sequence, false))
            .chain(
                fasta
                    .decoys
                    .par_iter()
                    .flat_map(|(protein, sequence)| trypsin.digest(protein, sequence, true)),
            )
            .collect::<HashSet<_>>();

        // From our set of unique peptide sequence, apply any modifications
        // and convert to [`TargetDecoy`] enum
        let mut target_decoys = targets
            .par_iter()
            .filter_map(|f| Peptide::try_from(f).ok())
            .flat_map(|peptide| peptide.apply(&self.variable_mods, &self.static_mods))
            .filter(|peptide| {
                peptide.monoisotopic >= self.peptide_min_mass
                    && peptide.monoisotopic <= self.peptide_max_mass
            })
            .map(
                |peptide| match peptide.protein.starts_with(&self.decoy_prefix) {
                    true => TargetDecoy::Decoy(peptide),
                    false => TargetDecoy::Target(peptide),
                },
            )
            .collect::<Vec<TargetDecoy>>();
        (&mut target_decoys).par_sort_unstable_by(|a, b| a.neutral().total_cmp(&b.neutral()));
        target_decoys
    }

    pub fn build(self) -> Result<IndexedDatabase, Box<dyn std::error::Error>> {
        let trypsin = Trypsin::new(
            self.missed_cleavages,
            self.peptide_min_len,
            self.peptide_max_len,
        );
        let mut fasta = Fasta::open(&self.fasta, &self.decoy_prefix)?;
        fasta.make_decoys(&self.decoy_prefix);

        let target_decoys = self.digest(&fasta, &trypsin);
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
                            // precursor_mz: peptide.neutral(),
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
        (&mut fragments).par_sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

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
                chunk.par_sort_unstable_by(|a, b| a.peptide_index.cmp(&b.peptide_index));
                min
            })
            .collect::<Vec<_>>();

        Ok(IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
            fasta,
            trypsin,
        })
    }
}

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize)]
#[repr(transparent)]
pub struct PeptideIx(pub u32);

// This is unsafe for use outside of this crate
impl Default for PeptideIx {
    fn default() -> Self {
        Self(u32::MAX)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub fragment_mz: f32,
    pub kind: Kind,
}

pub struct IndexedDatabase {
    pub peptides: Vec<TargetDecoy>,
    pub fragments: Vec<Theoretical>,
    pub(crate) min_value: Vec<f32>,
    bucket_size: usize,
    fasta: Fasta,
    trypsin: Trypsin,
}

impl IndexedDatabase {
    /// Create a new [`IndexedQuery`] for a specific [`ProcessedSpectrum`]
    ///
    /// All matches returned by the query will be within the specified tolerance
    /// parameters
    pub fn query(
        &self,
        precursor_mass: f32,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
        min_isotope_err: i8,
        max_isotope_err: i8,
    ) -> IndexedQuery<'_> {
        let (precursor_lo, precursor_hi) = precursor_tol.bounds(precursor_mass);

        let (pre_idx_lo, pre_idx_hi) = binary_search_slice(
            &self.peptides,
            |p, bounds| p.neutral().total_cmp(bounds),
            precursor_lo - (NEUTRON * max_isotope_err as f32).abs(),
            precursor_hi + (NEUTRON * min_isotope_err as f32).abs(),
        );

        IndexedQuery {
            db: self,
            precursor_mass,
            precursor_tol,
            fragment_tol,
            min_isotope_err,
            max_isotope_err,
            pre_idx_lo,
            pre_idx_hi,
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

pub struct IndexedQuery<'d> {
    db: &'d IndexedDatabase,
    precursor_mass: f32,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    min_isotope_err: i8,
    max_isotope_err: i8,
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

impl<'d> IndexedQuery<'d> {
    /// Search for a specified `fragment_mz` within the database
    pub fn page_search(&self, fragment_mz: f32) -> impl Iterator<Item = &Theoretical> {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(fragment_mz);
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.precursor_mass);

        // Locate the left and right page indices that contain matching fragments
        // Note that we need to multiply by `bucket_size` to transform these into
        // indices that can be used with `self.db.fragments`
        let (left_idx, right_idx) = binary_search_slice(
            &self.db.min_value,
            |min, bounds| min.total_cmp(bounds),
            fragment_lo,
            fragment_hi,
        );

        // It is absolutely critical that we do not cross page boundaries!
        // If we do, we can no longer rely on total ordering of peptide_index (precursor m/z)
        (left_idx..right_idx).flat_map(move |page| {
            let left_idx = page * self.db.bucket_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = ((page + 1) * self.db.bucket_size).min(self.db.fragments.len());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &&self.db.fragments[left_idx..right_idx];

            let (inner_left, inner_right) = binary_search_slice(
                slice,
                |frag, bounds| (frag.peptide_index.0 as usize).cmp(bounds),
                self.pre_idx_lo,
                self.pre_idx_hi,
            );

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
                let neutral = self.db[frag.peptide_index].neutral();
                (self.min_isotope_err..=self.max_isotope_err).any(|isotope_err| {
                    let delta = isotope_err as f32 * NEUTRON;
                    (neutral >= precursor_lo - delta) && (neutral <= precursor_hi - delta)
                }) && frag.fragment_mz >= fragment_lo
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
pub fn binary_search_slice<T, F, S>(slice: &[T], key: F, low: S, high: S) -> (usize, usize)
where
    F: Fn(&T, &S) -> Ordering,
{
    let left_idx = match slice.binary_search_by(|a| key(a, &low)) {
        Ok(idx) | Err(idx) => {
            let mut idx = idx.saturating_sub(1);
            while idx > 0 && key(&slice[idx], &low) != Ordering::Less {
                idx -= 1;
            }
            idx
        }
    };

    let right_idx = match slice[left_idx..].binary_search_by(|a| key(a, &high)) {
        Ok(idx) | Err(idx) => {
            let mut idx = idx + left_idx;
            while idx < slice.len() && key(&slice[idx], &high) != Ordering::Greater {
                idx = idx.saturating_add(1);
            }
            idx.min(slice.len())
        }
    };
    (left_idx, right_idx)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn binary_search_slice_smoke() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.75, 3.5);
        assert_eq!(bounds, (1, 6));
        assert!(data[bounds.0] <= 1.75);
        assert_eq!(&data[bounds.0..bounds.1], &[1.5, 2.0, 2.5, 3.0, 3.5]);

        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 0.0, 5.0);
        assert_eq!(bounds, (0, data.len()));
    }

    #[test]
    fn binary_search_slice_run() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0];
        let (left, right) = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.5, 3.25);
        assert!(data[left] <= 1.5);
        assert!(data[right] > 3.25);
        assert_eq!(
            &data[left..right],
            &[1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0]
        );
    }
}
