use crate::enzyme::{Enzyme, EnzymeParameters};
use crate::fasta::Fasta;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON};
use crate::peptide::Peptide;
use log::error;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::path::PathBuf;

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use
    pub missed_cleavages: Option<u8>,
    /// Minimum peptide length that will be fragmented
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    pub max_len: Option<usize>,
    pub cleave_at: Option<String>,
    pub restrict: Option<char>,
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            cleave_at: Some("KR".into()),
            restrict: Some('P'),
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        EnzymeParameters {
            missed_cleavages: en.missed_cleavages.unwrap_or(0),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enyzme: Enzyme::new(
                &en.cleave_at.unwrap_or_else(|| "KR".into()),
                en.restrict.or(Some('P')),
            ),
        }
    }
}

#[derive(Deserialize, Default)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// This parameter allows tuning of the internal search structure
    pub bucket_size: Option<usize>,

    pub enzyme: Option<EnzymeBuilder>,
    /// Minimum fragment m/z that will be stored in the database
    pub fragment_min_mz: Option<f32>,
    /// Maximum fragment m/z that will be stored in the database
    pub fragment_max_mz: Option<f32>,
    /// Minimum peptide monoisotopic mass that will be fragmented
    pub peptide_min_mass: Option<f32>,
    /// Maximum peptide monoisotopic mass that will be fragmented
    pub peptide_max_mass: Option<f32>,
    /// Minimum ion index to be generated: 1 will remove b1/y1 ions
    /// 2 will remove b1/b2/y1/y2 ions, etc
    pub min_ion_index: Option<usize>,
    /// Static modifications to add to matching amino acids
    pub static_mods: Option<HashMap<char, f32>>,
    /// Variable modifications to add to matching amino acids
    pub variable_mods: Option<HashMap<char, f32>>,
    /// Limit number of variable modifications on a peptide
    pub max_variable_mods: Option<usize>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<PathBuf>,
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
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            static_mods: Self::validate_mods(self.static_mods),
            variable_mods: Self::validate_mods(self.variable_mods),
            max_variable_mods: self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2),
            generate_decoys: self.generate_decoys.unwrap_or(true),
            fasta: self.fasta.expect("A fasta file must be provided!"),
        }
    }

    pub fn update_fasta(&mut self, fasta: PathBuf) {
        self.fasta = Some(fasta)
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    bucket_size: usize,
    pub enzyme: EnzymeBuilder,
    pub fragment_min_mz: f32,
    pub fragment_max_mz: f32,
    peptide_min_mass: f32,
    peptide_max_mass: f32,
    min_ion_index: usize,
    static_mods: HashMap<char, f32>,
    variable_mods: HashMap<char, f32>,
    max_variable_mods: usize,
    decoy_tag: String,
    generate_decoys: bool,
    pub fasta: PathBuf,
}

impl Parameters {
    fn digest(
        &self,
        fasta: &Fasta,
        enzyme: &EnzymeParameters,
    ) -> (Vec<Peptide>, HashMap<String, Vec<String>>) {
        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        let digests = fasta.digest(enzyme);

        let mods = self
            .variable_mods
            .iter()
            .map(|(a, b)| (*a, *b))
            .collect::<Vec<_>>();

        // From our set of unique peptide sequence, apply any modifications
        // and convert to [`TargetDecoy`] enum
        let mut target_decoys = digests
            .par_iter()
            .filter_map(|(digest, _)| Peptide::try_from(digest).ok())
            .flat_map(|peptide| peptide.apply(&mods, &self.static_mods, self.max_variable_mods))
            .filter(|peptide| {
                peptide.monoisotopic >= self.peptide_min_mass
                    && peptide.monoisotopic <= self.peptide_max_mass
            })
            .collect::<Vec<Peptide>>();

        // NB: Stable sorting here (and only here?) is critical to achieving determinism...
        // not totally sure why... Probably has to do with using PeptideIxs
        // as keys for scoring vectors
        (&mut target_decoys).par_sort_by(|a, b| a.monoisotopic.total_cmp(&b.monoisotopic));

        let peptide_graph = digests
            .into_par_iter()
            .map(|(digest, proteins)| (digest.sequence, proteins))
            .collect::<HashMap<String, Vec<String>>>();

        (target_decoys, peptide_graph)
    }

    // pub fn build(self) -> Result<IndexedDatabase, Box<dyn std::error::Error + Send + Sync + 'static>> {
    pub fn build(self) -> std::io::Result<IndexedDatabase> {
        let enzyme = self.enzyme.clone().into();
        let fasta = Fasta::open(&self.fasta, self.decoy_tag.clone(), self.generate_decoys)?;

        let (target_decoys, peptide_graph) = self.digest(&fasta, &enzyme);

        // Finally, perform in silico digest for our target sequences
        // Note that multiple charge states are actually handled by
        // [`SpectrumProcessor`] or during scoring - all theoretical
        // fragments are monoisotopic/uncharged
        let mut fragments = target_decoys
            .par_iter()
            .enumerate()
            .flat_map_iter(|(idx, peptide)| {
                // Generate both B and Y ions, then filter down to make sure that
                // theoretical fragments are within the search space
                IonSeries::new(peptide, Kind::B)
                    .enumerate()
                    .chain(IonSeries::new(peptide, Kind::Y).enumerate())
                    .filter(|(ion_idx, ion)| {
                        // Don't store b1, b2, y1, y2 ions for preliminary scoring
                        let ion_idx_filter = match ion.kind {
                            Kind::B => (ion_idx + 1) > self.min_ion_index,
                            Kind::Y => {
                                peptide.sequence.len().saturating_sub(1) - ion_idx
                                    > self.min_ion_index
                            }
                        };
                        ion_idx_filter
                            && ion.monoisotopic_mass >= self.fragment_min_mz
                            && ion.monoisotopic_mass <= self.fragment_max_mz
                    })
                    .map(move |(_, ion)| Theoretical {
                        peptide_index: PeptideIx(idx as u32),
                        fragment_mz: ion.monoisotopic_mass,
                    })
            })
            .collect::<Vec<_>>();

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
                chunk.par_sort_unstable_by(|a, b| a.peptide_index.cmp(&b.peptide_index));
                min
            })
            .collect::<Vec<_>>();

        let mut potential_mods = self
            .static_mods
            .iter()
            .map(|(&x, &y)| (x, y))
            .collect::<Vec<(char, f32)>>();
        for (resi, mass) in &self.variable_mods {
            match self.static_mods.get(resi) {
                Some(mass_) if mass_ == mass => {}
                _ => {
                    potential_mods.push((*resi, *mass));
                }
            }
        }

        Ok(IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
            peptide_graph,
            potential_mods,
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
}

pub struct IndexedDatabase {
    pub peptides: Vec<Peptide>,
    pub fragments: Vec<Theoretical>,
    pub(crate) min_value: Vec<f32>,
    /// Keep a list of potential (AA, mass) modifications for RT prediction
    pub potential_mods: Vec<(char, f32)>,
    pub bucket_size: usize,
    peptide_graph: HashMap<String, Vec<String>>,
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
            |p, bounds| p.monoisotopic.total_cmp(bounds),
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

    pub fn assign_proteins(&self, peptide: &Peptide) -> (usize, String) {
        let sequence_without_mods = peptide
            .sequence
            .iter()
            .map(|resi| match resi {
                crate::mass::Residue::Just(r) => r,
                crate::mass::Residue::Mod(r, _) => r,
            })
            .collect::<String>();

        match self.peptide_graph.get(&sequence_without_mods) {
            Some(proteins) => {
                let n = proteins.len();
                let proteins = proteins.join(";");
                (n, proteins)
            }
            None => {
                panic!(
                    "BUG: peptide sequence {} doesn't appear in peptide graph!",
                    &sequence_without_mods
                );
            }
        }
    }
}

impl std::ops::Index<PeptideIx> for IndexedDatabase {
    type Output = Peptide;

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
                let neutral = self.db[frag.peptide_index].monoisotopic;
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
