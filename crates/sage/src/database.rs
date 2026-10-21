use crate::enzyme::{group_digests, Enzyme, EnzymeParameters};
use crate::fasta::Fasta;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{
    DetailedMod, ModEntry, ModId, ModificationError, ModificationSpecificity, Modifications,
};
use crate::peptide::Peptide;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::hash::Hash;

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use
    pub missed_cleavages: Option<u8>,
    /// Minimum peptide length that will be fragmented
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    pub max_len: Option<usize>,
    pub cleave_at: Option<String>,
    pub restrict: Option<String>,
    pub c_terminal: Option<bool>,
    pub semi_enzymatic: Option<bool>,
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            semi_enzymatic: Some(false),
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        EnzymeParameters {
            missed_cleavages: en.missed_cleavages.unwrap_or(1),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enzyme: Enzyme::new(
                &en.cleave_at.unwrap_or_else(|| "KR".into()),
                &en.restrict.unwrap_or_else(|| "".into()),
                en.c_terminal.unwrap_or(true),
                en.semi_enzymatic.unwrap_or(false),
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
    /// Minimum peptide monoisotopic mass that will be fragmented
    pub peptide_min_mass: Option<f32>,
    /// Maximum peptide monoisotopic mass that will be fragmented
    pub peptide_max_mass: Option<f32>,
    /// Which kind of fragment ions to generate (a, b, c, x, y, z, z_radical)
    pub ion_kinds: Option<Vec<Kind>>,
    /// Minimum ion index to be generated: 1 will remove b1/y1 ions
    /// 2 will remove b1/b2/y1/y2 ions, etc
    pub min_ion_index: Option<usize>,
    /// Static modifications to add to matching amino acids.
    /// Each entry is either a bare mass (`57.0215`) or an object with `mass` and
    /// an optional ProForma `name` (`{"mass": 57.0215, "name": "UNIMOD:4"}`).
    pub static_mods: Option<HashMap<ModificationSpecificity, ModEntry>>,
    /// Variable modifications to add to matching amino acids.
    /// Each entry is either a bare mass (`15.9949`) or an object with `mass` and
    /// optional `max_count` and `name` fields
    /// (`{"mass": 15.9949, "max_count": 1, "name": "UNIMOD:35"}`).
    pub variable_mods: Option<HashMap<ModificationSpecificity, Vec<ModEntry>>>,
    /// Limit number of variable modifications on a peptide
    pub max_variable_mods: Option<usize>,
    /// Hard cap on the total peptide variants generated per input peptide,
    /// including its unmodified form. Values below 1 are normalized to 1.
    /// Variants with fewer PTMs are preferred (generated first).
    pub max_combinations: Option<usize>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
    /// Number of sequences to handle simultaneously when pre-filtering the db
    pub prefilter_chunk_size: Option<usize>,
    /// Pre-filter the database to minimize memory usage
    pub prefilter: Option<bool>,
    /// Pre-filter the database with a minimal amount of memory at the cost of speed
    pub prefilter_low_memory: Option<bool>,
}

impl Builder {
    pub fn make_parameters(self) -> Result<Parameters, ModificationError> {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        let modifications = ResolvedModifications::new(
            &self.static_mods.unwrap_or_default(),
            &self.variable_mods.unwrap_or_default(),
        )?;
        Ok(Parameters {
            bucket_size,
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            modifications,
            max_variable_mods: self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2),
            max_combinations: self.max_combinations.map(|x| x.max(1)),
            generate_decoys: self.generate_decoys.unwrap_or(true),
            fasta: self.fasta.expect("A fasta file must be provided!"),
            prefilter_chunk_size: self.prefilter_chunk_size.unwrap_or(0),
            prefilter: self.prefilter.unwrap_or(false),
            prefilter_low_memory: self.prefilter_low_memory.unwrap_or(true),
        })
    }

    pub fn update_fasta(&mut self, fasta: String) {
        self.fasta = Some(fasta)
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    pub bucket_size: usize,
    pub enzyme: EnzymeBuilder,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub ion_kinds: Vec<Kind>,
    pub min_ion_index: usize,
    /// Serialized as the `static_mods` and `variable_mods` it was built from
    #[serde(flatten)]
    pub modifications: ResolvedModifications,
    pub max_variable_mods: usize,
    pub max_combinations: Option<usize>,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub fasta: String,
    pub prefilter_chunk_size: usize,
    pub prefilter: bool,
    pub prefilter_low_memory: bool,
}

impl Parameters {
    pub fn auto_calculate_prefilter_chunk_size(&mut self, fasta: &Fasta) {
        const MAX_PEPS_PER_CHUNK: usize = 2usize.pow(23);
        self.prefilter_chunk_size = match self.prefilter_chunk_size {
            0 => {
                let enzyme = self.enzyme.clone().into();
                let total_unmodified_pep_count: usize = fasta.digest(&enzyme).len();
                let combination_factor = if self.max_variable_mods >= usize::BITS as usize {
                    usize::MAX
                } else {
                    1usize << self.max_variable_mods
                };
                let mut mod_count_estimate = (self.modifications.variable_specificities() + 1)
                    .saturating_mul(combination_factor);
                if let Some(max_combinations) = self.max_combinations {
                    mod_count_estimate = mod_count_estimate.min(max_combinations);
                }
                let chunk_count = mod_count_estimate.saturating_mul(total_unmodified_pep_count)
                    / MAX_PEPS_PER_CHUNK;
                fasta
                    .targets
                    .len()
                    .checked_div(chunk_count)
                    .unwrap_or(fasta.targets.len())
            }
            x => x,
        };
    }

    pub fn digest(&self, fasta: &Fasta) -> Vec<Peptide> {
        log::trace!("digesting fasta");
        let enzyme = self.enzyme.clone().into();
        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        let digests = fasta.digest(&enzyme);

        log::trace!("grouping digests");
        let start_num = digests.len();
        let digests = group_digests(digests);
        log::trace!(
            "grouped {} digests into {} groups",
            start_num,
            digests.len()
        );

        let mods = &self.modifications;

        let targets: DashSet<_, FnvBuildHasher> = DashSet::default();
        digests
            .par_iter()
            .filter(|digest| !digest.reference.decoy)
            .for_each(|digest| {
                targets.insert(digest.reference.sequence.clone().into_bytes());
            });

        log::trace!("modifying peptides");
        let mut target_decoys = digests
            .into_par_iter()
            .map(Peptide::try_from)
            .filter_map(Result::ok)
            .flat_map_iter(|peptide| {
                peptide
                    .apply(
                        &mods.variable_mods,
                        &mods.static_mods,
                        &mods.registry,
                        self.max_variable_mods,
                        self.max_combinations,
                    )
                    .into_iter()
                    .filter(|peptide| {
                        peptide.monoisotopic >= self.peptide_min_mass
                            && peptide.monoisotopic <= self.peptide_max_mass
                    })
                    .flat_map(|peptide| {
                        if self.generate_decoys {
                            vec![peptide.reverse(), peptide].into_iter()
                        } else {
                            vec![peptide].into_iter()
                        }
                    })
                    .filter(|peptide| !peptide.decoy || !targets.contains(&(peptide.sequence[..])))
            })
            .collect::<Vec<_>>();

        Self::reorder_peptides(&mut target_decoys);

        target_decoys
    }

    pub fn reorder_peptides(target_decoys: &mut Vec<Peptide>) {
        log::trace!("sorting and deduplicating peptides");

        let init_size = target_decoys.len();
        // This is equivalent to a stable sort
        target_decoys.par_sort_unstable_by(|a, b| {
            a.monoisotopic
                .total_cmp(&b.monoisotopic)
                .then_with(|| a.initial_sort(b))
        });
        target_decoys.dedup_by(|remove, keep| {
            if remove.monoisotopic == keep.monoisotopic
                && remove.sequence == keep.sequence
                && remove.modifications == keep.modifications
                && remove.nterm == keep.nterm
                && remove.cterm == keep.cterm
            {
                keep.proteins.extend(remove.proteins.iter().cloned());
                // When merging peptides from different Fastas,
                // decoys in one fasta might be targets in another
                keep.decoy &= remove.decoy;
                true
            } else {
                false
            }
        });

        target_decoys
            .par_iter_mut()
            .for_each(|peptide| peptide.proteins.sort_unstable());

        let num_dropped = init_size - target_decoys.len();
        log::trace!(
            "dropped {} t/d pairs, remaining {}",
            num_dropped,
            target_decoys.len(),
        );
    }

    pub fn build(self, fasta: Fasta) -> IndexedDatabase {
        let target_decoys = self.digest(&fasta);
        self.build_from_peptides(target_decoys)
    }

    pub fn build_from_peptides(self, target_decoys: Vec<Peptide>) -> IndexedDatabase {
        log::trace!("generating fragments");

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
                self.ion_kinds
                    .iter()
                    .flat_map(|kind| {
                        IonSeries::new(peptide, *kind, &self.modifications.registry).enumerate()
                    })
                    .filter(|(ion_idx, ion)| {
                        // Don't store b1, b2, y1, y2 ions for preliminary scoring

                        match ion.kind {
                            Kind::A | Kind::B | Kind::C => (ion_idx + 1) > self.min_ion_index,
                            Kind::X | Kind::Y | Kind::Z | Kind::ZRadical => {
                                peptide.sequence.len().saturating_sub(1) - ion_idx
                                    > self.min_ion_index
                            }
                        }
                    })
                    .map(move |(_, ion)| Theoretical {
                        peptide_index: PeptideIx(idx as u32),
                        fragment_mz: ion.monoisotopic_mass,
                    })
            })
            .collect::<Vec<_>>();
        log::trace!("finalizing index");

        // Sort all of our theoretical fragments by m/z, from low to high
        fragments.par_sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

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

        IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            modifications: self.modifications.registry,
            decoy_tag: self.decoy_tag,
        }
    }
}

/// Configured modifications, as ids into `registry`
#[derive(Clone, Debug, Default)]
pub struct ResolvedModifications {
    pub registry: Modifications,
    pub static_mods: Vec<(ModificationSpecificity, ModId)>,
    /// `(specificity, modification, per_mod_limit)`, in a stable order
    pub variable_mods: Vec<(ModificationSpecificity, ModId, Option<usize>)>,
}

impl ResolvedModifications {
    /// Number of distinct sites (e.g. `M`, `^`) with variable modifications
    fn variable_specificities(&self) -> usize {
        // `variable_mods` is sorted by specificity
        self.variable_mods
            .iter()
            .map(|(specificity, _, _)| specificity)
            .dedup()
            .count()
    }

    /// Serializable config entry for a registered modification
    fn entry(&self, id: ModId, max_count: Option<usize>) -> ModEntry {
        let modification = self.registry.get(id);
        match (&modification.name, max_count) {
            (None, None) => ModEntry::Mass(modification.mass),
            (name, max_count) => ModEntry::Detailed(DetailedMod {
                mass: modification.mass,
                max_count,
                name: name.clone(),
            }),
        }
    }

    /// Register the configured modifications. Ids depend only on the
    /// configuration, so every database built from the same parameters (e.g.
    /// each prefilter chunk) agrees on them.
    pub fn new(
        static_mods: &HashMap<ModificationSpecificity, ModEntry>,
        variable_mods: &HashMap<ModificationSpecificity, Vec<ModEntry>>,
    ) -> Result<Self, ModificationError> {
        let mut registry = Modifications::default();

        if let Some((specificity, _)) = static_mods
            .iter()
            .find(|(_, entry)| entry.max_count().is_some())
        {
            return Err(ModificationError::StaticMaxCount(*specificity));
        }

        let mut static_entries = static_mods.iter().collect::<Vec<_>>();
        static_entries.sort_unstable_by_key(|(specificity, _)| **specificity);
        let static_mods = static_entries
            .into_iter()
            .map(|(specificity, entry)| {
                Ok((*specificity, registry.insert(entry.mass(), entry.name())?))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Flatten variable modifications into a stable order. This matters when
        // `max_combinations` truncates variants: equivalent configurations must
        // retain the same variants regardless of randomized `HashMap` iteration.
        let mut variable_entries = variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries
                    .iter()
                    .enumerate()
                    .map(|(entry_order, entry)| (*specificity, entry_order, entry))
            })
            .collect::<Vec<_>>();
        variable_entries.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        let variable_mods = variable_entries
            .into_iter()
            .map(|(specificity, _, entry)| {
                let id = registry.insert(entry.mass(), entry.name())?;
                Ok((specificity, id, entry.max_count()))
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(ResolvedModifications {
            registry,
            static_mods,
            variable_mods,
        })
    }
}

impl Serialize for ResolvedModifications {
    /// Same shape as the `static_mods` and `variable_mods` config fields
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;

        let static_mods = self
            .static_mods
            .iter()
            .map(|&(specificity, id)| (specificity, self.entry(id, None)))
            .collect::<BTreeMap<_, _>>();
        let mut variable_mods = BTreeMap::<_, Vec<_>>::new();
        for &(specificity, id, max_count) in &self.variable_mods {
            variable_mods
                .entry(specificity)
                .or_default()
                .push(self.entry(id, max_count));
        }

        let mut state = serializer.serialize_struct("ResolvedModifications", 2)?;
        state.serialize_field("static_mods", &static_mods)?;
        state.serialize_field("variable_mods", &variable_mods)?;
        state.end()
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

#[derive(Default)]
pub struct IndexedDatabase {
    pub peptides: Vec<Peptide>,
    pub fragments: Vec<Theoretical>,
    pub ion_kinds: Vec<Kind>,
    pub min_value: Vec<f32>,
    pub modifications: Modifications,
    pub bucket_size: usize,
    pub generate_decoys: bool,
    pub decoy_tag: String,
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
    ) -> IndexedQuery<'_> {
        let (precursor_lo, precursor_hi) = precursor_tol.bounds(precursor_mass);

        let (pre_idx_lo, pre_idx_hi) = binary_search_slice(
            &self.peptides,
            |p, bounds| p.monoisotopic.total_cmp(bounds),
            precursor_lo,
            precursor_hi,
        );

        IndexedQuery {
            db: self,
            precursor_mass,
            precursor_tol,
            fragment_tol,
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

    pub fn serialize(&self) {
        use std::io::Write;
        let mut wtr = std::io::BufWriter::new(std::fs::File::create("fragments.bin").unwrap());
        for fragment in &self.fragments {
            let _ = wtr.write(&fragment.fragment_mz.to_le_bytes()).unwrap();
            let _ = wtr.write(&fragment.peptide_index.0.to_le_bytes()).unwrap();
        }
        wtr.flush().unwrap();

        let mut wtr = std::io::BufWriter::new(std::fs::File::create("peptides.csv").unwrap());
        writeln!(wtr, "peptide,proteins,monoisotopic,decoy").unwrap();
        for fragment in &self.peptides {
            writeln!(
                wtr,
                "{},{},{},{}",
                fragment.display(&self.modifications),
                fragment.proteins(&self.decoy_tag, self.generate_decoys),
                fragment.monoisotopic,
                fragment.decoy
            )
            .unwrap();
        }
        wtr.flush().unwrap();
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
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

impl IndexedQuery<'_> {
    /// Search for a specified `fragment_mz` within the database
    pub fn page_search(&self, mass: f32) -> impl Iterator<Item = &Theoretical> {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(mass);
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
                // This looks somewhat complicated, but it's a consequence of
                // how the `binary_search_slice` function works - it will return
                // the set of indices that maximally cover the desired range - the exact
                // `left` and `right` indices may be valid, or just outside of the range.
                // Anything interior of `left` and `right` is guaranteed to be within the
                // precursor tolerance, so we just need to check the edge cases
                //
                // Previously, a direct lookup to check the mass of the current fragment was
                // performed, but the pointer indirection + float comparison can slow down
                // open searches by as much as 2x!!
                // e.g. used to be `self.db[frag.peptide_index].monoisotopic >= precursor_lo`
                (frag.peptide_index.0 > self.pre_idx_lo as u32
                    || (frag.peptide_index.0 == self.pre_idx_lo as u32
                        && self.db[frag.peptide_index].monoisotopic >= precursor_lo))
                    && (frag.peptide_index.0 < self.pre_idx_hi as u32
                        || (frag.peptide_index.0 == self.pre_idx_hi as u32
                            && self.db[frag.peptide_index].monoisotopic <= precursor_hi))
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
/// * `slice[right] > high || right == slice.len()`
/// * `0 <= left <= right <= slice.len()`
#[inline]
pub fn binary_search_slice<T, F, S>(slice: &[T], key: F, low: S, high: S) -> (usize, usize)
where
    F: Fn(&T, &S) -> Ordering,
{
    let left_idx = slice
        .partition_point(|a| key(a, &low) == Ordering::Less)
        .saturating_sub(1);

    let right_idx =
        slice[left_idx..].partition_point(|a| key(a, &high) != Ordering::Greater) + left_idx;

    (left_idx, right_idx)
}

#[cfg(test)]
mod test {
    use std::sync::Arc;

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

    #[test]
    fn structured_variable_mod_config_round_trips() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "variable_mods": {
                "M": [15.9949],
                "K": [
                    {"mass": 42.0106, "max_count": 1},
                    {"mass": 14.0157}
                ]
            },
            "max_variable_mods": 2,
            "max_combinations": 0
        }))
        .unwrap();

        let params = builder.make_parameters().unwrap();
        assert_eq!(params.max_variable_mods, 2);
        assert_eq!(params.max_combinations, Some(1));

        let resolved = &params.modifications;
        let mods = &resolved.variable_mods;
        let mass = |id| resolved.registry.mass(Some(id));
        assert_eq!(mods.len(), 3);
        assert_eq!(mods[0].0, ModificationSpecificity::Residue(b'K'));
        assert!((mass(mods[0].1) - 42.0106).abs() < 1e-4);
        assert_eq!(mods[0].2, Some(1));
        assert_eq!(mods[1].0, ModificationSpecificity::Residue(b'K'));
        assert!((mass(mods[1].1) - 14.0157).abs() < 1e-4);
        assert_eq!(mods[1].2, None);
        assert_eq!(mods[2].0, ModificationSpecificity::Residue(b'M'));
        assert!((mass(mods[2].1) - 15.9949).abs() < 1e-4);
        assert_eq!(mods[2].2, None);

        let serialized = serde_json::to_value(params).unwrap();
        let k_entries = &serialized["variable_mods"]["K"];
        assert!(k_entries[0].is_object());
        assert_eq!(k_entries[0]["max_count"], 1);
        // `{"mass": 14.0157}` carries nothing beyond the mass, so it is written back
        // in the bare-mass form
        assert_eq!(k_entries[1], serde_json::json!(14.0157f32));
        assert!(serialized["variable_mods"]["M"][0].is_number());
    }

    #[test]
    fn digestion() {
        let fasta = r#"
        >sp|AAAAA
        MEWKLEQSMREQALLKAQLTQLK
        >sp|BBBBB
        RMEWKLEQSMREQALLKAQLTQLK
        "#;

        let fasta = Fasta::parse(fasta.into(), "rev_", false);

        // Make sure that FASTA parsed OK
        assert_eq!(
            fasta.targets,
            vec![
                (
                    Arc::from("sp|AAAAA".to_string()),
                    "MEWKLEQSMREQALLKAQLTQLK".into()
                ),
                (
                    Arc::from("sp|BBBBB".to_string()),
                    "RMEWKLEQSMREQALLKAQLTQLK".into()
                ),
            ]
        );

        let variable_mods = [(
            ModificationSpecificity::ProteinN(None),
            vec![ModEntry::Mass(42.0)],
        )]
        .into_iter()
        .collect();
        let params = Parameters {
            bucket_size: 128,
            enzyme: EnzymeBuilder {
                missed_cleavages: Some(1),
                min_len: Some(6),
                max_len: Some(10),
                ..Default::default()
            },
            peptide_min_mass: 150.0,
            peptide_max_mass: 5000.0,
            ion_kinds: vec![Kind::B, Kind::Y],
            min_ion_index: 2,
            modifications: ResolvedModifications::new(&HashMap::default(), &variable_mods).unwrap(),
            max_variable_mods: 2,
            max_combinations: None,
            decoy_tag: "rev_".into(),
            generate_decoys: false,
            fasta: "none".into(),
            prefilter: false,
            prefilter_chunk_size: 0,
            prefilter_low_memory: true,
        };

        let peptides = params.digest(&fasta);

        let expected = [
            "EQALLK",
            "LEQSMR",
            "AQLTQLK",
            "MEWKLEQSMR",
            "[+42]-MEWKLEQSMR",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

        let mods = params.modifications.registry.clone();
        let sequences = peptides
            .iter()
            .map(|p| p.display(&mods).to_string())
            .collect::<Vec<_>>();
        assert_eq!(expected, sequences);

        // All peptides are shared except for the protein N-term mod
        for peptide in &peptides[..4] {
            assert_eq!(peptide.proteins.len(), 2, "{:?}", peptide);
        }
        // Ensure that this mod is uniquely called as the first protein
        assert_eq!(
            peptides.last().unwrap().proteins,
            vec!["sp|AAAAA".to_string().into()]
        );
    }

    #[test]
    fn named_modifications() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": {
                "C": {"mass": 57.0215, "name": "UNIMOD:4"},
                "^": 229.1629
            },
            "variable_mods": {
                "M": [{"mass": 15.9949, "name": "UNIMOD:35"}],
                "[": [{"mass": 42.0106, "name": "Acetyl"}]
            },
            "generate_decoys": false,
            "peptide_min_mass": 0.0
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();

        // Registry ids are independent of HashMap iteration order
        let resolved = &params.modifications;
        assert_eq!(resolved.registry.len(), 4);

        // Serializes back to the config it was built from
        let serialized = serde_json::to_value(&params).unwrap();
        assert_eq!(
            serialized["static_mods"],
            serde_json::json!({"C": {"mass": 57.0215f32, "name": "UNIMOD:4"}, "^": 229.1629f32})
        );
        assert_eq!(
            serialized["variable_mods"],
            serde_json::json!({
                "M": [{"mass": 15.9949f32, "name": "UNIMOD:35"}],
                "[": [{"mass": 42.0106f32, "name": "Acetyl"}]
            })
        );
        let reparsed: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": serialized["static_mods"],
            "variable_mods": serialized["variable_mods"],
            "generate_decoys": false,
            "peptide_min_mass": 0.0
        }))
        .unwrap();
        let reparsed = reparsed.make_parameters().unwrap().modifications;
        assert_eq!(reparsed.registry, resolved.registry);
        assert_eq!(reparsed.static_mods, resolved.static_mods);
        assert_eq!(reparsed.variable_mods, resolved.variable_mods);

        let fasta = Fasta::parse(">sp|AAAAA\nMCPEPK".into(), "rev_", false);
        let db = params.build(fasta);
        let mut sequences = db
            .peptides
            .iter()
            .map(|p| p.display(&db.modifications).to_string())
            .collect::<Vec<_>>();
        sequences.sort();
        assert_eq!(
            sequences,
            vec![
                "[+229.1629]-MC[UNIMOD:4]PEPK",
                "[+229.1629]-M[UNIMOD:35]C[UNIMOD:4]PEPK",
                "[Acetyl]-MC[UNIMOD:4]PEPK",
                "[Acetyl]-M[UNIMOD:35]C[UNIMOD:4]PEPK",
            ]
        );
    }

    #[test]
    fn invalid_modification_config() {
        let parse = |mods: serde_json::Value| {
            let mut config = serde_json::json!({"fasta": "none"});
            config
                .as_object_mut()
                .unwrap()
                .extend(mods.as_object().unwrap().clone());
            serde_json::from_value::<Builder>(config)
        };

        // Rejected while parsing
        assert!(parse(
            serde_json::json!({"variable_mods": {"M": [{"mass": 15.9949, "name": "Ox]"}]}})
        )
        .is_err());
        assert!(parse(serde_json::json!({"variable_mods": {"Z": [15.9949]}})).is_err());
        assert!(parse(serde_json::json!({"static_mods": {"^Z": 42.0}})).is_err());

        // Rejected when building parameters
        let builder =
            parse(serde_json::json!({"static_mods": {"C": {"mass": 57.0215, "max_count": 1}}}));
        assert_eq!(
            builder.unwrap().make_parameters().err(),
            Some(ModificationError::StaticMaxCount(
                ModificationSpecificity::Residue(b'C')
            ))
        );

        let variable_mods = (0..256)
            .map(|i| (i as f32).into())
            .collect::<Vec<serde_json::Value>>();
        let builder = parse(serde_json::json!({"variable_mods": {"M": variable_mods}}));
        assert_eq!(
            builder.unwrap().make_parameters().err(),
            Some(ModificationError::TooMany)
        );
    }
}
