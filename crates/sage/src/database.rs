use crate::enzyme::{Enzyme, EnzymeParameters};
use crate::fasta::Fasta;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{validate_mods, validate_var_mods, ModificationSpecificity};
use crate::peptide::Peptide;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::HashMap;
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
    pub restrict: Option<char>,
    pub c_terminal: Option<bool>,
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            cleave_at: Some("KR".into()),
            restrict: Some('P'),
            c_terminal: Some(true),
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        EnzymeParameters {
            missed_cleavages: en.missed_cleavages.unwrap_or(1),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enyzme: Enzyme::new(
                &en.cleave_at.unwrap_or_else(|| "KR".into()),
                en.restrict,
                en.c_terminal.unwrap_or(true),
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
    /// Which kind of fragment ions to generate (a, b, c, x, y, z)
    pub ion_kinds: Option<Vec<Kind>>,
    /// Minimum ion index to be generated: 1 will remove b1/y1 ions
    /// 2 will remove b1/b2/y1/y2 ions, etc
    pub min_ion_index: Option<usize>,
    /// Static modifications to add to matching amino acids
    pub static_mods: Option<HashMap<String, f32>>,
    /// Variable modifications to add to matching amino acids
    pub variable_mods: Option<HashMap<String, crate::modification::ValueOrVec>>,
    /// Limit number of variable modifications on a peptide
    pub max_variable_mods: Option<usize>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
}

impl Builder {
    pub fn make_parameters(self) -> Parameters {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        Parameters {
            bucket_size,
            fragment_min_mz: self.fragment_min_mz.unwrap_or(150.0),
            fragment_max_mz: self.fragment_max_mz.unwrap_or(2000.0),
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            static_mods: validate_mods(self.static_mods),
            variable_mods: validate_var_mods(self.variable_mods),
            max_variable_mods: self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2),
            generate_decoys: self.generate_decoys.unwrap_or(true),
            fasta: self.fasta.expect("A fasta file must be provided!"),
        }
    }

    pub fn update_fasta(&mut self, fasta: String) {
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
    ion_kinds: Vec<Kind>,
    min_ion_index: usize,
    static_mods: HashMap<ModificationSpecificity, f32>,
    variable_mods: HashMap<ModificationSpecificity, Vec<f32>>,
    max_variable_mods: usize,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub fasta: String,
}

impl Parameters {
    fn digest(&self, fasta: &Fasta) -> Vec<Peptide> {
        log::trace!("digesting fasta");
        let enzyme = self.enzyme.clone().into();
        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        let digests = fasta.digest(&enzyme);

        let mods = self
            .variable_mods
            .iter()
            .flat_map(|(a, b)| b.iter().map(|b| (*a, *b)))
            .collect::<Vec<_>>();

        let targets: DashSet<_, FnvBuildHasher> = DashSet::default();
        digests
            .par_iter()
            .filter(|digest| !digest.decoy)
            .for_each(|digest| {
                targets.insert(digest.sequence.clone().into_bytes());
            });

        log::trace!("modifying peptides");
        let mut target_decoys = digests
            .into_par_iter()
            .map(Peptide::try_from)
            .filter_map(Result::ok)
            .flat_map_iter(|peptide| {
                peptide
                    .apply(&mods, &self.static_mods, self.max_variable_mods)
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

        log::trace!("sorting and deduplicating peptides");

        // This is equivalent to a stable sort
        target_decoys.par_sort_unstable_by(|a, b| {
            a.monoisotopic
                .total_cmp(&b.monoisotopic)
                .then_with(|| a.initial_sort(b))
        });
        target_decoys.dedup_by(|remove, keep| {
            if remove.sequence == keep.sequence
                && remove.modifications == keep.modifications
                && remove.nterm == keep.nterm
                && remove.cterm == keep.cterm
            {
                keep.proteins.extend(remove.proteins.iter().cloned());
                true
            } else {
                false
            }
        });

        target_decoys
            .par_iter_mut()
            .for_each(|peptide| peptide.proteins.sort_unstable());

        target_decoys
    }

    // pub fn build(self) -> Result<IndexedDatabase, Box<dyn std::error::Error + Send + Sync + 'static>> {
    pub fn build(self, fasta: Fasta) -> IndexedDatabase {
        let target_decoys = self.digest(&fasta);
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
                    .flat_map(|kind| IonSeries::new(peptide, *kind).enumerate())
                    .filter(|(ion_idx, ion)| {
                        // Don't store b1, b2, y1, y2 ions for preliminary scoring
                        let ion_idx_filter = match ion.kind {
                            Kind::A | Kind::B | Kind::C => (ion_idx + 1) > self.min_ion_index,
                            Kind::X | Kind::Y | Kind::Z => {
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

        let potential_mods = self
            .variable_mods
            .iter()
            .flat_map(|(a, b)| b.iter().map(|b| (*a, *b)))
            .collect::<Vec<(ModificationSpecificity, f32)>>();

        IndexedDatabase {
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            potential_mods,
            decoy_tag: self.decoy_tag,
        }
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
    pub ion_kinds: Vec<Kind>,
    pub(crate) min_value: Vec<f32>,
    /// Keep a list of potential (AA, mass) modifications for RT prediction
    pub potential_mods: Vec<(ModificationSpecificity, f32)>,
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
                fragment,
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

            // let (inner_left, inner_right) = binary_search_slice(
            //     slice,
            //     |frag, bounds| (frag.peptide_index.0 as usize).cmp(bounds),
            //     self.pre_idx_lo,
            //     self.pre_idx_hi,
            // );
            let (inner_left, inner_right) = interpolation_search_slice(
                slice,
                PeptideIx(self.pre_idx_lo as u32),
                PeptideIx(self.pre_idx_hi as u32),
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

fn interpolation_search(slice: &[Theoretical], peptide: PeptideIx) -> Result<usize, usize> {
    let mut low = 0;
    let mut high = slice.len().saturating_sub(1);

    if peptide > slice[slice.len() - 1].peptide_index {
        return Err(slice.len());
    } else if peptide < slice[0].peptide_index {
        return Err(0);
    }

    while low < high {
        // SAFETY: We know that `low` and `high` are always in bounds
        // and that we are not dividing by 0
        let slope = unsafe {
            (slice.get_unchecked(high).peptide_index.0 - slice.get_unchecked(low).peptide_index.0)
                as f32
                / (high - low) as f32
        };

        let diff = ((peptide.0 as i64 - slice[low].peptide_index.0 as i64) as f32 / slope) as isize;
        let mid = low as isize + diff;
        let mid = mid.clamp(low as isize, high as isize) as usize;

        // SAFETY: We know that `low`  and `high` are always in bounds
        // and that `mid` is clamped between these values.
        let cmp = unsafe { slice.get_unchecked(mid).peptide_index.cmp(&peptide) };

        if cmp == Ordering::Less {
            low = mid + 1;
        } else if cmp == Ordering::Greater {
            high = mid - 1;
        } else {
            return Ok(mid);
        }
    }
    return Err(low);
}

#[inline]
pub fn interpolation_search_slice(
    slice: &[Theoretical],
    low: PeptideIx,
    high: PeptideIx,
) -> (usize, usize) {
    // let left_idx = match slice.binary_search_by(|a| key(a, &low)) {
    let left_idx = match interpolation_search(slice, low) {
        Ok(idx) | Err(idx) => {
            let mut idx = idx.saturating_sub(4);
            while idx > 4 && slice[idx].peptide_index.cmp(&low) != Ordering::Less {
                idx -= 4;
            }
            idx
        }
    };

    // let right_idx = match slice[left_idx..].binary_search_by(|a| key(a, &high)) {
    let right_idx = match interpolation_search(&slice[left_idx..], high) {
        Ok(idx) | Err(idx) => {
            let mut idx = idx + left_idx;
            while idx < slice.len() - 4 && slice[idx].peptide_index.cmp(&high) != Ordering::Greater
            {
                idx = idx.saturating_add(4);
            }
            idx.min(slice.len())
        }
    };
    (left_idx, right_idx)
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
                    Arc::new("sp|AAAAA".to_string()),
                    "MEWKLEQSMREQALLKAQLTQLK".into()
                ),
                (
                    Arc::new("sp|BBBBB".to_string()),
                    "RMEWKLEQSMREQALLKAQLTQLK".into()
                ),
            ]
        );

        let params = Parameters {
            bucket_size: 128,
            enzyme: EnzymeBuilder {
                missed_cleavages: Some(1),
                min_len: Some(6),
                max_len: Some(10),
                ..Default::default()
            },
            fragment_min_mz: 100.0,
            fragment_max_mz: 1000.0,
            peptide_min_mass: 150.0,
            peptide_max_mass: 5000.0,
            ion_kinds: vec![Kind::B, Kind::Y],
            min_ion_index: 2,
            static_mods: HashMap::default(),
            variable_mods: [(ModificationSpecificity::ProteinN(None), vec![42.0])]
                .into_iter()
                .collect(),
            max_variable_mods: 2,
            decoy_tag: "rev_".into(),
            generate_decoys: false,
            fasta: "none".into(),
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

        let sequences = peptides.iter().map(|p| p.to_string()).collect::<Vec<_>>();
        assert_eq!(expected, sequences);

        // All peptides are shared except for the protein N-term mod
        for peptide in &peptides[..4] {
            assert_eq!(peptide.proteins.len(), 2);
        }
        // Ensure that this mod is uniquely called as the first protein
        assert_eq!(
            peptides.last().unwrap().proteins,
            vec!["sp|AAAAA".to_string().into()]
        );
    }
}
