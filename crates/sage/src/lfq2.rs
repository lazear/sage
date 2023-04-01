use crate::database::{binary_search_slice, IndexedDatabase, PeptideIx};
use crate::mass::{Tolerance, NEUTRON};
use crate::ml::{matrix::Matrix, retention_alignment::Alignment};
use crate::scoring::Feature;
use crate::spectrum::ProcessedSpectrum;
use dashmap::DashMap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cell::Cell;
use std::{collections::HashMap, hash::BuildHasherDefault};

/// Minimum normalized spectral angle required to integrate a peak
// const MIN_SPECTRAL_ANGLE: f64 = 0.70;
/// Retention time tolerance, in fraction of total run length, to search for
/// precursor ions
const RT_TOL: f32 = 0.0050;
/// Width of gaussian kernel used for smoothing intensities
const K_WIDTH: usize = 10;
/// Mass tolerance, in ppm, to seach for precursor ions
// const PPM_TOL: f32 = 5.0;
/// Number of equally spaced bins that will be used to integrate ions in (-RT_TOL, +RT_TOL)
const GRID_SIZE: usize = 100;
/// Number of isotopes to search for
const N_ISOTOPES: usize = 4;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum PeakScoringStrategy {
    RetentionTime,
    SpectralAngle,
    Hybrid,
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum IntegrationStrategy {
    Max,
    Sum,
}
#[derive(Copy, Clone, Debug, Deserialize, Serialize)]
pub struct LfqSettings {
    pub peak_scoring: PeakScoringStrategy,
    pub integration: IntegrationStrategy,
    pub spectral_angle: f64,
    pub ppm_tolerance: f32,
}

impl Default for LfqSettings {
    fn default() -> Self {
        Self {
            peak_scoring: PeakScoringStrategy::Hybrid,
            integration: IntegrationStrategy::Sum,
            spectral_angle: 0.70,
            ppm_tolerance: 5.0,
        }
    }
}

#[derive(Copy, Clone, Debug, Serialize)]
pub struct PrecursorRange {
    pub rt: f32,
    pub mass_lo: f32,
    pub mass_hi: f32,
    pub charge: u8,
    pub isotope: usize,
    pub peptide: PeptideIx,
    pub file_id: usize,
}

/// Create a data structure analogous to [`IndexedDatabase`] - instaed of
/// storing fragment masses binned by precursor mass, store MS1 precursors
/// binned by RT - This should enable rapid quantification as well
pub struct FeatureMap {
    ranges: Vec<PrecursorRange>,
    min_rts: Vec<f32>,
    bin_size: usize,
    settings: LfqSettings,
}

pub fn build_feature_map(settings: LfqSettings, features: &[Feature]) -> FeatureMap {
    let map: DashMap<PeptideIx, PrecursorRange, fnv::FnvBuildHasher> = DashMap::default();
    features
        .into_par_iter()
        .filter(|feat| feat.spectrum_q <= 0.01 && feat.peptide_q <= 0.05)
        .for_each(|feat| {
            // `features` is sorted by confidence, so just take the first entry
            if !map.contains_key(&feat.peptide_idx) {
                map.insert(
                    feat.peptide_idx,
                    PrecursorRange {
                        rt: feat.aligned_rt,
                        mass_lo: feat.calcmass,
                        mass_hi: 0.0,
                        peptide: feat.peptide_idx,
                        charge: feat.charge,
                        isotope: 0,
                        file_id: feat.file_id,
                    },
                );
            }
        });

    // Cartesian product of (observed RT ranges, mass) with charge and isotopes
    let mut ranges = map
        .into_par_iter()
        .flat_map_iter(|(_, range)| {
            (2..5).flat_map(move |charge| {
                (0..N_ISOTOPES).map(move |isotope| {
                    let mass = (range.mass_lo + isotope as f32 * NEUTRON) / charge as f32;
                    let (mass_lo, mass_hi) =
                        Tolerance::Ppm(-settings.ppm_tolerance, settings.ppm_tolerance)
                            .bounds(mass);

                    PrecursorRange {
                        mass_lo,
                        mass_hi,
                        charge,
                        isotope,
                        ..range
                    }
                })
            })
        })
        .collect::<Vec<_>>();

    // Essentially the same procedure for binning as the MS2 search engine
    // (see `database.rs` for explanation)
    ranges.par_sort_unstable_by(|a, b| a.rt.total_cmp(&b.rt));
    let min_rts = ranges
        .par_chunks_mut(16 * 1024)
        .map(|chunk| {
            // There should always be at least one item in the chunk!
            //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
            let min = chunk[0].rt;
            chunk.par_sort_unstable_by(|a, b| a.mass_lo.total_cmp(&b.mass_lo));
            min
        })
        .collect::<Vec<_>>();

    FeatureMap {
        ranges,
        min_rts,
        bin_size: 16 * 1024,
        settings,
    }
}

struct Query<'a> {
    ranges: &'a [PrecursorRange],
    page_lo: usize,
    page_hi: usize,
    bin_size: usize,
    min_rt: f32,
    max_rt: f32,
}

impl FeatureMap {
    fn rt_slice(&self, rt: f32, rt_tol: f32) -> Query<'_> {
        let (page_lo, page_hi) = binary_search_slice(
            &self.min_rts,
            |rt, x| rt.total_cmp(x),
            rt - rt_tol,
            rt + rt_tol,
        );

        Query {
            ranges: &self.ranges,
            page_lo,
            page_hi,
            bin_size: self.bin_size,
            max_rt: rt + rt_tol,
            min_rt: rt - rt_tol,
        }
    }

    /// Run label-free quantification module
    pub fn quantify(
        &self,
        db: &IndexedDatabase,
        spectra: &[ProcessedSpectrum],
        alignments: &[Alignment],
    ) -> HashMap<PeptideIx, Vec<f64>> {
        let scores: DashMap<PeptideIx, Grid, BuildHasherDefault<fnv::FnvHasher>> =
            DashMap::default();

        log::info!("tracing MS1 features");
        spectra
            .par_iter()
            .filter(|s| s.level == 1)
            .for_each(|spectrum| {
                let a = alignments[spectrum.file_id];
                let rt = (spectrum.scan_start_time / a.max_rt) * a.slope + a.intercept;
                let query = self.rt_slice(rt, RT_TOL);

                for peak in &spectrum.peaks {
                    for entry in query.mass_lookup(peak.mass) {
                        let mut grid = scores.entry(entry.peptide).or_insert_with(|| {
                            let p = &db[entry.peptide];
                            let dist = crate::isotopes::peptide_isotopes(p.carbons, p.sulfurs);
                            Grid::new(entry, RT_TOL, dist, alignments.len(), GRID_SIZE)
                        });

                        grid.add_entry(rt, entry.isotope, spectrum.file_id, peak.intensity);
                    }
                }
            });

        log::info!("integrating MS1 features");

        scores
            .into_par_iter()
            .filter_map(|(peptide_ix, mut grid)| {
                // MS1 ions have been added to any relevant grids, so we now
                // attempt to trace the peaks, find the best peak, and integrate
                // it across all of the files
                let mut traces = grid.summarize_traces();
                let data = traces.integrate(&self.settings)?;

                #[cfg(debug)]
                {
                    let writer =
                        std::fs::File::create(format!("dotp/{}.json", db[peptide_ix])).unwrap();
                    let writer = std::io::BufWriter::new(writer);
                    serde_json::to_writer_pretty(
                        writer,
                        &serde_json::json!({
                            "dist": &grid.distribution,
                            "data": &data,
                            "cosine": &traces.spectral_angle.data,
                            "dotp": &traces.dot_product.data,
                            "left": traces.left,
                            "right": traces.right,
                        }),
                    )
                    .unwrap();
                }

                Some((peptide_ix, data))
            })
            .collect::<HashMap<PeptideIx, _>>()
    }
}

pub struct Grid {
    rt_min: f32,
    rt_step: f32,
    files: usize,
    reference_file_id: usize,
    /// Relative theoretical isotopic abundances
    distribution: [f32; N_ISOTOPES],
    /// Matrix of summed intensities for each isotopic trace in each file, divided
    /// among equally spaced retention time bins. This is a [N_FILES * N_ISOTOPES, GRID_SIZE]
    /// sized matrix.
    ///
    /// Isotopic summed intensities are arranged in consequtive rows, ordered by file. The first
    /// N_ISOTOPE rows correspond to file 0, then the following N_ISOTOPE rows correspond to file 1, etc.
    /// Indexing into the matrix is done by [(n_file * N_ISOTOPES + isotope, rt_window)]
    pub matrix: Matrix,
}

pub struct Traces {
    /// Matrix of dot(MS1 ions, Grid.distribution). This collapses our N_FILES * N_ISOTOPES rows
    /// down to just N_FILES
    dot_product: Matrix,
    /// Matrix of spectral angles at each retention time for each file
    spectral_angle: Matrix,
    // These are just here for debugging purposes
    pub left: usize,
    pub right: usize,
    pub best_peak: Option<Peak>,
    reference_file_id: usize,
}

#[derive(Clone, Debug)]
pub struct Peak {
    /// Discretized retention time
    pub rt: usize,
    /// Intensity weighted normalized spectral angle
    pub spectral_angle: f64,
    // Index of parent cell in [`DisjointPeakSet`], this is essentially an intrusively
    // linked list, allowing us to efficienty merge local maxima
    parent: Cell<usize>,
}

impl Peak {
    /// Higher is better
    pub fn score(&self, strategy: PeakScoringStrategy) -> f64 {
        const CENTER: isize = GRID_SIZE as isize / 2;
        match strategy {
            PeakScoringStrategy::RetentionTime => {
                1.0 - ((self.rt as isize - CENTER).abs() as f64 / CENTER as f64)
            }
            PeakScoringStrategy::SpectralAngle => self.spectral_angle,
            PeakScoringStrategy::Hybrid => {
                let rt_score = self.score(PeakScoringStrategy::RetentionTime);
                (self.spectral_angle.powi(2) * rt_score.powi(2)).sqrt()
            }
        }
    }
}

#[derive(Default)]
struct DisjointPeakSet {
    peaks: Vec<Peak>,
}

impl DisjointPeakSet {
    /// Find the best peak in the set that has not been merged into a parent peak
    pub fn best_peak(self, strategy: PeakScoringStrategy) -> Option<Peak> {
        let mut best = None;
        for (idx, peak) in self.peaks.into_iter().enumerate() {
            if idx == peak.parent.get() {
                let x = best
                    .as_ref()
                    .map(|peak: &Peak| peak.score(strategy))
                    .unwrap_or(0.0);
                if peak.score(strategy) > x {
                    best = Some(peak);
                }
            }
        }
        best
    }

    /// Merge the worse scoring peak into the better scoring peak
    pub fn union(&mut self, left: usize, right: usize, strategy: PeakScoringStrategy) {
        let a = self.find_set(left);
        let b = self.find_set(right);
        if self.peak(left).score(strategy) > self.peak(right).score(strategy) {
            self.peaks[b].parent.set(a);
        } else {
            self.peaks[a].parent.set(b);
        }
    }

    /// Record a new local maxima, retuning it's index in the set
    pub fn singleton(&mut self, rt: usize, spectral_angle: f64) -> usize {
        let idx = self.peaks.len();
        self.peaks.push(Peak {
            rt,
            spectral_angle,
            parent: Cell::new(idx),
        });
        idx
    }

    pub fn peak(&self, idx: usize) -> &Peak {
        let ptr = self.find_set(idx);
        &self.peaks[ptr]
    }

    pub fn find_set(&self, idx: usize) -> usize {
        if idx != self.peaks[idx].parent.get() {
            self.peaks[idx]
                .parent
                .set(self.find_set(self.peaks[idx].parent.get()));
        }
        self.peaks[idx].parent.get()
    }
}

impl Traces {
    /// Calculate and apply time warping factors
    fn warp(&mut self) {
        let time_warps = self.find_time_warps(&self.dot_product, 75);
        Self::apply_time_warps(&mut self.spectral_angle, &time_warps);
        Self::apply_time_warps(&mut self.dot_product, &time_warps);
    }

    /// Find time warping offsets for each file that maximize the dot product
    /// with the most intense run
    ///
    /// * Use the LC-MS run with the most confident PSM for a peptide as the reference run
    /// * For each LC-MS run, find the time warping shif that maximizes dot product
    ///   with the `reference` run
    fn find_time_warps(&self, matrix: &Matrix, slack: isize) -> Vec<isize> {
        let reference = matrix.row_slice(self.reference_file_id);
        let mut offsets = vec![0; matrix.rows];

        for (row, offset) in offsets.iter_mut().enumerate() {
            let run = matrix.row_slice(row);

            let mut best_offset = (0, 0.0);
            for offset in -slack..=slack {
                let mut dot = 0.0;
                for (i, ref_int) in reference.iter().enumerate() {
                    let j = i as isize + offset;
                    if j >= 0 && j < run.len() as isize {
                        dot += ref_int * run[j as usize];
                    }
                }

                if dot >= best_offset.1 {
                    best_offset = (offset, dot);
                }
            }
            *offset = best_offset.0;
        }
        offsets
    }

    /// Perform local Correlation Optimization Warping
    fn apply_time_warps(matrix: &mut Matrix, time_warps: &[isize]) {
        for (row, warp) in time_warps.iter().enumerate() {
            let run = matrix.row_slice_mut(row);
            let mut shifted = vec![0.0; run.len()];
            for (i, val) in shifted.iter_mut().enumerate() {
                let j = i as isize + warp;
                if j >= 0 && j < run.len() as isize {
                    *val = run[j as usize];
                }
            }
            run.copy_from_slice(&shifted);
        }
    }

    fn find_peak(
        spectral_angle: &[f64],
        threshold: f64,
        strategy: PeakScoringStrategy,
    ) -> Option<Peak> {
        let mut set = DisjointPeakSet::default();

        // Locate points in RT space that are local maxima in both
        // spectral angle and intensity dimensions
        for index in 1..spectral_angle.len().saturating_sub(1) {
            let sa = spectral_angle[index];
            if sa >= threshold && sa >= spectral_angle[index - 1] && sa >= spectral_angle[index + 1]
            {
                set.singleton(index, sa);
            }
        }

        // Proceed from left to right, attempting to merge peaks if there is no valley between them
        'outer: for idx in 1..set.peaks.len() {
            let left = set.peak(idx - 1);
            let right = set.peak(idx);

            if right.rt - left.rt > GRID_SIZE / 4 {
                continue;
            }

            let sa_valley = left.spectral_angle.max(right.spectral_angle) * 0.75;

            for sa in &spectral_angle[left.rt..right.rt] {
                if *sa <= sa_valley {
                    continue 'outer;
                }
            }

            set.union(idx - 1, idx, strategy);
        }

        set.best_peak(strategy)
    }

    /// Align and integrate MS1 traces across files
    ///
    /// * Calculate time warping factors for each file, aligning them so that
    ///   the correlation between them is maximized (we just do a dot product)
    /// * Locate the retention time window corresponding to the maximum average
    ///   angle observed across all of the files
    /// * Integrate all of the MS1 traces within said window, returning a vector
    ///   of length `n_files` containing the summed MS1 intensities
    pub fn integrate(&mut self, settings: &LfqSettings) -> Option<Vec<f64>> {
        self.warp();

        // Calculate the average spectral angle at each point in the retention time grid
        // weighted by intensities
        let mut score = Vec::with_capacity(self.spectral_angle.cols);
        for col in 0..self.spectral_angle.cols {
            let mut summed_int = 0.0;
            let mut weighted = 0.0;
            for (sa, dotp) in self.spectral_angle.col(col).zip(self.dot_product.col(col)) {
                weighted += sa * dotp;
                summed_int += dotp;
            }
            score.push(weighted / summed_int);
        }

        // TODO: Instead of this fallback, maybe try determining which individual file contains the best trace,
        // and using that file as the reference run for peak selection
        let best = Traces::find_peak(&score, settings.spectral_angle, settings.peak_scoring)?;

        // Find peak boundaries
        let mut left = best.rt - 1;
        let mut right = best.rt + 1;

        let threshold = best.spectral_angle * 0.75;

        // Don't let peaks extend more than 20 RT bins to either side
        while left > best.rt.saturating_sub(20) && score[left] >= threshold {
            left -= 1;
        }

        while right < score.len().saturating_sub(1).min(best.rt + 20) && score[right] >= threshold {
            right += 1;
        }

        self.left = left;
        self.right = right;
        self.best_peak = Some(best);

        // Actually perform integration
        let mut areas = Vec::with_capacity(self.dot_product.rows);

        for file in 0..self.dot_product.rows {
            let area = match settings.integration {
                IntegrationStrategy::Sum => self.dot_product.row_slice(file)[left..right]
                    .iter()
                    .sum::<f64>(),
                IntegrationStrategy::Max => self.dot_product.row_slice(file)[left..right]
                    .iter()
                    .fold(0.0f64, |acc, x| acc.max(*x)),
            };

            areas.push(area);
        }

        Some(areas)
    }
}

impl Grid {
    pub fn new(
        entry: &PrecursorRange,
        rt_tol: f32,
        distribution: [f32; 4],
        files: usize,
        grid_size: usize,
    ) -> Grid {
        let matrix = Matrix::new(vec![0.0; grid_size * files * 4], files * 4, grid_size);
        let rt_step = (rt_tol * 2.0) / (grid_size) as f32;

        Grid {
            rt_min: entry.rt - rt_tol,
            rt_step,
            distribution,
            matrix,
            files,
            reference_file_id: entry.file_id,
        }
    }

    pub fn add_entry(&mut self, spectrum_rt: f32, isotope: usize, file_id: usize, intensity: f32) {
        let bin_lo = ((spectrum_rt - self.rt_min) / self.rt_step).floor() as usize;
        let bin_lo = bin_lo.min(self.matrix.cols - 1);
        let bin_hi = (bin_lo + 1).min(self.matrix.cols - 1);

        let bin_lo_rt = bin_lo as f32 * self.rt_step + self.rt_min;
        // what fraction [0.0, 1.0] of the way are we to the higher bin?
        let interp = (spectrum_rt - bin_lo_rt) / self.rt_step;

        self.matrix[(file_id * N_ISOTOPES + isotope, bin_lo)] +=
            ((1.0 - interp) * intensity) as f64;
        self.matrix[(file_id * N_ISOTOPES + isotope, bin_hi)] += (interp * intensity) as f64;
    }

    /// Combine individual isotopic traces across files into aligned, summed
    /// MS1 traces for each file
    ///
    /// * Perform gaussian smoothing on summed intensities
    /// * Calculate normalized spectral angle for observed isotopic distribution
    ///   relative to theoretical distribution
    pub fn summarize_traces(&mut self) -> Traces {
        let k = gaussian_kernel(0.5, K_WIDTH);

        let mut spectral_angle = Matrix::new(
            vec![0.0; self.files * self.matrix.cols],
            self.files,
            self.matrix.cols,
        );

        let mut dot_product = spectral_angle.clone();

        // square root of the summed squared relative abundances of theoretical
        // isotopic distribution
        let ss_dist = self
            .distribution
            .iter()
            .map(|x| x.powi(2))
            .sum::<f32>()
            .sqrt() as f64;

        for file in 0..self.files {
            let mut summed_squared_intensities = vec![0.0; self.matrix.cols];
            for isotope in 0..N_ISOTOPES {
                let convolved = convolve(self.matrix.row_slice(file * N_ISOTOPES + isotope), &k);
                for (col, intensity) in convolved.iter().enumerate() {
                    spectral_angle[(file, col)] += intensity * self.distribution[isotope] as f64;
                    summed_squared_intensities[col] += intensity.powi(2);
                }
                self.matrix
                    .row_slice_mut(file * N_ISOTOPES + isotope)
                    .copy_from_slice(&convolved);
            }

            for (col, ss) in summed_squared_intensities.iter().enumerate() {
                let dot = spectral_angle[(file, col)];
                let similarity = dot / (ss.sqrt() * ss_dist);

                // Calculate the normalized spectral angle
                spectral_angle[(file, col)] = 1.0 - 2.0 * similarity.acos() / std::f64::consts::PI;
                dot_product[(file, col)] = dot;
            }
        }

        Traces {
            dot_product,
            spectral_angle,
            left: 0,
            right: self.matrix.cols - 1,
            reference_file_id: self.reference_file_id,
            best_peak: None,
        }
    }
}

/// Create a symmetrical gaussian kernel of given standard deviation and length
fn gaussian_kernel(sigma: f64, len: usize) -> Vec<f64> {
    let step = 2.0 / (len - 1) as f64;
    let constant = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());

    let mut kernel = (0..len)
        .map(|i| {
            let x = i as f64 * step - 1.0;
            constant * (-0.5 * (x / sigma).powi(2)).exp()
        })
        .collect::<Vec<_>>();

    let sum = kernel.iter().sum::<f64>();
    kernel.iter_mut().for_each(|x| *x /= sum);
    kernel
}

/// Convolve a signal with a symmetrical kernel
/// - This should behave the same as `np.convolve(..., mode='same')`
fn convolve(slice: &[f64], kernel: &[f64]) -> Vec<f64> {
    // Middle index of kernel
    let n = kernel.len() - (kernel.len() / 2);

    (0..slice.len())
        .map(|idx| {
            // If idx < kernel.len(), take only some subset of the kernel (at least half)
            let k = &kernel[kernel.len().saturating_sub(n + idx)..];
            // If idx < kernel.len(), then start from 0
            let w = &slice[idx.saturating_sub(n - 1)..];
            // Dot product
            w.iter().zip(k).fold(0.0, |acc, (x, y)| acc + x * y)
        })
        .collect()
}

impl<'a> Query<'a> {
    pub fn mass_lookup(&self, mass: f32) -> impl Iterator<Item = &PrecursorRange> {
        (self.page_lo..self.page_hi).flat_map(move |page| {
            let left_idx = page * self.bin_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = (left_idx + self.bin_size).min(self.ranges.len());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &self.ranges[left_idx..right_idx];

            let (inner_left, inner_right) = binary_search_slice(
                slice,
                |frag, bounds| frag.mass_lo.total_cmp(bounds),
                mass - 0.1,
                mass + 0.1,
            );

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
                frag.rt <= self.max_rt
                    && frag.rt >= self.min_rt
                    && mass >= frag.mass_lo
                    && mass <= frag.mass_hi
            })
        })
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn disjoint_set() {
        let mut set = DisjointPeakSet::default();

        let a = set.singleton(10, 0.75);
        let b = set.singleton(12, 0.90);
        let c = set.singleton(13, 0.75);
        let d = set.singleton(25, 0.90);

        set.union(a, b, PeakScoringStrategy::SpectralAngle);
        set.union(b, c, PeakScoringStrategy::SpectralAngle);

        assert_eq!(set.peak(a).rt, 12);
        assert_eq!(set.peak(b).rt, 12);
        assert_eq!(set.peak(c).rt, 12);
        assert_eq!(set.peak(d).rt, 25);

        let peaks = set
            .peaks
            .iter()
            .enumerate()
            .filter(|(idx, peak)| peak.parent.get() == *idx)
            .map(|(_, peak)| peak)
            .collect::<Vec<_>>();

        assert_eq!(peaks[0].rt, 12);
        assert_eq!(peaks[1].rt, 25);
    }
}
