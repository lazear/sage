use dashmap::DashMap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::hash::BuildHasherDefault;

use crate::database::{binary_search_slice, PeptideIx};
use crate::mass::{Tolerance, NEUTRON};
use crate::scoring::Feature;
use crate::spectrum::ProcessedSpectrum;

pub fn quantify(features: &mut [Feature], spectra: &[ProcessedSpectrum]) {
    let index = LfqIndex::new(features);
    index.quantify(spectra, features)
}

/// Create a data structure analogous to [`IndexedDatabase`] - instaed of
/// storing fragment masses binned by precursor mass, store MS1 precursors
/// binned by RT - This should enable rapid quantification as well
pub struct LfqIndex {
    entries: Vec<LfqEntry>,
    min_rts: Vec<f32>,
    bin_size: usize,
}

#[derive(Copy, Clone, Debug)]
struct LfqEntry {
    rt: f32,
    charge: u8,
    isotope: i32,

    // Pre-computed m/z bounds
    mz_tol_lo: f32,
    mz_tol_hi: f32,

    peptide: PeptideIx,
}

struct Query<'a> {
    entries: &'a [LfqEntry],
    page_lo: usize,
    page_hi: usize,
    bin_size: usize,

    min_rt: f32,
    max_rt: f32,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
struct Quant {
    charge: u8,
    mass: f32,
    intensity: f32,
    isotope: i32,
    rt: f32,
    #[serde(skip_serializing, skip_deserializing)]
    peptide_ix: PeptideIx,
}

impl LfqIndex {
    pub fn new(features: &[Feature]) -> Self {
        let min_charge = 2u8;
        let max_charge = features.iter().map(|feat| feat.charge).max().unwrap_or(4);

        let mut entries = features
            .par_iter()
            .flat_map(|feat| {
                (min_charge..=max_charge)
                    .par_bridge()
                    .flat_map_iter(move |charge| {
                        (0..=6).map(move |isotopomer| {
                            let mono_mass = (feat.expmass - feat.isotope_error
                                + isotopomer as f32 * NEUTRON)
                                / charge as f32;
                            let (mz_tol_lo, mz_tol_hi) =
                                Tolerance::Ppm(-10.0, 10.0).bounds(mono_mass);

                            LfqEntry {
                                rt: feat.rt,
                                peptide: feat.peptide_idx,
                                isotope: isotopomer,
                                charge,
                                mz_tol_lo,
                                mz_tol_hi,
                            }
                        })
                    })
            })
            .collect::<Vec<_>>();

        entries.par_sort_unstable_by(|a, b| a.rt.total_cmp(&b.rt));
        let min_rts = entries
            .par_chunks_mut(1024)
            .map(|chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].rt;
                chunk.par_sort_unstable_by(|a, b| a.mz_tol_lo.total_cmp(&b.mz_tol_lo));
                min
            })
            .collect::<Vec<_>>();

        LfqIndex {
            entries,
            min_rts,
            bin_size: 1024,
        }
    }

    fn rt_slice(&self, rt: f32, rt_tol: f32) -> Query<'_> {
        let (page_lo, page_hi) = binary_search_slice(
            &self.min_rts,
            |rt, x| rt.total_cmp(x),
            rt - rt_tol,
            rt + rt_tol,
        );

        Query {
            entries: &self.entries,
            page_lo,
            page_hi,
            bin_size: self.bin_size,
            max_rt: rt + rt_tol,
            min_rt: rt - rt_tol,
        }
    }

    pub fn quantify(&self, spectra: &[ProcessedSpectrum], features: &mut [Feature]) {
        log::trace!("LFQ");

        let scores: DashMap<PeptideIx, Vec<_>, BuildHasherDefault<fnv::FnvHasher>> =
            DashMap::default();
        spectra
            .par_iter()
            .filter(|s| s.level == 1)
            .for_each(|spectrum| {
                let query = self.rt_slice(spectrum.scan_start_time, 1.0);
                for peak in &spectrum.peaks {
                    let mut seen = HashSet::new();
                    for entry in query.mass_lookup(peak.mass) {
                        if seen.insert(entry.peptide) {
                            let quant = Quant {
                                rt: spectrum.scan_start_time,
                                charge: entry.charge,
                                mass: peak.mass,
                                intensity: peak.intensity,
                                isotope: entry.isotope,
                                peptide_ix: entry.peptide,
                            };
                            scores.entry(quant.peptide_ix).or_default().push(quant);
                        }
                    }
                }
            });

        // Perform gaussian smoothing on MS1 intensities
        let kernel = gaussian_kernel(0.5, 10);

        // Iterate through list of quantification points, considering
        // all MS1 XICs for a given peptide (e.g. aggregating across PSMs)
        let map: HashMap<PeptideIx, Area> = scores
            .into_par_iter()
            .filter(|(_, scores)| scores.len() > 2)
            .map(|(peptide, mut scores)| {
                scores.par_sort_unstable_by(|a, b| a.rt.total_cmp(&b.rt));
                let rt_max = scores.last().unwrap().rt;
                let rt_min = scores.first().unwrap().rt;
                let area = integrate(&scores, &kernel, rt_min - 0.25, rt_max + 0.25);
                (peptide, area)
            })
            .collect();

        features.par_iter_mut().for_each(|feature| {
            if let Some(area) = map.get(&feature.peptide_idx) {
                feature.ms1_intensity = area.integrated_area;
                feature.ms1_apex = area.apex;
                feature.ms1_apex_rt = area.rt;
            }
        })
    }
}

pub struct Area {
    apex: f32,
    integrated_area: f32,
    // fwhm: f32,
    rt: f32,
}

fn integrate(points: &[Quant], kernel: &[f32], rt_min: f32, rt_max: f32) -> Area {
    let grid_size = 100usize;
    let step = (rt_max - rt_min) / (grid_size - 1) as f32;

    // Calculate apex MS1 intensity across all isotopologues/charge states
    let mut apex = vec![0.0f32; grid_size];
    for point in points {
        let ix = ((point.rt - rt_min) / step).floor() as usize;
        apex[ix] += point.intensity;
        // apex[ix] = apex[ix].max(point.intensity);
    }

    let smoothed = convolve(&apex, kernel);

    // Calculate KDE-smoothed density of MS1 ions
    let rt = points.iter().map(|pt| pt.rt as f64).collect::<Vec<_>>();
    let kde = crate::ml::kde::Kde::new(&rt);

    // Scaled the smoothed MS1 intensites by MS1 ion density
    let density = (0..grid_size)
        .map(|ix| {
            let rt = (ix as f32 * step + rt_min) as f64;
            kde.pdf(rt) as f32 * smoothed[ix]
        })
        .collect::<Vec<_>>();

    let (left, peak, right) = find_peak(&density);

    let rt = peak as f32 * step + rt_min;
    // let fwhm = (right - left) as f32 * step;
    let (apex, integrated_area) = apex[left..right]
        .iter()
        .fold((0.0f32, 0.0), |(apex, area), &x| (apex.max(x), area + x));

    Area {
        apex,
        integrated_area,
        rt,
        // fwhm,
    }
}

/// Create a symmetrical gaussian kernel of given standard deviation and length
fn gaussian_kernel(sigma: f32, len: usize) -> Vec<f32> {
    let step = 2.0 / (len - 1) as f32;
    let constant = 1.0 / (sigma * (2.0 * std::f32::consts::PI).sqrt());

    let mut kernel = (0..len)
        .map(|i| {
            let x = i as f32 * step - 1.0;
            constant * (-0.5 * (x / sigma).powi(2)).exp()
        })
        .collect::<Vec<_>>();

    let sum = kernel.iter().sum::<f32>();
    kernel.iter_mut().for_each(|x| *x /= sum);
    kernel
}
/// Convolve a signal with a symmetrical kernel
/// - This should behave the same as `np.convolve(..., mode='same')`
fn convolve(slice: &[f32], kernel: &[f32]) -> Vec<f32> {
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

/// Find the first peak with global maximum intensity
/// Return the left, center, and right indices of the peak, where left and right
/// correspond to full-width at half-maximum
fn find_peak(slice: &[f32]) -> (usize, usize, usize) {
    let (idx, height) =
        slice
            .iter()
            .enumerate()
            .fold((0, 0.0), |(peak_idx, peak_height), (idx, x)| {
                if *x > peak_height {
                    (idx, *x)
                } else {
                    (peak_idx, peak_height)
                }
            });

    let (mut left, mut right) = (idx, idx);

    while left > 0 {
        if slice[left] <= height / 2.0 {
            break;
        }
        left -= 1
    }
    while right < slice.len() {
        if slice[right] <= height / 2.0 {
            break;
        }
        right += 1
    }

    (left, idx, right)
}

impl<'a> Query<'a> {
    pub fn mass_lookup(&self, mass: f32) -> impl Iterator<Item = &LfqEntry> {
        (self.page_lo..self.page_hi).flat_map(move |page| {
            let left_idx = page * self.bin_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = (left_idx + self.bin_size).min(self.entries.len());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &self.entries[left_idx..right_idx];

            let (inner_left, inner_right) = binary_search_slice(
                slice,
                |frag, bounds| frag.mz_tol_lo.total_cmp(bounds),
                mass - 1.0,
                mass + 1.0,
            );

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
                frag.rt <= self.max_rt
                    && frag.rt >= self.min_rt
                    && mass >= frag.mz_tol_lo
                    && mass <= frag.mz_tol_hi
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convolve_simple() {
        let array = [1., 2., 3., 5., 4., 6., 4., 7., 8.];
        let kernel = [1., 2., 3., 2., 1.];
        let result = convolve(&array, &kernel);

        let expected = [10., 19., 28., 37., 41., 46., 50., 51., 42.];
        assert_eq!(result, expected);
    }

    #[test]
    fn mk_gaussian_kernel() {
        let kern = gaussian_kernel(0.5, 10);
        let expected = [
            0.024612451,
            0.054237682,
            0.09809815,
            0.14562434,
            0.17742734,
            0.17742734,
            0.14562432,
            0.09809815,
            0.054237682,
            0.024612451,
        ];
        assert_eq!(kern, expected)
    }

    #[test]
    fn pick_peak() {
        let data = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1];

        let (left, apex, right) = find_peak(&data);

        assert_eq!(left, 2);
        assert_eq!(apex, 5);
        assert_eq!(right, 9);
    }

    #[test]
    fn integrate_peak() {
        let json = include_str!("../../../tests/lfq_test.json");
        let pts: Vec<Quant> = serde_json::from_str(&json).unwrap();
        let kernel = gaussian_kernel(0.5, 10);

        let area = integrate(&pts, &kernel, 68.0, 70.0);
        assert!((area.rt - 69.0).abs() < 3.1 / 99.0, "{}", area.rt);
        assert!((area.apex - 477762660.0).abs() < 1.0, "{}", area.apex);
    }
}
