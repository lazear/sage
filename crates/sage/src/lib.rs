pub mod database;
pub mod enzyme;
pub mod fasta;
pub mod fdr;
pub mod heap;
pub mod ion_series;
pub mod isotopes;
pub mod lfq;
pub mod mass;
pub mod ml;
pub mod modification;
pub mod peptide;
pub mod scoring;
pub mod spectrum;
pub mod tmt;

fn group_indices_by<X, F: Fn(&X) -> K, K: Ord>(my_vec: &[X], key_fn: F) -> Vec<Vec<usize>> {
    use itertools::Itertools;
    (0..my_vec.len())
        .sorted_by_key(|&i| key_fn(&my_vec[i]))
        .into_iter()
        .group_by(|&i| key_fn(&my_vec[i]))
        .into_iter()
        .map(|(_, group)| group.into_iter().collect())
        .collect()
}
