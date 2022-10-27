use crate::scoring::Percolator;

/// Assign q_values in place to a set of PSMs, returning the number of PSMs
/// q <= 0.01
///
/// # Invariants
/// * `scores` must be sorted in descending order (e.g. best PSM is first)
pub fn spectrum_q_value(scores: &mut [Percolator]) -> usize {
    // FDR Calculation:
    // * Sort by score, descending
    // * Estimate FDR
    // * Calculate q-value

    let mut decoy = 1;
    let mut target = 0;

    for score in scores.iter_mut() {
        match score.label == -1 {
            true => decoy += 1,
            false => target += 1,
        }
        score.q_value = decoy as f32 / target as f32;
    }

    // Reverse slice, and calculate the cumulative minimum
    let mut q_min = 1.0f32;
    let mut passing = 0;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.q_value);
        score.q_value = q_min;
        if q_min <= 0.01 {
            passing += 1;
        }
    }
    passing
}
