use rayon::prelude::*;
use sage_core::spectrum::ProcessedSpectrum;
use sage_core::{scoring::Feature, tmt::TmtQuant};

#[derive(Default)]
pub struct SageResults {
    pub ms1: Vec<ProcessedSpectrum>,
    pub features: Vec<Feature>,
    pub quant: Vec<TmtQuant>,
}

impl FromParallelIterator<SageResults> for SageResults {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = SageResults>,
    {
        par_iter
            .into_par_iter()
            .reduce(SageResults::default, |mut acc, x| {
                acc.features.extend(x.features);
                acc.quant.extend(x.quant);
                acc.ms1.extend(x.ms1);
                acc
            })
    }
}

impl FromIterator<SageResults> for SageResults {
    fn from_iter<I>(par_iter: I) -> Self
    where
        I: IntoIterator<Item = SageResults>,
    {
        par_iter
            .into_iter()
            .fold(SageResults::default(), |mut acc, x| {
                acc.features.extend(x.features);
                acc.quant.extend(x.quant);
                acc.ms1.extend(x.ms1);
                acc
            })
    }
}
