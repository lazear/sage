use rayon::prelude::*;
use sage_core::spectrum::ProcessedSpectrum;
use sage_core::{scoring::Feature, tmt::TmtQuant};

#[derive(Default)]
pub struct SageResults {
    pub ms1: Vec<ProcessedSpectrum>,
    pub features: Vec<Feature>,
    pub quant: Vec<TmtQuant>,
}

impl SageResults {
    fn fold(mut self, other: SageResults) -> Self {
        self.features.extend(other.features);
        self.quant.extend(other.quant);
        self.ms1.extend(other.ms1);
        self
    }
}

impl FromParallelIterator<SageResults> for SageResults {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = SageResults>,
    {
        par_iter
            .into_par_iter()
            .reduce(SageResults::default, SageResults::fold)
    }
}

impl FromIterator<SageResults> for SageResults {
    fn from_iter<I>(par_iter: I) -> Self
    where
        I: IntoIterator<Item = SageResults>,
    {
        par_iter
            .into_iter()
            .fold(SageResults::default(), SageResults::fold)
    }
}
