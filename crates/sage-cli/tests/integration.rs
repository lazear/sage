use sage_core::database::Builder;
use sage_core::mass::Tolerance;
use sage_core::scoring::Scorer;
use sage_core::spectrum::SpectrumProcessor;

#[test]
fn integration() -> anyhow::Result<()> {
    let mut builder = Builder::default();
    builder.update_fasta("../../tests/Q99536.fasta".into());

    let database = builder.make_parameters().build()?;

    let spectra = sage_core::read_mzml("../../tests/LQSRPAAPPAPGPGQLTLR.mzML", None)?;
    assert_eq!(spectra.len(), 1);

    let sp = SpectrumProcessor::new(100, 0.0, 1500.0, true, 0);
    let processed = sp.process(spectra[0].clone());
    assert!(processed.peaks.len() <= 300);

    let scorer = Scorer {
        db: &database,
        precursor_tol: Tolerance::Ppm(-50.0, 50.0),
        fragment_tol: Tolerance::Ppm(-10.0, 10.0),
        min_matched_peaks: 4,
        min_isotope_err: -1,
        max_isotope_err: 3,
        max_fragment_charge: Some(1),
        min_fragment_mass: 0.0,
        max_fragment_mass: 1500.0,
        chimera: false,
    };

    let psm = scorer.score(&processed, 1);
    assert_eq!(psm.len(), 1);
    assert_eq!(psm[0].peptide, "LQSRPAAPPAPGPGQLTLR");
    assert_eq!(psm[0].matched_peaks, 21);

    Ok(())
}
