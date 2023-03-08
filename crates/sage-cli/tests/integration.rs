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

    let scorer = Scorer::new(
        &database,
        Tolerance::Ppm(-50.0, 50.0),
        Tolerance::Ppm(-10.0, 10.0),
        -1,
        3,
        None,
        0.0,
        1500.0,
        false,
    );

    let psm = scorer.score(&processed, 1);
    assert_eq!(psm.len(), 1);
    assert_eq!(psm[0].peptide, "LQSRPAAPPAPGPGQLTLR");
    assert_eq!(psm[0].matched_peaks, 28);

    Ok(())
}
