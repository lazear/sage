use sage::database::{binary_search_slice, PeptideIx, Theoretical};
use sage::fasta::{Digest, Trypsin};
use sage::ion_series::{IonSeries, Kind};
use sage::mass::Tolerance;
use sage::peptide::Peptide;
use sage::spectrum::SpectrumProcessor;
use std::collections::HashMap;
use std::io::BufReader;

const SEQUENCE: &'static str = "
MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYD
KVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIA
VGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVL
FDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDY
HTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMA
LARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDS
VWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN";

#[test]
/// Confirm that a known match is good!
/// Data from PXD001468 - searched with 10ppm precursor and fragment tolerance
/// Top hit after FDR refinement using [mokapot](https://github.com/wfondrie/mokapot)
pub fn peptide_id() -> Result<(), Box<dyn std::error::Error + Send + Sync + 'static>> {
    let spectra = sage::mzml::MzMlReader::read_ms2("tests/LQSRPAAPPAPGPGQLTLR.mzML")?;
    assert_eq!(spectra.len(), 1);

    let sp = SpectrumProcessor::new(100, 2, 1500.0);
    let processed = sp.process(spectra[0].clone()).unwrap();
    assert!(processed.peaks.len() <= 300);

    let sequence = SEQUENCE.split_whitespace().collect::<String>();
    let mut peptides = Trypsin::new(0, 5, 50)
        .digest("Q99536", &sequence)
        .into_iter()
        .map(|dig| Peptide::try_from(&dig))
        .collect::<Result<Vec<Peptide>, _>>()
        .expect("this better parse!");

    peptides.sort_by(|a, b| a.monoisotopic.total_cmp(&b.monoisotopic));

    // assert_eq!(peptides.len(), 52);

    // let peptides
    let mut hit_index = 0;
    let mut fragments = peptides
        .iter()
        .enumerate()
        .flat_map(|(idx, peptide)| {
            if peptide.to_string() == "LQSRPAAPPAPGPGQLTLR" {
                hit_index = idx;
            }
            IonSeries::new(peptide, Kind::B)
                .chain(IonSeries::new(peptide, Kind::Y))
                .map(move |ion| Theoretical {
                    peptide_index: PeptideIx::for_testing_only_seriously_though(idx),
                    fragment_mz: ion.monoisotopic_mass,
                    kind: ion.kind,
                })
        })
        .collect::<Vec<Theoretical>>();

    let s = serde_json::to_string_pretty(&fragments)?;
    std::fs::write("ions.json", s)?;
    std::fs::write(
        "peptides.json",
        serde_json::to_string_pretty(
            &peptides
                .iter()
                .map(|p| p.to_string())
                .collect::<Vec<String>>(),
        )?,
    )?;

    fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

    // Track scores for all peptides in this protein
    let mut scores = HashMap::new();

    for (mz, _) in processed.peaks {
        let (low, high) = Tolerance::Ppm(-10.0, 10.0).bounds(mz);
        let (i, j) = binary_search_slice(&fragments, |f, x| f.fragment_mz.total_cmp(x), low, high);
        for fragment in fragments[i..j]
            .iter()
            .filter(|frag| frag.fragment_mz >= low && frag.fragment_mz <= high)
        {
            eprintln!(
                "m/z {} matched {} ({})",
                mz,
                fragment.fragment_mz,
                (mz - fragment.fragment_mz).abs()
            );
            let entry = scores.entry(fragment.peptide_index).or_insert(0);
            *entry += 1;
        }
    }

    assert_eq!(
        scores[&PeptideIx::for_testing_only_seriously_though(hit_index)],
        27
    );

    Ok(())
}

#[test]
// We use a funky strategy to simulate charge states (see [`SpectrumProcessor`])
// Confirm that we see the right ID's!
pub fn confirm_charge_state_simulation(
) -> Result<(), Box<dyn std::error::Error + Send + Sync + 'static>> {
    let spectra = sage::mzml::MzMlReader::read_ms2("tests/LQSRPAAPPAPGPGQLTLR.mzML")?;
    assert_eq!(spectra.len(), 1);

    let sp = SpectrumProcessor::new(100, 2, 1500.0);
    let processed = sp.process(spectra[0].clone()).unwrap();
    assert!(processed.peaks.len() <= 300);

    let peptide = Peptide::try_from(&Digest {
        protein: "Q99536",
        sequence: "LQSRPAAPPAPGPGQLTLR".into(),
    })
    .unwrap();

    // Manually generate charge state 2
    let mut fragments = IonSeries::new(&peptide, Kind::B)
        .chain(IonSeries::new(&peptide, Kind::Y))
        .map(move |ion| Theoretical {
            peptide_index: PeptideIx::for_testing_only_seriously_though(0),
            fragment_mz: ion.monoisotopic_mass / 2.0,
            kind: ion.kind,
        })
        .collect::<Vec<Theoretical>>();

    fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));
    assert_eq!(fragments.len(), 36);

    let (matched_b, matched_y) = match_peaks(&fragments, &processed.peaks, 10.0);
    assert_eq!(matched_b + matched_y, 8);

    Ok(())
}

fn match_peaks(fragments: &[Theoretical], peaks: &[(f32, f32)], ppm: f32) -> (usize, usize) {
    let mut matched_b = 0;
    let mut matched_y = 0;
    for (mz, int) in peaks {
        let (low, high) = Tolerance::Ppm(-ppm, ppm).bounds(*mz);
        let (i, j) = binary_search_slice(&fragments, |f, x| f.fragment_mz.total_cmp(x), low, high);
        for fragment in fragments[i..j]
            .iter()
            .filter(|frag| frag.fragment_mz >= low && frag.fragment_mz <= high)
        {
            eprintln!(
                "m/z {} matched {} ({}) int {}",
                mz,
                fragment.fragment_mz,
                (mz - fragment.fragment_mz).abs(),
                int
            );
            match fragment.kind {
                Kind::B => matched_b += 1,
                Kind::Y => matched_y += 1,
            }
        }
    }
    (matched_b, matched_y)
}
