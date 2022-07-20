use carina::database::{binary_search_slice, PeptideIx, Theoretical};
use carina::fasta::{Digest, Trypsin};
use carina::ion_series::{IonSeries, Kind};
use carina::mass::Tolerance;
use carina::peptide::Peptide;
use carina::spectrum::{Spectrum, SpectrumProcessor};
use std::collections::HashMap;

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
pub fn peptide_id() -> Result<(), Box<dyn std::error::Error>> {
    let data: &str = include_str!("LQSRPAAPPAPGPGQLTLR.json");
    let spectrum: Spectrum = serde_json::from_str(data)?;

    let sp = SpectrumProcessor::new(100, 2, 1500.0);
    let processed = sp.process(spectrum);
    assert!(processed.peaks.len() <= 300);

    let sequence = SEQUENCE.split_whitespace().collect::<String>();
    let peptides = Trypsin::new(false, 1, 5, 50)
        .digest("Q99536", &sequence)
        .into_iter()
        .map(|dig| Peptide::try_from(&dig))
        .collect::<Result<Vec<Peptide>, _>>()
        .expect("this better parse!");

    assert_eq!(peptides.len(), 52);

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
                    precursor_mz: peptide.monoisotopic,
                    fragment_mz: ion.monoisotopic_mass,
                    kind: ion.kind,
                })
        })
        .collect::<Vec<Theoretical>>();

    fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

    // Track scores for all peptides in this protein
    let mut scores = HashMap::new();

    for (mz, _) in processed.peaks {
        let (low, high) = Tolerance::Ppm(10.0).bounds(mz);
        let (i, j) = binary_search_slice(&fragments, |f| f.fragment_mz, low, high);
        for fragment in fragments[i..j]
            .iter()
            .filter(|frag| frag.fragment_mz >= low && frag.fragment_mz <= high)
        {
            // eprintln!("m/z {} matched {} ({}) int {}", mz, fragment.fragment_mz, (mz - fragment.fragment_mz).abs(), int);
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
pub fn confirm_charge_state_simulation() -> Result<(), Box<dyn std::error::Error>> {
    let data: &str = include_str!("LQSRPAAPPAPGPGQLTLR.json");
    let spectrum: Spectrum = serde_json::from_str(data)?;

    // charge state = 1, e.g. don't simulate higher states
    let sp = SpectrumProcessor::new(100, 1, 1500.0);
    let processed = sp.process(spectrum);
    assert!(processed.peaks.len() <= 300);

    let peptide = Peptide::try_from(&Digest {
        protein: "Q99536",
        sequence: "LQSRPAAPPAPGPGQLTLR".into(),
        reversed: false,
    })
    .unwrap();

    // Manually generate charge state 2
    let mut fragments = IonSeries::new(&peptide, Kind::B)
        .chain(IonSeries::new(&peptide, Kind::Y))
        .map(move |ion| Theoretical {
            peptide_index: PeptideIx::for_testing_only_seriously_though(0),
            precursor_mz: peptide.monoisotopic,
            fragment_mz: ion.monoisotopic_mass / 2.0,
            kind: ion.kind,
        })
        .collect::<Vec<Theoretical>>();

    fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));
    assert_eq!(fragments.len(), 36);

    let mut matched_b = 0;
    let mut matched_y = 0;

    for (mz, int) in processed.peaks {
        let (low, high) = Tolerance::Ppm(10.0).bounds(mz);
        let (i, j) = binary_search_slice(&fragments, |f| f.fragment_mz, low, high);
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

    assert_eq!(matched_b + matched_y, 8);

    Ok(())
}
