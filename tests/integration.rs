use sage::database::{binary_search_slice, PeptideIx, Theoretical};
use sage::fasta::Trypsin;
use sage::ion_series::{IonSeries, Kind};
use sage::mass::Tolerance;
use sage::peptide::Peptide;
use sage::spectrum::{Peak, SpectrumProcessor};
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
pub fn peptide_id() -> Result<(), Box<dyn std::error::Error + Send + Sync + 'static>> {
    let spectra = sage::mzml::MzMlReader::read_ms2("tests/LQSRPAAPPAPGPGQLTLR.mzML")?;
    assert_eq!(spectra.len(), 1);

    let sp = SpectrumProcessor::new(100, 0.0, 1500.0, true);
    let processed = sp.process(spectra[0].clone());
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
                    peptide_index: PeptideIx(idx as u32),
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

    for Peak { mass: mz, .. } in processed.peaks {
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
        scores[&PeptideIx(hit_index as u32)],
        26
    );

    Ok(())
}
