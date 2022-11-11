use regex::Regex;

use crate::mass::VALID_AA;

#[derive(Clone, PartialOrd, Ord, Debug, Default)]
/// An enzymatic digest
///
/// # Important invariant about [`Digest`]:
/// * two digests are equal if and only if their sequences are equal
///   i.e., decoy status is ignored for equality and hashing
pub struct Digest {
    /// Is this a decoy peptide?
    pub decoy: bool,
    /// Cleaved peptide sequence
    pub sequence: String,
    /// Missed cleavages
    pub missed_cleavages: u8,
    /// Is this an N-terminal peptide of the protein?
    pub n_terminal: bool,
}

impl Digest {
    /// Generate an internal decoy sequence by reversing the sequence
    /// inside the first and last amino acids
    pub fn reverse(&self) -> Self {
        if self.decoy {
            return self.clone();
        }

        let mut sequence = self.sequence.chars().rev().collect::<Vec<char>>();
        let n = sequence.len().saturating_sub(1);
        sequence.swap(0, n);

        Digest {
            decoy: true,
            sequence: sequence.into_iter().collect(),
            missed_cleavages: self.missed_cleavages,
            n_terminal: self.n_terminal,
        }
    }
}

impl PartialEq for Digest {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence
    }
}

impl Eq for Digest {
    fn assert_receiver_is_total_eq(&self) {}
}

impl std::hash::Hash for Digest {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sequence.hash(state);
    }
}

pub struct EnzymeParameters {
    /// Number of missed cleavages to produce
    pub missed_cleavages: u8,
    /// Inclusive
    pub min_len: usize,
    /// Inclusive
    pub max_len: usize,
    pub enyzme: Enzyme,
}

pub struct Enzyme {
    // Skip cleaving if the site is followed matching this AA
    skip_suffix: Option<char>,
    // Regex for matching cleavage sites
    regex: Regex,
}

impl Enzyme {
    pub fn new(cleave: &str, skip_suffix: Option<char>) -> Self {
        assert!(
            cleave.chars().all(|x| VALID_AA.contains(&x)),
            "Enzyme cleavage sequence contains non-amino acid characters: {}",
            cleave
        );

        assert!(
            skip_suffix.map(|x| VALID_AA.contains(&x)).unwrap_or(true),
            "Enzyme cleavage restriction is non-amino acid character: {}",
            skip_suffix.unwrap(),
        );

        Enzyme {
            skip_suffix,
            regex: Regex::new(&format!("[{}]", cleave)).unwrap(),
        }
    }
}

impl EnzymeParameters {
    fn cleavage_sites(&self, sequence: &str) -> Vec<std::ops::Range<usize>> {
        let mut ranges = Vec::new();
        let mut left = 0;
        for mat in self.enyzme.regex.find_iter(sequence) {
            let right = mat.end();
            if let Some(skip) = self.enyzme.skip_suffix {
                if right < sequence.len() && sequence[right..].starts_with(skip) {
                    continue;
                }
            }
            ranges.push(left..right);
            left = right;
        }
        ranges.push(left..sequence.len());
        ranges
    }

    pub fn digest(&self, sequence: &str) -> Vec<Digest> {
        let mut digests = Vec::new();
        let sites = self.cleavage_sites(sequence);
        for cleavage in 1..=(1 + self.missed_cleavages) {
            // Generate missed cleavages
            for win in sites.windows(cleavage as usize) {
                let sequence = &sequence[win[0].start..win[cleavage as usize - 1].end];
                let len = sequence.len();
                if len >= self.min_len && len <= self.max_len {
                    digests.push(Digest {
                        sequence: sequence.into(),
                        missed_cleavages: cleavage - 1,
                        decoy: false,
                        n_terminal: win[0].start == 0,
                    });
                }
            }
        }
        digests
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn hash_digest() {
        let tryp = EnzymeParameters {
            min_len: 2,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", Some('P')),
        };

        let sequence = "MADEEKMADEEK";
        let expected = vec!["MADEEK", "MADEEK"];

        let mut observed = tryp.digest(sequence);

        // Make sure digest worked!
        assert_eq!(
            expected,
            observed.iter().map(|d| &d.sequence).collect::<Vec<_>>()
        );

        assert_eq!(observed[0].sequence, observed[1].sequence);

        // Make sure hashing a digest works
        let set = observed.drain(..).collect::<HashSet<_>>();
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn trypsin() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec!["MADEEK", "LPPGWEK", "MSR", "SSGR", "VYYFNHITNASQWERPSGN"];

        let tryp = EnzymeParameters {
            min_len: 2,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", Some('P')),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn trypsin_missed_cleavage() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "R",
            "MSR",
            "SSGR",
            "VYYFNHITNASQWERPSGN",
            "MADEEKLPPGWEK",
            "LPPGWEKR",
            "RMSR",
            "MSRSSGR",
            "SSGRVYYFNHITNASQWERPSGN",
        ];

        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: 50,
            missed_cleavages: 1,
            enyzme: Enzyme::new("KR", Some('P')),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn trypsin_missed_cleavage_2() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "R",
            "MSR",
            "SSGR",
            "VYYFNHITNASQWERPSGN",
            "MADEEKLPPGWEK",
            "LPPGWEKR",
            "RMSR",
            "MSRSSGR",
            "SSGRVYYFNHITNASQWERPSGN",
            "MADEEKLPPGWEKR",
            "LPPGWEKRMSR",
            "RMSRSSGR",
            "MSRSSGRVYYFNHITNASQWERPSGN",
        ];

        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: 50,
            missed_cleavages: 2,
            enyzme: Enzyme::new("KR", Some('P')),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_trypsin_pro() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "MSR",
            "SSGR",
            "VYYFNHITNASQWER",
            "PSGN",
        ];

        let tryp = EnzymeParameters {
            min_len: 2,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", None),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_chymotrypsin_pro() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
        let expected = vec![
            "MADEEKL",
            "PPGW",
            "EKRMSRSSGRVY",
            "Y",
            "F",
            "NHITNASQW",
            "ERPSGNW",
        ];

        let tryp = EnzymeParameters {
            min_len: 1,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("FYWL", None),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }
}
