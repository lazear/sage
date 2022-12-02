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
    pub position: Position,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Position {
    Nterm,
    Cterm,
    Full,
    Internal,
}

impl Default for Position {
    fn default() -> Self {
        Self::Internal
    }
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
            position: self.position,
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
    pub enyzme: Option<Enzyme>,
}

pub struct Enzyme {
    // Skip cleaving if the site is followed matching this AA
    skip_suffix: Option<char>,
    // Regex for matching cleavage sites
    regex: Regex,
    // Cleave at c-terminal?
    c_terminal: bool,
}

impl Enzyme {
    pub fn new(cleave: &str, skip_suffix: Option<char>, c_terminal: bool) -> Option<Self> {
        assert!(
            cleave.chars().all(|x| VALID_AA.contains(&x)) || cleave == "$",
            "Enzyme cleavage sequence contains non-amino acid characters: {}",
            cleave
        );

        assert!(
            skip_suffix.map(|x| VALID_AA.contains(&x)).unwrap_or(true),
            "Enzyme cleavage restriction is non-amino acid character: {}",
            skip_suffix.unwrap(),
        );

        // At this point, cleave can be three things: empty, "$", or a string of valid AA's
        match cleave {
            "" => None,
            "$" => Some(Enzyme {
                regex: Regex::new("$").unwrap(),
                skip_suffix: None,
                // Allowing this to be set to false could cause unexpected behavior
                c_terminal: true,
            }),
            _ => Some(Enzyme {
                regex: Regex::new(&format!("[{}]", cleave)).unwrap(),
                skip_suffix,
                c_terminal,
            }),
        }
    }

    fn cleavage_sites(&self, sequence: &str) -> Vec<std::ops::Range<usize>> {
        let mut ranges = Vec::new();
        let mut left = 0;
        for mat in self.regex.find_iter(sequence) {
            let right = match self.c_terminal {
                true => mat.end(),
                false => mat.start(),
            };
            if let Some(skip) = self.skip_suffix {
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
}

impl EnzymeParameters {
    fn cleavage_sites(&self, sequence: &str) -> Vec<std::ops::Range<usize>> {
        match &self.enyzme {
            Some(enzyme) => enzyme.cleavage_sites(sequence),
            None => {
                // Perform a non-specific digest
                let mut v = Vec::new();
                for len in self.min_len..=self.max_len {
                    for i in 0..=sequence.len().saturating_sub(len) {
                        v.push(i..i + len)
                    }
                }
                v
            }
        }
    }

    pub fn digest(&self, sequence: &str) -> Vec<Digest> {
        let n = sequence.len();
        let mut digests = Vec::new();
        let sites = self.cleavage_sites(sequence);
        // Allowing missed_cleavages with non-specific digest causes OOB panics
        // in the below indexing code
        let missed_cleavages = match self.enyzme {
            None => 0,
            _ => self.missed_cleavages,
        };

        for cleavage in 1..=(1 + missed_cleavages) {
            // Generate missed cleavages
            for win in sites.windows(cleavage as usize) {
                // dbg!(&win);
                let start = win[0].start;
                let end = win[cleavage as usize - 1].end;
                // if start >= sequence.len() || end >= sequence.len() {
                //     continue;
                // }

                // let sequence = &sequence[win[0].start..win[cleavage as usize - 1].end];
                let sequence = match sequence.get(start..end) {
                    Some(sequence) => sequence,
                    None => continue,
                };

                let len = sequence.len();

                let position = match (start == 0, end == n) {
                    (true, true) => Position::Full,
                    (true, false) => Position::Nterm,
                    (false, true) => Position::Cterm,
                    (false, false) => Position::Internal,
                };

                if len >= self.min_len && len <= self.max_len && len > 0 {
                    digests.push(Digest {
                        sequence: sequence.into(),
                        missed_cleavages: cleavage - 1,
                        decoy: false,
                        position,
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
            enyzme: Enzyme::new("KR", Some('P'), true),
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
            enyzme: Enzyme::new("KR", Some('P'), true),
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
            enyzme: Enzyme::new("KR", Some('P'), true),
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
            enyzme: Enzyme::new("KR", Some('P'), true),
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
            enyzme: Enzyme::new("KR", None, true),
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
    fn test_asp_n() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
        let expected = vec!["MA", "DEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW"];

        let tryp = EnzymeParameters {
            min_len: 1,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("D", None, false),
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
            enyzme: Enzyme::new("FYWL", None, true),
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
    fn nonspecific_digest_5() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";

        let expected = sequence
            .as_bytes()
            .windows(5)
            .flat_map(std::str::from_utf8)
            .collect::<Vec<_>>();

        let tryp = EnzymeParameters {
            min_len: 5,
            max_len: 5,
            missed_cleavages: 0,
            enyzme: None,
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
    fn nonspecific_digest_5_7() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";

        let expected = (5..=7)
            .flat_map(|window| {
                sequence
                    .as_bytes()
                    .windows(window)
                    .flat_map(std::str::from_utf8)
            })
            .collect::<Vec<_>>();

        let tryp = EnzymeParameters {
            min_len: 5,
            max_len: 7,
            missed_cleavages: 0,
            enyzme: Enzyme::new("", None, true),
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
    fn no_digest() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
        let expected = vec![sequence];

        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: usize::MAX,
            missed_cleavages: 0,
            enyzme: Enzyme::new("$", None, true),
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
