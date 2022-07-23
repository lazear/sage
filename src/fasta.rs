use std::collections::BTreeMap;
use std::io;
use std::path::Path;

pub struct Fasta {
    pub proteins: BTreeMap<String, String>,
}

impl Fasta {
    /// Open and parse a fasta file
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Fasta> {
        let buf = std::fs::read_to_string(path)?;

        let mut map = BTreeMap::new();
        let mut iter = buf.as_str().lines();
        let mut last_id = iter.next().unwrap();
        let mut s = String::new();

        for line in iter {
            if line.is_empty() {
                continue;
            }
            if let Some(id) = line.strip_prefix('>') {
                if !s.is_empty() {
                    if last_id.starts_with("Reverse") {
                        s.clear();
                        continue;
                    }
                    let acc = last_id.split('|').nth(1).unwrap().into();
                    map.insert(acc, std::mem::take(&mut s));
                }
                last_id = id;
            } else {
                s.push_str(line);
            }
        }

        Ok(Fasta { proteins: map })
    }
}

pub struct Trypsin {
    miss_cleavage: u8,
    min_len: usize,
    max_len: usize,
}

#[derive(Clone, PartialOrd, Ord)]
/// A tryptic digest
///
/// # Important invariant about [`Digest`]:
/// * two digests are equal if and only if their sequences are equal
///   i.e., protein ID is ignored for equality and hashing
pub struct Digest<'s> {
    /// Parent protein ID
    pub protein: &'s str,
    /// Tryptic peptide sequence
    pub sequence: String,
}

impl<'s> PartialEq for Digest<'s> {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence
    }
}

impl<'s> Eq for Digest<'s> {
    fn assert_receiver_is_total_eq(&self) {}
}

impl<'s> std::hash::Hash for Digest<'s> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sequence.hash(state);
    }
}

impl Trypsin {
    fn inner<'s>(&self, sequence: &'s str) -> Vec<&'s str> {
        let mut digests = Vec::new();
        let mut left = 0;
        for (right, ch) in sequence.chars().enumerate() {
            match ch {
                'K' | 'R' => {
                    if right + 1 < sequence.len() && sequence[right + 1..].starts_with('P') {
                        continue;
                    }
                    digests.push(&sequence[left..=right]);
                    left = right + 1;
                }
                _ => (),
            }
        }
        digests.push(&sequence[left..]);
        digests
    }

    /// Generate a series of tryptic digests for a given `sequence`
    pub fn digest<'f>(&self, protein: &'f str, sequence: &str) -> Vec<Digest<'f>> {
        let mut digests = Vec::new();
        let peptides = self.inner(sequence);
        for cleavage in 1..=(1 + self.miss_cleavage) {
            // Generate missed cleavages
            for win in peptides.windows(cleavage as usize) {
                let len: usize = win.iter().map(|w| w.len()).sum();
                if len >= self.min_len && len <= self.max_len {
                    let sequence = win.concat();
                    digests.push(Digest { protein, sequence })
                }
            }
        }
        digests
    }

    /// Create a new [`Trypsin`] struct, which will use the specified parameters
    /// for all in silico digests
    pub fn new(miss_cleavage: u8, min_len: usize, max_len: usize) -> Self {
        Self {
            miss_cleavage,
            max_len,
            min_len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Trypsin;
    use std::collections::HashSet;

    #[test]
    fn hash_digest() {
        let trypsin = Trypsin::new(0, 2, 50);
        let sequence = "MADEEKMADEEK";
        let expected = vec!["MADEEK", "MADEEK"];

        let mut observed = trypsin.digest("A".into(), sequence.into());

        // Make sure digest worked!
        assert_eq!(
            expected,
            observed.iter().map(|d| &d.sequence).collect::<Vec<_>>()
        );

        assert_eq!(observed[0].sequence, observed[1].sequence);

        observed[1].protein = "A";
        // Make sure hashing a digest works
        let set = observed.drain(..).collect::<HashSet<_>>();
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn digest() {
        let trypsin = Trypsin::new(0, 2, 50);
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec!["MADEEK", "LPPGWEK", "MSR", "SSGR", "VYYFNHITNASQWERPSGN"];
        // assert_eq!(super::digest(sequence, false), expected);
        assert_eq!(
            trypsin
                .digest("".into(), sequence.into())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[test]
    fn reverse() {
        let trypsin = Trypsin::new(0, 2, 50);
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN"
            .chars()
            .rev()
            .collect::<String>();
        let expected = vec![
            "NGSPR",
            "EWQSANTIHNFYYVR",
            "GSSR",
            "SMR",
            "EWGPPLK",
            "EEDAM",
        ];

        assert_eq!(
            trypsin
                .digest("".into(), &sequence)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[test]
    fn digest_missed_cleavage() {
        let trypsin = Trypsin::new(1, 0, 50);
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
        assert_eq!(
            trypsin
                .digest("".into(), sequence.into())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[test]
    fn digest_missed_cleavage_2() {
        let trypsin = Trypsin::new(2, 0, 50);
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
        assert_eq!(
            trypsin
                .digest("".into(), sequence.into())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>(),
            expected
        );
    }
}
