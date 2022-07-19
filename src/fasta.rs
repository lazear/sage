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
    reverse: bool,
    miss_cleavage: u8,
    min_len: usize,
    max_len: usize,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Digest<'s> {
    /// Parent protein ID
    pub protein: &'s str,
    /// Tryptic peptide sequence
    pub sequence: String,
    /// Reversed sequence?
    pub reversed: bool,
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

    fn digest_one_dir<'f>(
        &self,
        protein: &'f str,
        sequence: &str,
        reversed: bool,
        digests: &mut Vec<Digest<'f>>,
    ) {
        let peptides = self.inner(sequence);
        for cleavage in 1..=(1 + self.miss_cleavage) {
            // Generate missed cleavages
            for win in peptides.windows(cleavage as usize) {
                let len: usize = win.iter().map(|w| w.len()).sum();
                if len >= self.min_len && len <= self.max_len {
                    let sequence = win.concat();
                    digests.push(Digest {
                        protein,
                        sequence,
                        reversed,
                    })
                }
            }
        }
    }

    /// Generate a series of tryptic digests for a given `sequence`
    pub fn digest<'f>(&self, protein: &'f str, sequence: &str) -> Vec<Digest<'f>> {
        let mut digests = Vec::new();
        self.digest_one_dir(protein, sequence, false, &mut digests);

        if self.reverse {
            let sequence = sequence.chars().rev().collect::<String>();
            self.digest_one_dir(protein, &sequence, true, &mut digests);
        }
        digests
    }

    /// Create a new [`Trypsin`] struct, which will use the specified parameters
    /// for all in silico digests
    pub fn new(reverse: bool, miss_cleavage: u8, min_len: usize, max_len: usize) -> Self {
        Self {
            reverse,
            miss_cleavage,
            max_len,
            min_len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Trypsin;

    #[test]
    fn digest() {
        let trypsin = Trypsin::new(false, 0, 2, 50);
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
    fn digest_missed_cleavage() {
        let trypsin = Trypsin::new(false, 1, 0, 50);
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
        let trypsin = Trypsin::new(false, 2, 0, 50);
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
