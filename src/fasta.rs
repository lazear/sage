use std::collections::HashMap;
use std::io;
use std::path::Path;

pub struct Fasta {
    pub targets: Vec<(String, String)>,
    pub decoys: Vec<(String, String)>,
}

impl Fasta {
    /// Open and parse a fasta file
    pub fn open<P: AsRef<Path>>(path: P, decoy_prefix: &str) -> io::Result<Fasta> {
        let buf = std::fs::read_to_string(path)?;

        let mut targets = Vec::new();
        let mut decoys = Vec::new();
        let mut last_id = "";
        let mut s = String::new();

        for line in buf.as_str().lines() {
            if line.is_empty() {
                continue;
            }
            if let Some(id) = line.strip_prefix('>') {
                if !s.is_empty() {
                    let acc: String = last_id.split_ascii_whitespace().next().unwrap().into();
                    let seq = std::mem::take(&mut s);
                    match acc.starts_with(decoy_prefix) {
                        true => decoys.push((acc, seq)),
                        false => targets.push((acc, seq)),
                    }
                }
                last_id = id;
            } else {
                s.push_str(line);
            }
        }

        if !s.is_empty() {
            let acc: String = last_id.split_ascii_whitespace().next().unwrap().into();
            match acc.starts_with(decoy_prefix) {
                true => decoys.push((acc, s)),
                false => targets.push((acc, s)),
            }
        }

        Ok(Fasta { targets, decoys })
    }

    pub fn make_decoys(&mut self, decoy_prefix: &str) {
        if !self.decoys.is_empty() {
            return;
        }
        let decoys = self
            .targets
            .iter()
            .map(|(p, s)| (format!("{}{}", decoy_prefix, p), s.chars().rev().collect()))
            .collect::<Vec<_>>();

        self.decoys = decoys;
    }

    pub fn digest(&self, trypsin: &Trypsin) -> HashMap<Digest, Vec<String>> {
        let mut targets: HashMap<Digest, Vec<String>> = HashMap::new();
        for (acc, seq) in &self.targets {
            for digest in trypsin.digest(seq, false) {
                targets.entry(digest).or_default().push(acc.clone());
            }
        }
        let mut decoys: HashMap<Digest, Vec<String>> = HashMap::new();
        for (acc, seq) in &self.decoys {
            for digest in trypsin.digest(seq, true) {
                decoys.entry(digest).or_default().push(acc.clone());
            }
        }
        // Overwrite decoy entries
        for (k, v) in targets {
            decoys.insert(k, v);
        }
        decoys
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
pub struct Digest {
    pub decoy: bool,
    /// Tryptic peptide sequence
    pub sequence: String,
    /// Missed cleavages
    pub missed_cleavages: u8,
    pub idx: u16,
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

impl Trypsin {
    fn inner(&self, sequence: &str) -> Vec<std::ops::Range<usize>> {
        let mut digests = Vec::new();
        let mut left = 0;
        for (right, ch) in sequence.chars().enumerate() {
            match ch {
                'K' | 'R' => {
                    if right + 1 < sequence.len() && sequence[right + 1..].starts_with('P') {
                        continue;
                    }
                    digests.push(left..right + 1);
                    left = right + 1;
                }
                _ => (),
            }
        }
        digests.push(left..sequence.len());
        digests
    }

    /// Generate a series of tryptic digests for a given `sequence`
    pub fn digest(&self, sequence: &str, is_decoy: bool) -> Vec<Digest> {
        let mut digests = Vec::new();
        let mut idx = 0;
        let peptides = self.inner(sequence);
        for cleavage in 1..=(1 + self.miss_cleavage) {
            // Generate missed cleavages
            for win in peptides.windows(cleavage as usize) {
                let sequence = &sequence[win[0].start..win[cleavage as usize - 1].end];
                let len = sequence.len();
                if len >= self.min_len && len <= self.max_len {
                    digests.push(Digest {
                        // protein,
                        sequence: sequence.into(),
                        missed_cleavages: cleavage - 1,
                        idx,
                        decoy: is_decoy,
                    });
                    idx += 1;
                }
            }
        }
        if is_decoy {
            let n = digests.len() as u16;
            digests.iter_mut().for_each(|dig| dig.idx = n - dig.idx - 1);
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

        let mut observed = trypsin.digest(sequence.into(), false);

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
    fn digest() {
        let trypsin = Trypsin::new(0, 2, 50);
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec!["MADEEK", "LPPGWEK", "MSR", "SSGR", "VYYFNHITNASQWERPSGN"];
        assert_eq!(
            trypsin
                .digest(sequence.into(), false)
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
                .digest(&sequence, true)
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
                .digest(sequence.into(), false)
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
                .digest(sequence.into(), false)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>(),
            expected
        );
    }

    /// Check that picked-peptide approach will match forward and reverse peptides
    #[test]
    fn digest_index() {
        let trypsin = Trypsin::new(0, 2, 50);

        let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let rev = fwd.chars().rev().collect::<String>();

        let fwd = trypsin.digest(fwd, false);
        let mut rev = trypsin.digest(&rev, true);
        rev.sort_by(|a, b| a.idx.cmp(&b.idx));

        for (f, r) in fwd.iter().zip(rev.iter()) {
            assert_eq!(f.idx, r.idx);
            let r_ = r.sequence[..r.sequence.len() - 1]
                .chars()
                .rev()
                .skip(1)
                .collect::<String>();
            let f_ = &f.sequence[1..f.sequence.len() - 1];
            assert!(f_.contains(&r_), "{} {}", f_, r_);
        }
    }
}
