use std::collections::HashMap;
use std::io;
use std::path::Path;

use rayon::prelude::*;

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
        let mut iter = buf.as_str().lines();
        let mut last_id = iter.next().unwrap();
        let mut s = String::new();

        for line in iter {
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
            let acc: String = last_id.split('|').nth(1).unwrap().into();
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

    pub fn peptide_graph(&self, trypsin: &Trypsin) -> HashMap<&str, Vec<&str>> {
        let digests = self
            .targets
            .par_iter()
            .flat_map(|(acc, seq)| trypsin.digest(acc, seq, false))
            .collect::<Vec<_>>();

        let mut map: HashMap<&str, Vec<&str>> = HashMap::new();
        for digest in digests {
            map.entry(digest.sequence).or_default().push(digest.protein);
        }

        map
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
    pub sequence: &'s str,
    /// Missed cleavages
    pub missed_cleavages: u8,
    pub idx: u16,
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
    pub fn digest<'f>(&self, protein: &'f str, sequence: &'f str, is_decoy: bool) -> Vec<Digest<'f>> {
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
                        protein,
                        sequence,
                        missed_cleavages: cleavage - 1,
                        idx,
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

pub struct Kmer {
    pub size: usize,
}

impl Kmer {
    pub fn digest<'a>(&self, protein: &'a str, sequence: &'a str) -> Vec<Digest<'a>> {
        let mut digests = Vec::new();
        for i in 0..sequence.len().saturating_sub(self.size - 1) {
            digests.push(Digest {
                protein,
                sequence: &sequence[i..i + self.size],
                missed_cleavages: 0,
                idx: i as u16,
            })
        }
        digests
    }
}

#[cfg(test)]
mod tests {
    use super::{Kmer, Trypsin};
    use std::collections::HashSet;

    #[test]
    fn hash_digest() {
        let trypsin = Trypsin::new(0, 2, 50);
        let sequence = "MADEEKMADEEK";
        let expected = vec!["MADEEK", "MADEEK"];

        let mut observed = trypsin.digest("A".into(), sequence.into(), false);

        // Make sure digest worked!
        assert_eq!(
            expected,
            observed.iter().map(|d| d.sequence).collect::<Vec<_>>()
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
                .digest("".into(), sequence.into(), false)
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
                .digest("".into(), &sequence, true)
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
                .digest("".into(), sequence.into(), false)
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
                .digest("".into(), sequence.into(), false)
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[test]
    fn kmer() {
        let kmer = Kmer { size: 8 };
        let sequence = "MADEEKLPPGWEKRMSRS";
        let expected = [
            "MADEEKLP", "ADEEKLPP", "DEEKLPPG", "EEKLPPGW", "EKLPPGWE", "KLPPGWEK", "LPPGWEKR",
            "PPGWEKRM", "PGWEKRMS", "GWEKRMSR", "WEKRMSRS",
        ];
        assert_eq!(
            kmer.digest("".into(), sequence.into())
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
        let rev = fwd.chars()
            .rev()
            .collect::<String>();
        
        let mut fwd = trypsin.digest("", fwd, false);
        let mut rev = trypsin.digest("", &rev, true);
        rev.sort_by(|a, b| a.idx.cmp(&b.idx));

        for (f, r) in fwd.iter().zip(rev.iter()) {
            assert_eq!(f.idx, r.idx);
            let r_ = r.sequence[1..r.sequence.len()-1].chars().rev().collect::<String>();
            let f_ = &f.sequence[1..f.sequence.len()-1];
            assert_eq!(f_, r_);
        }

    }
}
