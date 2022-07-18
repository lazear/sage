use std::collections::BTreeMap;
use std::io;
use std::path::Path;

pub struct Fasta {
    pub proteins: BTreeMap<String, String>,
}

impl Fasta {
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
            if line.starts_with('>') {
                if !s.is_empty() {
                    if last_id.starts_with("Reverse") {
                        s.clear();
                        continue;
                    }
                    let acc = last_id.split('|').nth(1).unwrap().into();
                    map.insert(acc, std::mem::take(&mut s));
                    last_id = &line[1..];
                } else {
                    last_id = &line[1..];
                }
            } else {
                s.push_str(line);
            }
        }

        Ok(Fasta { proteins: map })
    }
}

pub struct Trypsin {
    reverse: bool,
    miss_cleavage: bool,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Digest<'s> {
    pub protein: &'s str,
    pub sequence: String,
    pub reversed: bool,
}

impl<'s> std::hash::Hash for Digest<'s> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sequence.hash(state);
    }
}

impl Trypsin {
    fn inner<'f>(
        &self,
        protein: &'f str,
        sequence: &str,
        reversed: bool,
        peptides: &mut Vec<Digest<'f>>,
    ) {
        let mut left = 0;
        let mut left_ = 0;
        for (right, ch) in sequence.chars().enumerate() {
            match ch {
                'K' | 'R' => {
                    if right + 1 < sequence.len() && sequence[right + 1..].starts_with('P') {
                        continue;
                    }
                    peptides.push(Digest {
                        protein,
                        sequence: sequence[left..=right].to_string(),
                        reversed,
                    });
                    if left_ != left && self.miss_cleavage {
                        peptides.push(Digest {
                            protein,
                            sequence: sequence[left_..=right].to_string(),
                            reversed,
                        });
                    }
                    left_ = left;
                    left = right + 1;
                }
                _ => (),
            }
        }
        peptides.push(Digest {
            protein,
            sequence: sequence[left..].to_string(),
            reversed,
        });
        if left_ != left && self.miss_cleavage {
            peptides.push(Digest {
                protein,
                sequence: sequence[left_..].to_string(),
                reversed,
            });
        }
    }

    pub fn digest<'f>(&self, protein: &'f str, sequence: &str) -> Vec<Digest<'f>> {
        let mut peptides = Vec::new();
        self.inner(protein, sequence, false, &mut peptides);
        if self.reverse {
            let reverse = sequence.chars().rev().collect::<String>();
            self.inner(protein, &reverse, true, &mut peptides);
        }

        peptides
    }

    pub fn new(reverse: bool, miss_cleavage: bool) -> Self {
        Self {
            reverse,
            miss_cleavage,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Trypsin;

    #[test]
    fn digest() {
        let trypsin = Trypsin::new(false, false);
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "R",
            "MSR",
            "SSGR",
            "VYYFNHITNASQWERPSGN",
        ];
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
        let trypsin = Trypsin::new(false, true);
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "MADEEKLPPGWEK",
            "R",
            "LPPGWEKR",
            "MSR",
            "RMSR",
            "SSGR",
            "MSRSSGR",
            "VYYFNHITNASQWERPSGN",
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
}
