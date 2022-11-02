use fnv::FnvHashMap;
use std::io;
use std::path::Path;

use crate::enzyme::{Digest, EnzymeParameters};

pub struct Fasta {
    pub targets: Vec<(String, String)>,
    decoy_tag: String,
    // Should we ignore decoys in the fasta database
    // and generate them internally?
    generate_decoys: bool,
}

impl Fasta {
    /// Open and parse a fasta file
    pub fn open<P: AsRef<Path>, S: Into<String>>(
        path: P,
        decoy_tag: S,
        generate_decoys: bool,
    ) -> io::Result<Fasta> {
        let decoy_tag = decoy_tag.into();
        let buf = std::fs::read_to_string(path)?;

        let mut targets = Vec::new();
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
                    if !acc.contains(&decoy_tag) || !generate_decoys {
                        targets.push((acc, seq));
                    }
                }
                last_id = id;
            } else {
                s.push_str(line);
            }
        }

        if !s.is_empty() {
            let acc: String = last_id.split_ascii_whitespace().next().unwrap().into();
            if !acc.contains(&decoy_tag) || !generate_decoys {
                targets.push((acc, s));
            }
        }

        Ok(Fasta {
            targets,
            decoy_tag,
            generate_decoys,
        })
    }

    pub fn digest(&self, enzyme: &EnzymeParameters) -> FnvHashMap<Digest, Vec<String>> {
        let mut targets: FnvHashMap<Digest, Vec<String>> = FnvHashMap::default();
        let mut decoys: FnvHashMap<Digest, Vec<String>> = FnvHashMap::default();
        for (acc, seq) in &self.targets {
            for mut digest in enzyme.digest(seq) {
                if self.generate_decoys {
                    decoys
                        .entry(digest.reverse())
                        .or_default()
                        .push(format!("{}{}", self.decoy_tag, acc));
                    targets.entry(digest).or_default().push(acc.clone());
                } else if acc.contains(&self.decoy_tag) {
                    digest.decoy = true;
                    decoys.entry(digest).or_default().push(acc.clone());
                } else {
                    targets.entry(digest).or_default().push(acc.clone());
                }
            }
        }

        // Overwrite any decoys that have the same sequence as a target
        for (k, v) in targets {
            decoys.insert(k, v);
        }
        decoys
    }
}
