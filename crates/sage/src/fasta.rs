use dashmap::DashMap;
use rayon::prelude::*;
use std::hash::BuildHasherDefault;

use crate::enzyme::{Digest, EnzymeParameters};

pub struct Fasta {
    pub targets: Vec<(String, String)>,
    decoy_tag: String,
    // Should we ignore decoys in the fasta database
    // and generate them internally?
    generate_decoys: bool,
}

impl Fasta {
    // Parse a string into a fasta database
    pub fn parse<S: Into<String>>(contents: String, decoy_tag: S, generate_decoys: bool) -> Fasta {
        let decoy_tag = decoy_tag.into();

        let mut targets = Vec::new();
        let mut last_id = "";
        let mut s = String::new();

        for line in contents.as_str().lines() {
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

        Fasta {
            targets,
            decoy_tag,
            generate_decoys,
        }
    }

    pub fn digest(
        &self,
        enzyme: &EnzymeParameters,
    ) -> DashMap<Digest, Vec<String>, BuildHasherDefault<fnv::FnvHasher>> {
        let targets: DashMap<Digest, Vec<String>, BuildHasherDefault<fnv::FnvHasher>> =
            DashMap::default();
        let decoys: DashMap<Digest, Vec<String>, BuildHasherDefault<fnv::FnvHasher>> =
            DashMap::default();

        self.targets.par_iter().for_each(|(acc, seq)| {
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
        });

        // Overwrite any decoys that have the same sequence as a target
        targets.into_par_iter().for_each(|(k, v)| {
            // If we don't remove existing decoy entries, we will end up with
            // target PSMs having `decoy=true` and `label=1`, because HashMap::insert
            // does not modify the existing entry (and we define digests to be equal
            // ignoring decoy status)
            if decoys.contains_key(&k) {
                decoys.remove(&k);
            }
            decoys.insert(k, v);
        });

        decoys
    }
}
