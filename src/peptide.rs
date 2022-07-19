use crate::{
    fasta::Digest,
    mass::{Mass, Residue, H2O, VALID_AA},
};

#[derive(Debug, Clone)]
pub enum TargetDecoy {
    Target(Peptide),
    Decoy(Peptide),
}

impl TargetDecoy {
    pub fn neutral(&self) -> f32 {
        match self {
            TargetDecoy::Target(p) | TargetDecoy::Decoy(p) => p.monoisotopic,
        }
    }

    pub fn peptide(&self) -> &Peptide {
        match self {
            TargetDecoy::Target(p) | TargetDecoy::Decoy(p) => p,
        }
    }

    pub fn label(&self) -> i32 {
        match self {
            TargetDecoy::Target(_) => 1,
            TargetDecoy::Decoy(_) => -1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Peptide {
    pub protein: String,
    pub sequence: Vec<Residue>,
    pub monoisotopic: f32,
}

impl Peptide {
    pub fn set_nterm_mod(&mut self, m: f32) {
        if let Some(resi) = self.sequence.first_mut() {
            // Don't overwrite an already modified amino acid!
            match *resi {
                Residue::Just(c) => {
                    self.monoisotopic += m;
                    *resi = Residue::Mod(c, m);
                }
                _ => {}
            };
        }
    }

    pub fn static_mod(&mut self, target: char, m: f32) {
        for resi in self.sequence.iter_mut() {
            // Don't overwrite an already modified amino acid!
            match resi {
                Residue::Just(c) if *c == target => {
                    self.monoisotopic += m;
                    *resi = Residue::Mod(target, m);
                }
                _ => {}
            }
        }
    }
}

impl<'a> TryFrom<&Digest<'a>> for Peptide {
    type Error = char;

    fn try_from(value: &Digest) -> Result<Self, Self::Error> {
        let mut sequence = Vec::with_capacity(value.sequence.len());
        let mut monoisotopic = H2O;

        for c in value.sequence.chars() {
            if !VALID_AA.contains(&c) {
                return Err(c);
            }
            monoisotopic += c.monoisotopic();
            sequence.push(Residue::Just(c));
        }

        Ok(Peptide {
            sequence,
            monoisotopic,
            protein: value.protein.to_string(),
        })
    }
}

impl std::fmt::Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(
            &self
                .sequence
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(""),
        )
    }
}
