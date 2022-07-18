use crate::{
    fasta::Digest,
    mass::{Mass, Modification, Residue, H2O},
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
    pub fn set_nterm_mod(&mut self, m: Modification) {
        if let Some(resi) = self.sequence.first_mut() {
            *resi = Residue::Mod(Box::new(resi.clone()), m);
            self.monoisotopic += m.monoisotopic();
        }
    }

    pub fn static_mod(&mut self, target: &Residue, m: Modification) {
        for resi in self.sequence.iter_mut() {
            if resi == target {
                *resi = Residue::Mod(Box::new(resi.clone()), m);
                self.monoisotopic += m.monoisotopic();
            }
        }
    }

    pub fn variable_mod(&self, target: &Residue, m: Modification) -> Vec<Peptide> {
        let mut v = Vec::new();
        if self.sequence.contains(target) {
            let mut modded = self.clone();
            modded.static_mod(target, m);
            v.push(modded);
        }
        v
    }
}

impl<'a> TryFrom<&Digest<'a>> for Peptide {
    type Error = char;

    fn try_from(value: &Digest) -> Result<Self, Self::Error> {
        let mut sequence = Vec::with_capacity(value.sequence.len());
        let mut monoisotopic = H2O;

        for c in value.sequence.chars() {
            let r = Residue::try_from(c)?;
            monoisotopic += r.monoisotopic();
            sequence.push(r);
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
