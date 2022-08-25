use std::collections::HashMap;

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

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Peptide {
    pub protein: String,
    pub sequence: Vec<Residue>,
    pub nterm: Option<f32>,
    pub monoisotopic: f32,
}

impl Peptide {
    fn set_nterm_mod(&mut self, m: f32) {
        if self.nterm.is_none() {
            self.nterm = Some(m);
            self.monoisotopic += m;
        }
    }

    /// Apply a static modification to a peptide in-place
    pub fn static_mod(&mut self, residue: char, mass: f32) {
        if residue == '^' {
            return self.set_nterm_mod(mass);
        }

        for resi in self.sequence.iter_mut() {
            // Don't overwrite an already modified amino acid!
            match resi {
                Residue::Just(c) if *c == residue => {
                    self.monoisotopic += mass;
                    *resi = Residue::Mod(residue, mass);
                }
                _ => {}
            }
        }
    }

    /// Create an iterator that will produce all singly-modified peptides with
    /// the given variable modification
    pub fn variable_mod(&self, residue: char, mass: f32) -> VariableMod<'_> {
        VariableMod {
            peptide: self,
            index: 0,
            residue,
            mass,
        }
    }

    /// Apply variable modifications, then static modifications to a peptide
    pub fn apply(
        mut self,
        variable_mods: &HashMap<char, f32>,
        static_mods: &HashMap<char, f32>,
    ) -> Vec<Peptide> {
        if variable_mods.is_empty() {
            for (resi, mass) in static_mods {
                self.static_mod(*resi, *mass);
            }
            vec![self]
        } else {
            let mut peptides =
                variable_mods
                    .iter()
                    .fold(vec![self], |acc, (resi, mass)| {
                        acc.iter()
                            .flat_map(|peptide| peptide.variable_mod(*resi, *mass))
                            .collect()
                    });
            peptides.iter_mut().for_each(|p| {
                for (resi, mass) in static_mods {
                    p.static_mod(*resi, *mass);
                }
            });
            peptides
        }
    }
}

pub struct VariableMod<'a> {
    peptide: &'a Peptide,
    index: usize,
    residue: char,
    mass: f32,
}

impl<'a> Iterator for VariableMod<'a> {
    type Item = Peptide;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index == 0 && self.residue == '^' {
            self.index = self.peptide.sequence.len();
            if self.peptide.nterm.is_none() {
                let mut modified = self.peptide.clone();
                modified.set_nterm_mod(self.mass);
                return Some(modified);
            }
        }

        while self.index < self.peptide.sequence.len() {
            match self.peptide.sequence[self.index] {
                Residue::Just(r) if r == self.residue => {
                    let mut modified = self.peptide.clone();
                    modified.sequence[self.index] = Residue::Mod(r, self.mass);
                    modified.monoisotopic += self.mass;
                    self.index += 1;
                    return Some(modified);
                }
                _ => self.index += 1,
            }
        }
        if self.index == self.peptide.sequence.len() {
            self.index += 1;
            return Some(self.peptide.clone());
        }
        None
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
            nterm: None,
            protein: value.protein.to_string(),
        })
    }
}

impl std::fmt::Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = self.nterm {
            write!(f, "[{}]", m)?;
        }
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_variable_mods() {
        let variable_mods = [('M', 16.), ('C', 57.)];
        let peptide = Peptide::try_from(&Digest {
            protein: "",
            sequence: "GCMGCMG",
        })
        .unwrap();

        let expected = vec![
            "GC(57)M(16)GCMG",
            "GCM(16)GC(57)MG",
            "GCM(16)GCMG",
            "GC(57)MGCM(16)G",
            "GCMGC(57)M(16)G",
            "GCMGCM(16)G",
            "GC(57)MGCMG",
            "GCMGC(57)MG",
            "GCMGCMG",
        ];

        // Fold variable mods over initial peptide, producing permutations in order
        // e.g., the first peptides yielded should contain Met mods at the first Met
        let peptides: Vec<String> = variable_mods
            .iter()
            .fold(vec![peptide], |acc, (resi, mass)| {
                acc.iter()
                    .flat_map(|peptide| peptide.variable_mod(*resi, *mass))
                    .collect()
            })
            .into_iter()
            .map(|peptide| peptide.to_string())
            .collect();

        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_nterm() {
        let variable_mods = [('^', 42.), ('M', 16.)];
        let peptide = Peptide::try_from(&Digest {
            protein: "",
            sequence: "GCMGCMG",
        })
        .unwrap();

        let expected = vec![
            "[42]GCM(16)GCMG",
            "[42]GCMGCM(16)G",
            "[42]GCMGCMG",
            "GCM(16)GCMG",
            "GCMGCM(16)G",
            "GCMGCMG",
        ];

        // Fold variable mods over initial peptide, producing permutations in order
        // e.g., the first peptides yielded should contain Met mods at the first Met
        let peptides: Vec<String> = variable_mods
            .iter()
            .fold(vec![peptide], |acc, (resi, mass)| {
                acc.iter()
                    .flat_map(|peptide| peptide.variable_mod(*resi, *mass))
                    .collect()
            })
            .into_iter()
            .map(|peptide| peptide.to_string())
            .collect();

        assert_eq!(peptides, expected);
    }

    #[test]
    fn apply_mods() {
        let peptide = Peptide::try_from(&Digest {
            protein: "",
            sequence: "AACAACAA",
        })
        .unwrap();

        let expected = vec![
            "AAC(30)AAC(57)AA",
            "AAC(57)AAC(30)AA",
            "AAC(57)AAC(57)AA",
        ]; 

        let mut static_mods = HashMap::new();
        static_mods.insert('C', 57.0);

        let mut variable_mods = HashMap::new();
        variable_mods.insert('C', 30.0);

        let actual = peptide.apply(&variable_mods, &static_mods)
            .into_iter().map(|p| p.to_string()).collect::<Vec<_>>();
        assert_eq!(actual, expected);


    }
}
