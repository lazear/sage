use std::{
    collections::HashMap,
    fmt::{Display, Write},
    str::FromStr,
};

use serde::{
    de::{self, value::MapAccessDeserializer, MapAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};

/// A variable modification with optional per-peptide occurrence limit.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct VariableModification {
    pub mass: f32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_count: Option<usize>,
}

/// A variable modification entry may use the existing bare-mass syntax or the
/// extensible object syntax.
#[derive(Clone, Debug, Serialize, PartialEq)]
#[serde(untagged)]
pub enum VarModEntry {
    Mass(f32),
    Detailed(VariableModification),
}

impl<'de> Deserialize<'de> for VarModEntry {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct VarModEntryVisitor;

        impl VarModEntryVisitor {
            fn mass<E>(value: f64) -> Result<VarModEntry, E>
            where
                E: de::Error,
            {
                if !value.is_finite() || value.abs() > f32::MAX as f64 {
                    return Err(E::custom("variable modification mass is out of range"));
                }
                Ok(VarModEntry::Mass(value as f32))
            }
        }

        impl<'de> Visitor<'de> for VarModEntryVisitor {
            type Value = VarModEntry;

            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str(
                    "a modification mass or an object containing mass and optional max_count",
                )
            }

            fn visit_f64<E>(self, value: f64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value)
            }

            fn visit_f32<E>(self, value: f32) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value as f64)
            }

            fn visit_i64<E>(self, value: i64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value as f64)
            }

            fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Self::mass(value as f64)
            }

            fn visit_map<A>(self, map: A) -> Result<Self::Value, A::Error>
            where
                A: MapAccess<'de>,
            {
                VariableModification::deserialize(MapAccessDeserializer::new(map))
                    .map(VarModEntry::Detailed)
            }
        }

        deserializer.deserialize_any(VarModEntryVisitor)
    }
}

impl VarModEntry {
    pub fn mass(&self) -> f32 {
        match self {
            VarModEntry::Mass(m) => *m,
            VarModEntry::Detailed(modification) => modification.mass,
        }
    }

    pub fn max_count(&self) -> Option<usize> {
        match self {
            VarModEntry::Mass(_) => None,
            VarModEntry::Detailed(modification) => modification.max_count,
        }
    }
}

use crate::mass::VALID_AA;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModificationSpecificity {
    PeptideN(Option<u8>),
    PeptideC(Option<u8>),
    ProteinN(Option<u8>),
    ProteinC(Option<u8>),
    Residue(u8),
}

impl Display for ModificationSpecificity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let r = match self {
            ModificationSpecificity::PeptideN(r) => {
                f.write_char('^')?;
                *r
            }
            ModificationSpecificity::PeptideC(r) => {
                f.write_char('$')?;
                *r
            }
            ModificationSpecificity::ProteinN(r) => {
                f.write_char('[')?;
                *r
            }
            ModificationSpecificity::ProteinC(r) => {
                f.write_char(']')?;
                *r
            }
            ModificationSpecificity::Residue(r) => Some(*r),
        };

        if let Some(r) = r {
            f.write_char(r as char)?;
        }

        Ok(())
    }
}

impl Serialize for ModificationSpecificity {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InvalidModification {
    Empty,
    InvalidResidue(char),
    TooLong(String),
}

impl FromStr for ModificationSpecificity {
    type Err = InvalidModification;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() > 2 {
            return Err(InvalidModification::TooLong(s.into()));
        }
        if let Some(rest) = s.strip_prefix('^') {
            return Ok(ModificationSpecificity::PeptideN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('$') {
            return Ok(ModificationSpecificity::PeptideC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('[') {
            return Ok(ModificationSpecificity::ProteinN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix(']') {
            return Ok(ModificationSpecificity::ProteinC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        match s.chars().next() {
            Some(c) => {
                if VALID_AA.contains(&(c as u8)) {
                    Ok(ModificationSpecificity::Residue(c as u8))
                } else {
                    Err(InvalidModification::InvalidResidue(c))
                }
            }
            None => Err(InvalidModification::Empty),
        }
    }
}

pub fn validate_mods(input: Option<HashMap<String, f32>>) -> HashMap<ModificationSpecificity, f32> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, mass) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(m, mass);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!("Invalid modification string: unrecognized residue ({})", c)
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

pub fn validate_var_mods(
    input: Option<HashMap<String, Vec<VarModEntry>>>,
) -> HashMap<ModificationSpecificity, Vec<VarModEntry>> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, entries) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(m, entries);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Skipping invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!(
                        "Skipping invalid modification string: unrecognized residue ({})",
                        c
                    )
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Skipping invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn parse_modifications() {
        use InvalidModification::*;
        use ModificationSpecificity::*;
        assert_eq!("[".parse::<ModificationSpecificity>(), Ok(ProteinN(None)));
        assert_eq!(
            "[M".parse::<ModificationSpecificity>(),
            Ok(ProteinN(Some(b'M')))
        );
        assert_eq!(
            "]M".parse::<ModificationSpecificity>(),
            Ok(ProteinC(Some(b'M')))
        );
        assert_eq!("M".parse::<ModificationSpecificity>(), Ok(Residue(b'M')));
        assert_eq!(
            "Z".parse::<ModificationSpecificity>(),
            Err(InvalidResidue('Z'))
        );
    }

    #[test]
    fn var_mod_entry_bare_mass() {
        let entry = VarModEntry::Mass(15.9949);
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.max_count(), None);
    }

    #[test]
    fn var_mod_entry_detailed_with_limit() {
        let entry = VarModEntry::Detailed(VariableModification {
            mass: 15.9949,
            max_count: Some(1),
        });
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.max_count(), Some(1));
    }

    #[test]
    fn deserialize_var_mod_entries() {
        let entries: Vec<VarModEntry> = serde_json::from_str(
            r#"[15.9949, {"mass": 42.0106, "max_count": 1}, {"mass": 14.0157}]"#,
        )
        .unwrap();

        assert_eq!(entries.len(), 3);
        assert!((entries[0].mass() - 15.9949).abs() < 1e-4);
        assert_eq!(entries[0].max_count(), None);
        assert!((entries[1].mass() - 42.0106).abs() < 1e-4);
        assert_eq!(entries[1].max_count(), Some(1));
        assert!((entries[2].mass() - 14.0157).abs() < 1e-4);
        assert_eq!(entries[2].max_count(), None);

        let serialized = serde_json::to_value(&entries).unwrap();
        assert!(serialized[0].is_number());
        assert!(serialized[1].is_object());
        assert_eq!(serialized[1]["max_count"], 1);
        assert!(serialized[2].is_object());
        assert!(serialized[2].get("max_count").is_none());

        let round_tripped: Vec<VarModEntry> = serde_json::from_value(serialized).unwrap();
        assert_eq!(entries, round_tripped);
    }

    #[test]
    fn reject_unsupported_var_mod_shapes_and_fields() {
        assert!(serde_json::from_str::<Vec<VarModEntry>>(r#"[[42.0106, 1]]"#).is_err());
        assert!(serde_json::from_str::<Vec<VarModEntry>>(
            r#"[{"mass": 42.0106, "max_count": 1, "name": "Acetyl"}]"#
        )
        .is_err());
        assert!(serde_json::from_str::<Vec<VarModEntry>>(
            r#"[{"mass": 42.0106, "max_counts": 1}]"#
        )
        .is_err());
    }

    #[test]
    fn validate_var_mods_mixed() {
        use ModificationSpecificity::*;
        // Mix bare masses and detailed entries
        let mut raw = HashMap::new();
        raw.insert(
            "M".to_string(),
            vec![
                VarModEntry::Mass(15.9949),
                VarModEntry::Detailed(VariableModification {
                    mass: 15.9949,
                    max_count: Some(1),
                }),
            ],
        );
        raw.insert(
            "C".to_string(),
            vec![VarModEntry::Detailed(VariableModification {
                mass: 57.0215,
                max_count: Some(2),
            })],
        );
        let result = validate_var_mods(Some(raw));

        let m_entries = result.get(&Residue(b'M')).unwrap();
        assert_eq!(m_entries.len(), 2);
        assert!((m_entries[0].mass() - 15.9949).abs() < 1e-4);
        assert_eq!(m_entries[0].max_count(), None);
        assert!((m_entries[1].mass() - 15.9949).abs() < 1e-4);
        assert_eq!(m_entries[1].max_count(), Some(1));

        let c_entries = result.get(&Residue(b'C')).unwrap();
        assert_eq!(c_entries.len(), 1);
        assert!((c_entries[0].mass() - 57.0215).abs() < 1e-4);
        assert_eq!(c_entries[0].max_count(), Some(2));
    }

    #[test]
    fn validate_var_mods_invalid_residue_skipped() {
        let mut raw = HashMap::new();
        raw.insert("Z".to_string(), vec![VarModEntry::Mass(15.9949)]);
        raw.insert("M".to_string(), vec![VarModEntry::Mass(15.9949)]);
        let result = validate_var_mods(Some(raw));
        // Z is invalid — only M should survive
        assert_eq!(result.len(), 1);
        assert!(result.contains_key(&ModificationSpecificity::Residue(b'M')));
    }
}
