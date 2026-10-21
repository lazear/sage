use std::{
    fmt::{Display, Write},
    num::NonZeroU8,
    str::FromStr,
    sync::Arc,
};

use serde::{
    de::{self, value::MapAccessDeserializer, MapAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};

/// A modification with an optional ProForma label and, for variable
/// modifications, an optional per-peptide occurrence limit.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct DetailedMod {
    pub mass: f32,
    /// Ignored for static modifications
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_count: Option<usize>,
    /// ProForma label reported instead of the mass, e.g. `UNIMOD:35`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<ModName>,
}

/// A modification entry may use the existing bare-mass syntax or the
/// extensible object syntax.
#[derive(Clone, Debug, Serialize, PartialEq)]
#[serde(untagged)]
pub enum ModEntry {
    Mass(f32),
    Detailed(DetailedMod),
}

impl<'de> Deserialize<'de> for ModEntry {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ModEntryVisitor;

        impl ModEntryVisitor {
            fn mass<E>(value: f64) -> Result<ModEntry, E>
            where
                E: de::Error,
            {
                if !value.is_finite() || value.abs() > f32::MAX as f64 {
                    return Err(E::custom("modification mass is out of range"));
                }
                Ok(ModEntry::Mass(value as f32))
            }
        }

        impl<'de> Visitor<'de> for ModEntryVisitor {
            type Value = ModEntry;

            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str(
                    "a modification mass or an object containing mass and optional max_count/name",
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
                DetailedMod::deserialize(MapAccessDeserializer::new(map)).map(ModEntry::Detailed)
            }
        }

        deserializer.deserialize_any(ModEntryVisitor)
    }
}

impl ModEntry {
    pub fn mass(&self) -> f32 {
        match self {
            ModEntry::Mass(m) => *m,
            ModEntry::Detailed(modification) => modification.mass,
        }
    }

    pub fn max_count(&self) -> Option<usize> {
        match self {
            ModEntry::Mass(_) => None,
            ModEntry::Detailed(modification) => modification.max_count,
        }
    }

    pub fn name(&self) -> Option<&ModName> {
        match self {
            ModEntry::Mass(_) => None,
            ModEntry::Detailed(modification) => modification.name.as_ref(),
        }
    }
}

/// A ProForma modification label, e.g. `UNIMOD:35`: non-empty, and without
/// square brackets so that it can be written as `M[UNIMOD:35]`
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ModName(Arc<str>);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InvalidModName(pub String);

impl Display for InvalidModName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "invalid modification name {:?}: must be non-empty and cannot contain '[' or ']'",
            self.0
        )
    }
}

impl std::error::Error for InvalidModName {}

impl FromStr for ModName {
    type Err = InvalidModName;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() || s.contains(['[', ']']) {
            return Err(InvalidModName(s.into()));
        }
        Ok(ModName(Arc::from(s)))
    }
}

impl ModName {
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for ModName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl Serialize for ModName {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.0)
    }
}

impl<'de> Deserialize<'de> for ModName {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(de::Error::custom)
    }
}

/// Index into a [`Modifications`] registry. Only meaningful together with the
/// registry that issued it - in practice, the [`crate::database::IndexedDatabase`]
/// that a peptide belongs to.
///
/// `Option<ModId>` is a single byte, which caps a registry at 255 distinct
/// modifications. Switch to `NonZeroU16` if that is ever too few.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ModId(NonZeroU8);

impl ModId {
    fn index(self) -> usize {
        self.0.get() as usize - 1
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Modification {
    pub mass: f32,
    pub name: Option<ModName>,
}

impl Display for Modification {
    /// ProForma label: the name if one was configured, otherwise the mass delta
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.name {
            Some(name) => name.fmt(f),
            None => write!(f, "{:+}", self.mass),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ModificationError {
    TooMany,
    /// `max_count` only applies to variable modifications
    StaticMaxCount(ModificationSpecificity),
}

impl Display for ModificationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ModificationError::TooMany => write!(
                f,
                "at most {} distinct modifications are supported",
                u8::MAX
            ),
            ModificationError::StaticMaxCount(specificity) => write!(
                f,
                "static modification {} cannot have a max_count",
                specificity
            ),
        }
    }
}

impl std::error::Error for ModificationError {}

/// Every distinct (mass, name) modification used by a database
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Modifications {
    entries: Vec<Modification>,
}

impl Modifications {
    /// Return the id for this modification, registering it if it is new
    pub fn insert(
        &mut self,
        mass: f32,
        name: Option<&ModName>,
    ) -> Result<ModId, ModificationError> {
        let index = match self
            .entries
            .iter()
            .position(|m| m.mass == mass && m.name.as_ref() == name)
        {
            Some(index) => index,
            None if self.entries.len() >= u8::MAX as usize => {
                return Err(ModificationError::TooMany)
            }
            None => {
                self.entries.push(Modification {
                    mass,
                    name: name.cloned(),
                });
                self.entries.len() - 1
            }
        };

        // `index < u8::MAX`, so `index + 1` is a valid non-zero u8
        Ok(ModId(NonZeroU8::new(index as u8 + 1).unwrap()))
    }

    pub fn get(&self, id: ModId) -> &Modification {
        &self.entries[id.index()]
    }

    /// Mass delta for a (possibly unmodified) site
    pub fn mass(&self, id: Option<ModId>) -> f32 {
        id.map_or(0.0, |id| self.get(id).mass)
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Modification> {
        self.entries.iter()
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

impl Display for InvalidModification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InvalidModification::Empty => write!(f, "invalid modification string: empty"),
            InvalidModification::InvalidResidue(c) => write!(
                f,
                "invalid modification string: unrecognized residue ({})",
                c
            ),
            InvalidModification::TooLong(s) => {
                write!(f, "invalid modification string: {} is too long", s)
            }
        }
    }
}

impl std::error::Error for InvalidModification {}

impl FromStr for ModificationSpecificity {
    type Err = InvalidModification;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut chars = s.chars();
        let first = chars.next().ok_or(InvalidModification::Empty)?;
        let rest = chars.as_str();

        let residue = |c: char| match c.is_ascii() && VALID_AA.contains(&(c as u8)) {
            true => Ok(c as u8),
            false => Err(InvalidModification::InvalidResidue(c)),
        };
        // An optional residue following a terminal prefix, e.g. `^M`
        let terminal_residue = || match rest.chars().collect::<Vec<_>>()[..] {
            [] => Ok(None),
            [c] => residue(c).map(Some),
            _ => Err(InvalidModification::TooLong(s.into())),
        };

        match first {
            '^' => terminal_residue().map(ModificationSpecificity::PeptideN),
            '$' => terminal_residue().map(ModificationSpecificity::PeptideC),
            '[' => terminal_residue().map(ModificationSpecificity::ProteinN),
            ']' => terminal_residue().map(ModificationSpecificity::ProteinC),
            c if rest.is_empty() => residue(c).map(ModificationSpecificity::Residue),
            _ => Err(InvalidModification::TooLong(s.into())),
        }
    }
}

impl<'de> Deserialize<'de> for ModificationSpecificity {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(de::Error::custom)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashMap;

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
        assert_eq!("^".parse::<ModificationSpecificity>(), Ok(PeptideN(None)));
        assert_eq!(
            "$K".parse::<ModificationSpecificity>(),
            Ok(PeptideC(Some(b'K')))
        );
        assert_eq!("".parse::<ModificationSpecificity>(), Err(Empty));
        // Residues after a terminal prefix are validated too
        assert_eq!(
            "^Z".parse::<ModificationSpecificity>(),
            Err(InvalidResidue('Z'))
        );
        assert_eq!(
            "MM".parse::<ModificationSpecificity>(),
            Err(TooLong("MM".into()))
        );
        assert_eq!(
            "^MM".parse::<ModificationSpecificity>(),
            Err(TooLong("^MM".into()))
        );
        assert_eq!(
            "^é".parse::<ModificationSpecificity>(),
            Err(InvalidResidue('é'))
        );
    }

    #[test]
    fn var_mod_entry_bare_mass() {
        let entry = ModEntry::Mass(15.9949);
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.max_count(), None);
    }

    #[test]
    fn var_mod_entry_detailed_with_limit() {
        let entry = ModEntry::Detailed(DetailedMod {
            mass: 15.9949,
            max_count: Some(1),
            name: None,
        });
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.max_count(), Some(1));
    }

    #[test]
    fn deserialize_var_mod_entries() {
        let entries: Vec<ModEntry> = serde_json::from_str(
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

        let round_tripped: Vec<ModEntry> = serde_json::from_value(serialized).unwrap();
        assert_eq!(entries, round_tripped);
    }

    #[test]
    fn reject_unsupported_var_mod_shapes_and_fields() {
        assert!(serde_json::from_str::<Vec<ModEntry>>(r#"[[42.0106, 1]]"#).is_err());
        assert!(serde_json::from_str::<Vec<ModEntry>>(
            r#"[{"mass": 42.0106, "max_count": 1, "label": "Acetyl"}]"#
        )
        .is_err());
        assert!(
            serde_json::from_str::<Vec<ModEntry>>(r#"[{"mass": 42.0106, "max_counts": 1}]"#)
                .is_err()
        );
    }

    #[test]
    fn deserialize_mod_map() {
        use ModificationSpecificity::*;
        let mods: HashMap<ModificationSpecificity, Vec<ModEntry>> = serde_json::from_str(
            r#"{"M": [15.9949, {"mass": 15.9949, "max_count": 1}], "^": [42.0106]}"#,
        )
        .unwrap();
        assert_eq!(mods[&Residue(b'M')].len(), 2);
        assert_eq!(mods[&Residue(b'M')][1].max_count(), Some(1));
        assert_eq!(mods[&PeptideN(None)], vec![ModEntry::Mass(42.0106)]);

        // Invalid keys are rejected while parsing the config, not skipped
        for key in ["Z", "^Z", "MM", ""] {
            let json = format!(r#"{{"{}": [15.9949], "M": [15.9949]}}"#, key);
            let err =
                serde_json::from_str::<HashMap<ModificationSpecificity, Vec<ModEntry>>>(&json)
                    .unwrap_err();
            assert!(
                err.to_string().contains("invalid modification string"),
                "{}",
                err
            );
        }
    }

    #[test]
    fn deserialize_named_mods() {
        let var: Vec<ModEntry> = serde_json::from_str(
            r#"[{"mass": 15.9949, "name": "UNIMOD:35", "max_count": 2}, 15.9949]"#,
        )
        .unwrap();
        assert_eq!(var[0].name().map(ModName::as_str), Some("UNIMOD:35"));
        assert_eq!(var[0].max_count(), Some(2));
        assert_eq!(var[1].name(), None);

        let fixed: HashMap<String, ModEntry> =
            serde_json::from_str(r#"{"C": 57.0215, "^": {"mass": 229.1629, "name": "TMT6plex"}}"#)
                .unwrap();
        assert_eq!(fixed["C"], ModEntry::Mass(57.0215));
        assert_eq!(fixed["^"].name().map(ModName::as_str), Some("TMT6plex"));
        assert!((fixed["^"].mass() - 229.1629).abs() < 1e-4);

        // Bare masses still serialize as numbers
        let serialized = serde_json::to_value(&fixed).unwrap();
        assert!(serialized["C"].is_number());
        assert_eq!(serialized["^"]["name"], "TMT6plex");
    }

    #[test]
    fn registry_dedups_and_labels() {
        let mut mods = Modifications::default();
        let name: ModName = "UNIMOD:35".parse().unwrap();
        let ox = mods.insert(15.9949, Some(&name)).unwrap();
        let ox_unnamed = mods.insert(15.9949, None).unwrap();
        let ph = mods.insert(79.9663, None).unwrap();
        assert_eq!(mods.insert(15.9949, Some(&name)), Ok(ox));
        assert_ne!(ox, ox_unnamed);
        assert_eq!(mods.len(), 3);

        assert_eq!(mods.get(ox).to_string(), "UNIMOD:35");
        assert_eq!(mods.get(ph).to_string(), "+79.9663");
        assert_eq!(mods.mass(Some(ph)), 79.9663);
        assert_eq!(mods.mass(None), 0.0);

        assert_eq!(std::mem::size_of::<Option<ModId>>(), 1);
    }

    #[test]
    fn reject_bad_names() {
        for name in ["", "UNIMOD:35]", "[x"] {
            assert_eq!(name.parse::<ModName>(), Err(InvalidModName(name.into())));
            let json = format!(r#"{{"mass": 1.0, "name": "{}"}}"#, name);
            let err = serde_json::from_str::<ModEntry>(&json).unwrap_err();
            assert!(
                err.to_string().contains("invalid modification name"),
                "{}",
                err
            );
        }
    }

    #[test]
    fn registry_capacity() {
        let mut mods = Modifications::default();
        for i in 0..u8::MAX {
            mods.insert(i as f32, None).unwrap();
        }
        assert_eq!(mods.insert(1000.0, None), Err(ModificationError::TooMany));
        // Existing entries can still be looked up once full
        assert!(mods.insert(1.0, None).is_ok());
    }
}
