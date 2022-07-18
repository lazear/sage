pub const H2O: f32 = 18.010565;
pub const PROTON: f32 = 1.007_276_4;
pub const NH3: f32 = 17.026548;

#[derive(Copy, Clone)]
pub enum Tolerance {
    Ppm(f32),
    Th(f32),
}

impl Tolerance {
    pub fn bounds(&self, center: f32) -> (f32, f32) {
        match self {
            Tolerance::Ppm(ppm) => {
                let delta = center * ppm / 1_000_000.0;
                (center - delta, center + delta)
            }
            Tolerance::Th(th) => (center - th, center + th),
        }
    }
}

pub trait Mass {
    fn monoisotopic(&self) -> f32;
}

#[derive(Hash, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Residue {
    Ala,
    Arg,
    Asn,
    Asp,
    Cys,
    Glu,
    Gln,
    Gly,
    His,
    Ile,
    Leu,
    Lys,
    Met,
    Phe,
    Pro,
    Ser,
    Thr,
    Trp,
    Tyr,
    Val,
    Mod(Box<Residue>, Modification),
}

#[derive(Hash, Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Modification {
    Carbamidomethyl,
    Oxidation,
    Acetylation,
}

impl Mass for Modification {
    fn monoisotopic(&self) -> f32 {
        match self {
            Modification::Carbamidomethyl => 57.021464,
            Modification::Oxidation => 15.994915,
            Modification::Acetylation => 42.010565,
        }
    }
}

impl Mass for Residue {
    fn monoisotopic(&self) -> f32 {
        use Residue::*;
        match self {
            Ala => 71.037_12,
            Arg => 156.101_1,
            Asn => 114.042_93,
            Asp => 115.026_94,
            Cys => 103.009_186,
            Glu => 129.042_59,
            Gln => 128.058_58,
            Gly => 57.021_465,
            His => 137.058_91,
            Ile => 113.084_06,
            Leu => 113.084_06,
            Lys => 128.094_96,
            Met => 131.040_48,
            Phe => 147.068_42,
            Pro => 97.052_765,
            Ser => 87.032_03,
            Thr => 101.047_676,
            Trp => 186.079_32,
            Tyr => 163.063_32,
            Val => 99.068_41,
            Mod(resi, modi) => resi.monoisotopic() + modi.monoisotopic(),
        }
    }
}

impl TryFrom<char> for Residue {
    type Error = char;
    fn try_from(value: char) -> Result<Self, Self::Error> {
        use Residue::*;
        let r = match value {
            'A' => Ala,
            'R' => Arg,
            'N' => Asn,
            'D' => Asp,
            'C' => Cys,
            'E' => Glu,
            'Q' => Gln,
            'G' => Gly,
            'H' => His,
            'I' => Ile,
            'L' => Leu,
            'K' => Lys,
            'M' => Met,
            'F' => Phe,
            'P' => Pro,
            'S' => Ser,
            'T' => Thr,
            'W' => Trp,
            'Y' => Tyr,
            'V' => Val,
            _ => return Err(value),
        };
        Ok(r)
    }
}

impl std::fmt::Display for Residue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use Residue::*;
        match self {
            Ala => write!(f, "A"),
            Arg => write!(f, "R"),
            Asn => write!(f, "N"),
            Asp => write!(f, "D"),
            Cys => write!(f, "C"),
            Glu => write!(f, "E"),
            Gln => write!(f, "Q"),
            Gly => write!(f, "G"),
            His => write!(f, "H"),
            Ile => write!(f, "I"),
            Leu => write!(f, "L"),
            Lys => write!(f, "K"),
            Met => write!(f, "M"),
            Phe => write!(f, "F"),
            Pro => write!(f, "P"),
            Ser => write!(f, "S"),
            Thr => write!(f, "T"),
            Trp => write!(f, "W"),
            Tyr => write!(f, "Y"),
            Val => write!(f, "V"),
            Mod(r, m) => write!(f, "{}({})", r, m.monoisotopic()),
        }
    }
}
