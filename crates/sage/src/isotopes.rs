// Manually unroll the convolution loop
fn convolve(a: &[f32; 4], b: &[f32; 4]) -> [f32; 4] {
    [
        a[0] * b[0],
        a[0] * b[1] + a[1] * b[0],
        a[0] * b[2] + a[1] * b[1] + a[2] * b[0],
        a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0],
        // a[0] * b[4] + a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + a[4] * b[0],
    ]
}

fn carbon_isotopes(count: u16) -> [f32; 4] {
    let lambda = count as f32 * 0.011;
    let mut c13 = [0.0; 4];

    let fact = [1, 1, 2, 6];
    for k in 0..4 {
        c13[k] = lambda.powi(k as i32) * f32::exp(-lambda) / fact[k] as f32;
    }
    c13
}

fn sulfur_isotopes(count: u16) -> [f32; 4] {
    let lambda33 = count as f32 * 0.0076;
    let lambda35 = count as f32 * 0.044;
    let mut s33 = [0.0; 4];
    let s35 = [
        lambda35.powi(0) * f32::exp(-lambda35),
        0.0,
        lambda35.powi(1) * f32::exp(-lambda35),
        0.0,
        // lambda35.powi(2) * f32::exp(-lambda35) / 2.0,
    ];

    let fact = [1, 1, 2, 6];
    for k in 0..4 {
        s33[k] = lambda33.powi(k as i32) * f32::exp(-lambda33) / fact[k] as f32;
    }

    convolve(&s33, &s35)
}

pub fn peptide_isotopes(carbons: u16, sulfurs: u16) -> [f32; 4] {
    let c = carbon_isotopes(carbons);
    let s = sulfur_isotopes(sulfurs);
    let mut c = convolve(&c, &s);
    let max = c[0].max(c[1]).max(c[2]);
    c.iter_mut().for_each(|val| *val /= max);
    c
}

#[cfg(test)]
mod tests {
    use super::peptide_isotopes;

    #[test]
    fn smoke_isotopes() {
        let iso = peptide_isotopes(60, 5);
        let mut expected = [0.3972, 0.2824, 0.1869, 0.0846];
        expected.iter_mut().for_each(|val| *val /= 0.3972);

        let matched = iso.iter().zip(expected).all(|(a, b)| (a - b).abs() <= 0.02);

        assert!(matched, "{:?} {:?}", iso, expected);
    }
}
