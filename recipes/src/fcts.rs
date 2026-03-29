use std::f64::consts;

const NCOF: usize = 28;
const COF: [f64; 28] = [
    -1.3026537197817094,
    6.4196979235649026e-1,
    1.9476473204185836e-2,
    -9.561514786808631e-3,
    -9.46595344482036e-4,
    3.66839497852761e-4,
    4.2523324806907e-5,
    -2.0278578112534e-5,
    -1.624290004647e-6,
    1.303655835580e-6,
    1.5626441722e-8,
    -8.5238095915e-8,
    6.529054439e-9,
    5.059343495e-9,
    -9.91364156e-10,
    -2.27365122e-10,
    9.6467911e-11,
    2.394038e-12,
    -6.886027e-12,
    8.94487e-13,
    3.13092e-13,
    -1.12708e-13,
    3.81e-16,
    7.106e-15,
    -1.523e-15,
    -9.4e-17,
    1.21e-16,
    -2.8e-17,
];

pub fn erf(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 - erfccheb(x)
    } else {
        erfccheb(-x) - 1.0
    }
}

pub fn erfc(x: f64) -> f64 {
    if x >= 0.0 {
        erfccheb(x)
    } else {
        2.0 - erfccheb(-x)
    }
}

fn erfccheb(z: f64) -> f64 {
    if z < 0.0 {
        panic!("erfccheb requires nonnegative argument");
    }
    let t = 2.0 / (2.0 + z);
    let ty = 4.0 * t - 2.0;
    let mut d = 0.0;
    let mut dd = 0.0;
    for j in (1..NCOF).rev() {
        let tmp = d;
        d = ty * d - dd + COF[j];
        dd = tmp;
    }
    t * f64::exp(-z * z + 0.5 * (COF[0] + ty * d) - dd)
}

pub fn inverfc(p: f64) -> f64 {
    if p >= 2.0 {
        return -100.0;
    }
    if p <= 0.0 {
        return 100.0;
    }
    let pp = if p < 1.0 { p } else { 2.0 - p };
    let t = f64::sqrt(-2.0 * f64::ln(pp / 2.0));
    let mut x = -0.70711 * ((2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t);
    for _ in 0..2 {
        let err = erfc(x) - pp;
        x += err / (1.12837916709551257 * f64::exp(-(x * x)) - x * err);
    }
    if p < 1.0 { x } else { -x }
}

pub fn inverf(p: f64) -> f64 {
    inverfc(1.0 - p)
}

pub fn gaussian_cumulative(x: f64) -> f64 {
    0.5 * (1.0 + erf(x * consts::FRAC_1_SQRT_2))
}

pub fn gaussian_cumulative_inverse(p: f64) -> f64 {
    consts::SQRT_2 * inverf(2.0 * p - 1.0)
}

#[cfg(test)]
mod tests {
    use crate::fcts::{erf, erfc, inverf, gaussian_cumulative, gaussian_cumulative_inverse};

    #[test]
    fn test_erf_zero() {
        assert!((erf(0.0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_erfc_zero() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_erf_symmetry() {
        assert!((erf(-1.0) + erf(1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_inverf_erf() {
        let x0 = 0.5;
        let p = erf(x0);
        let x1 = inverf(p);
        assert!((x0 - x1).abs() < 1e-6);
    }

    #[test]
    fn test_inv_gaussian_cumulative() {
        let x0 = 0.5;
        let p = gaussian_cumulative(x0);
        let x1 = gaussian_cumulative_inverse(p);
        assert!((x0 - x1).abs() < 1e-6);
    }

    #[test]
    fn test_gaussian_cumulative_3_stddev() {
        let expected = [
            ( 0.0, 0.0 ),
            ( 0.5, 38.2924922548026 ),
            ( 1.0, 68.2689492137086 ),
            ( 2.0, 95.4499736103642 ),
            ( 3.0, 99.7300203936740 ),
            ];

        for xp in expected {
            let x = xp.0;
            let p = 1.0 - 2.0*gaussian_cumulative(-x);
            let q = xp.1;
            assert!((p * 100.0 - q).abs() < 5e-14);
        }
    }
}