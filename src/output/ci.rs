// ci.rs - Missingness confidence interval for pairwise distances (opt-in).
//
// A cgDist distance rests only on the loci a pair shares; the unobserved loci
// would add more distance. Under MCAR missingness we model the unobserved
// distance M and report a normalized point estimate plus an interval.
//
// Decomposition (validated in cgPosterior/poc):
//   over the S shared loci we have d (distance), h (# differing loci),
//   q2 (sum of squared per-locus contributions); U = total - S unobserved.
//   K ~ BetaBinomial(U, 0.5+h, 0.5+(S-h))          differing unobserved loci
//   each contributes ~ (mean r=d/h, var s2=q2/h - r^2)
//   D = d + M
// The interval is the UNION of the exact Beta-Binomial count interval (carries
// the discrete small-distance / outbreak regime) and the normal compound
// interval (carries the contribution variance in the mid/large regime). This
// deterministic union matches a gold Monte-Carlo posterior predictive across
// regimes. Where both fail (few differing loci + high missingness) coverage is
// information-limited, not fixable — flagged via `reliable = false`.
//
// Self-contained: ln_gamma (Lanczos) and normal_ppf (Acklam) implemented here so
// no numeric dependency is added to the published tool.

/// Lanczos approximation of ln Γ(x) for x > 0.
// Coefficients are kept verbatim from the published tables.
#[allow(clippy::excessive_precision)]
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // reflection: Γ(x)Γ(1-x) = π / sin(πx)
        std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + G + 0.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Inverse standard-normal CDF (Acklam's rational approximation), |err| < 1.2e-9.
#[allow(clippy::excessive_precision)]
fn normal_ppf(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1, 2.209_460_984_245_205e2, -2.759_285_104_469_687e2,
        1.383_577_518_672_690e2, -3.066_479_806_614_716e1, 2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1, 1.615_858_368_580_409e2, -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1, -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3, -3.223_964_580_411_365e-1, -2.400_758_277_161_838,
        -2.549_732_539_343_734, 4.374_664_141_464_968, 2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3, 3.224_671_290_700_398e-1, 2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    let plow = 0.024_25;
    let phigh = 1.0 - plow;
    if p < plow {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= phigh {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// Lower/upper quantiles (as counts) of BetaBinomial(u, a, b) at central `level`.
/// Iterates the PMF via its closed-form consecutive ratio (one ln_gamma-heavy
/// evaluation for pmf(0), then cheap multiplications).
fn betabinom_interval(u: u64, a: f64, b: f64, level: f64) -> (f64, f64) {
    if u == 0 {
        return (0.0, 0.0);
    }
    let uf = u as f64;
    let qlo = (1.0 - level) / 2.0;
    let qhi = 1.0 - qlo;
    // ln pmf(0) = lnΓ(u+b) + lnΓ(a+b) - lnΓ(a+u+b) - lnΓ(b)
    let ln_p0 = ln_gamma(uf + b) + ln_gamma(a + b) - ln_gamma(a + uf + b) - ln_gamma(b);
    let mut p = ln_p0.exp();
    let mut cdf = p;
    let mut klo: Option<f64> = if cdf >= qlo { Some(0.0) } else { None };
    let mut khi = uf;
    for k in 0..u {
        let kf = k as f64;
        // ratio pmf(k+1)/pmf(k) = (u-k)/(k+1) * (k+a)/(u-k-1+b)
        let ratio = ((uf - kf) / (kf + 1.0)) * ((kf + a) / (uf - kf - 1.0 + b));
        p *= ratio;
        cdf += p;
        if klo.is_none() && cdf >= qlo {
            klo = Some(kf + 1.0);
        }
        if cdf >= qhi {
            khi = kf + 1.0;
            break;
        }
    }
    (klo.unwrap_or(0.0), khi)
}

/// Per-pair missingness CI. Returns (dist_norm, ci_low, ci_high, reliable).
/// `d` = observed distance, `h` = # differing shared loci, `q2` = sum of squared
/// per-locus contributions over shared loci, `shared` = co-present loci,
/// `total` = schema loci, `level` e.g. 0.95.
pub fn pair_ci(
    d: usize,
    h: usize,
    q2: u64,
    shared: usize,
    total: usize,
    level: f64,
) -> (f64, f64, f64, bool) {
    let d = d as f64;
    if h == 0 || shared == 0 || total <= shared {
        return (d, d, d, true);
    }
    let u = (total - shared) as u64;
    let uf = u as f64;
    let hf = h as f64;
    let sf = shared as f64;
    let a = 0.5 + hf;
    let b = 0.5 + (sf - hf);
    let r = d / hf;
    let s2 = (q2 as f64 / hf - r * r).max(0.0);
    // exact Beta-Binomial count interval
    let (klo, khi) = betabinom_interval(u, a, b, level);
    let bb_lo = d + r * klo;
    let bb_hi = d + r * khi;
    // normal compound interval
    let ek = uf * a / (a + b);
    let var_k = uf * a * b * (a + b + uf) / ((a + b).powi(2) * (a + b + 1.0));
    let em = r * ek;
    let var_m = ek * s2 + var_k * r * r;
    let z = normal_ppf(0.5 + level / 2.0);
    let se = var_m.max(0.0).sqrt();
    let norm_lo = d + em - z * se;
    let norm_hi = d + em + z * se;
    // union
    let lo = d.max(bb_lo.min(norm_lo));
    let hi = bb_hi.max(norm_hi);
    let dist_norm = d + em;
    // reliability: coverage degrades when few differing loci coincide with high
    // missingness (information limit). Conservative heuristic flag.
    let missing_frac = 1.0 - sf / total as f64;
    let reliable = !(h < 10 && missing_frac > 0.15);
    (dist_norm, lo, hi, reliable)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ln_gamma_known_values() {
        // Γ(5)=24 -> ln24; Γ(0.5)=sqrt(pi) -> 0.5*ln(pi)
        assert!((ln_gamma(5.0) - 24.0_f64.ln()).abs() < 1e-9);
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-9);
    }

    #[test]
    fn normal_ppf_known_values() {
        assert!((normal_ppf(0.5)).abs() < 1e-9);
        assert!((normal_ppf(0.975) - 1.959_963_984_540_054).abs() < 1e-6);
    }

    #[test]
    fn ci_brackets_and_orders() {
        // full coverage -> point interval
        let (dn, lo, hi, rel) = pair_ci(10, 4, 30, 100, 100, 0.95);
        assert_eq!((lo, hi), (dn, dn));
        assert!(rel);
        // with missingness the interval widens above d and stays >= d
        let (dn2, lo2, hi2, _) = pair_ci(10, 4, 30, 80, 100, 0.95);
        assert!(lo2 >= 10.0 - 1e-9 && hi2 > dn2 - 1e-9 && dn2 >= 10.0);
        assert!(hi2 > lo2);
    }

    #[test]
    fn low_info_flagged_unreliable() {
        // few differing loci + high missingness -> flagged
        let (_, _, _, rel) = pair_ci(6, 3, 14, 900, 1748, 0.95);
        assert!(!rel);
    }
}
