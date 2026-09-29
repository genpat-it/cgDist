// banded_fuzz.rs - Adversarial exactness test: certified banded alignment vs
// parasail (the reference), on randomly generated hard cases: tandem repeats
// and homopolymers with indels (many co-optimal alignments), N / lowercase
// symbols, tiny sequences, large length differences, unrelated sequences, and
// all three scoring presets. Any certified result must equal parasail's
// (SNPs, InDel events, InDel bases, score) exactly.
//
//   cargo run --release --example banded_fuzz -- [cases] [seed]

use cgdist::core::alignment::compute_alignment_stats;
use cgdist::core::banded::{align_certified, Scoring};
use parasail_rs::{Aligner, Matrix};
use rayon::prelude::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    fn chance(&mut self, p_per_mille: u64) -> bool {
        self.below(1000) < p_per_mille
    }
}

fn base(rng: &mut Rng) -> u8 {
    b"ACGT"[rng.below(4) as usize]
}

fn make_seq(rng: &mut Rng, len: usize) -> Vec<u8> {
    let mut v = Vec::with_capacity(len);
    while v.len() < len {
        match rng.below(10) {
            // tandem repeat unit 1..6 repeated
            0..=2 => {
                let unit: Vec<u8> = (0..1 + rng.below(6)).map(|_| base(rng)).collect();
                for _ in 0..2 + rng.below(12) {
                    v.extend_from_slice(&unit);
                }
            }
            3 => {
                let b = base(rng);
                for _ in 0..3 + rng.below(15) {
                    v.push(b);
                }
            }
            _ => {
                for _ in 0..1 + rng.below(30) {
                    v.push(base(rng));
                }
            }
        }
    }
    v.truncate(len);
    v
}

fn mutate(rng: &mut Rng, s: &[u8], rate_pm: u64) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len() + 16);
    let mut i = 0;
    while i < s.len() {
        if rng.chance(rate_pm) {
            match rng.below(6) {
                0 | 1 => out.push(base(rng)), // substitution
                2 => {
                    // insertion
                    out.push(s[i]);
                    for _ in 0..1 + rng.below(8) {
                        out.push(base(rng));
                    }
                }
                3 => i += rng.below(8) as usize, // deletion
                4 => out.push(if rng.chance(500) {
                    b'N'
                } else {
                    s[i].to_ascii_lowercase()
                }),
                _ => {
                    // duplicate a short segment (repeat expansion)
                    let l = 1 + rng.below(6) as usize;
                    let e = (i + l).min(s.len());
                    out.extend_from_slice(&s[i..e]);
                    out.push(s[i]);
                }
            }
        } else {
            out.push(s[i]);
        }
        i += 1;
    }
    if out.is_empty() {
        out.push(base(rng));
    }
    out
}

fn reference(q: &[u8], r: &[u8], sc: &Scoring) -> ((usize, usize, usize), i32) {
    let m = Matrix::create(b"ACGT", sc.match_score, sc.mismatch).unwrap();
    let a = Aligner::new()
        .matrix(m)
        .gap_open(sc.gap_open)
        .gap_extend(sc.gap_extend)
        .global()
        .use_trace()
        .build(); // striped, saturating: the original production kernel
    let res = a.align(Some(q), r).unwrap();
    let tb = res.get_traceback_strings(q, r).unwrap();
    (
        compute_alignment_stats(&tb.query, &tb.reference),
        res.get_score(),
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cases: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(200_000);
    let seed: u64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
    let presets = [
        (
            "dna",
            Scoring {
                match_score: 2,
                mismatch: -1,
                gap_open: 5,
                gap_extend: 2,
            },
        ),
        (
            "dna-strict",
            Scoring {
                match_score: 3,
                mismatch: -2,
                gap_open: 8,
                gap_extend: 3,
            },
        ),
        (
            "dna-permissive",
            Scoring {
                match_score: 1,
                mismatch: 0,
                gap_open: 3,
                gap_extend: 1,
            },
        ),
    ];
    let results: Vec<(usize, bool, bool)> = (0..cases)
        .into_par_iter()
        .map(|c| {
            let mut rng = Rng(seed.wrapping_mul(1_000_003).wrapping_add(c));
            let preset = (c % 3) as usize;
            let sc = presets[preset].1;
            let len = match rng.below(10) {
                0 => 1 + rng.below(12) as usize,
                1 => 20 + rng.below(80) as usize,
                2 => 1500 + rng.below(3000) as usize,
                _ => 200 + rng.below(1300) as usize,
            };
            let a = make_seq(&mut rng, len);
            let b = match rng.below(10) {
                0 => {
                    // unrelated
                    let l = 1 + rng.below(len as u64 * 2) as usize;
                    make_seq(&mut rng, l)
                }
                1 => {
                    // large length difference
                    let mut t = mutate(&mut rng, &a, 20);
                    let cut = rng.below(t.len() as u64 / 2 + 1) as usize;
                    if rng.chance(500) {
                        t.drain(..cut);
                    } else {
                        let extra = make_seq(&mut rng, cut + 1);
                        t.extend_from_slice(&extra);
                    }
                    if t.is_empty() {
                        t.push(b'A');
                    }
                    t
                }
                x => mutate(
                    &mut rng,
                    &a,
                    [2, 5, 10, 20, 40, 80, 150, 300][(x - 2) as usize],
                ),
            };
            let (q, r) = if rng.chance(500) { (a, b) } else { (b, a) };
            let fraction = if rng.chance(100) { 1.0 } else { 0.5 };
            match align_certified(&q, &r, &sc, fraction) {
                None => (preset, false, true),
                Some(res) => {
                    let got = ((res.snps, res.indel_events, res.indel_bases), res.score);
                    let want = reference(&q, &r, &sc);
                    if got != want {
                        eprintln!(
                            "MISMATCH preset={} lens={}/{} got={:?} want={:?}\nq={}\nr={}",
                            presets[preset].0,
                            q.len(),
                            r.len(),
                            got,
                            want,
                            String::from_utf8_lossy(&q),
                            String::from_utf8_lossy(&r)
                        );
                    }
                    (preset, got == want, false)
                }
            }
        })
        .collect();
    let mut wrong = 0;
    for (pi, (name, _)) in presets.iter().enumerate() {
        let tot = results.iter().filter(|r| r.0 == pi).count();
        let unc = results.iter().filter(|r| r.0 == pi && r.2).count();
        let ok = results.iter().filter(|r| r.0 == pi && r.1).count();
        let bad = tot - unc - ok;
        wrong += bad;
        println!("{name:<15} cases {tot:>8}  certified&identical {ok:>8}  uncertified(fallback) {unc:>7}  WRONG {bad}");
    }
    println!("TOTAL WRONG: {wrong}");
    if wrong > 0 {
        std::process::exit(1);
    }
}
