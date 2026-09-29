// banded_exhaustive.rs - Exhaustive verification of the certified banded
// aligner over EVERY pair of sequences up to a given length over a given
// alphabet, for all scoring presets and every band width.
//
// For each pair (q, r) and preset it checks, against parasail's original
// production function (nw_trace_striped_sat) as the reference:
//   * the full-width band (which needs no certificate) equals the reference,
//     for all three implementations (SIMD, scalar, row-major reference);
//   * for every narrower band, the three implementations agree on whether
//     the band is certified, and every certified result equals the reference;
//   * align_certified (the production entry point) equals the reference;
//   * parasail's scan 16-bit kernel (production fallback) equals the reference;
//   * the lower bound never exceeds the optimal score.
//
//   cargo run --release --example banded_exhaustive -- <alphabet> <max_len>
//   e.g. ACGTN 5    ACGT 6    ACGTa 4

use cgdist::core::alignment::compute_alignment_stats;
use cgdist::core::banded::{align_certified, verify, BandedStats, Scoring};
use parasail_rs::{Aligner, Matrix};
use rayon::prelude::*;
use std::cell::RefCell;
use std::sync::atomic::{AtomicU64, Ordering};

const PRESETS: [Scoring; 3] = [
    Scoring {
        match_score: 2,
        mismatch: -1,
        gap_open: 5,
        gap_extend: 2,
    },
    Scoring {
        match_score: 3,
        mismatch: -2,
        gap_open: 8,
        gap_extend: 3,
    },
    Scoring {
        match_score: 1,
        mismatch: 0,
        gap_open: 3,
        gap_extend: 1,
    },
];

thread_local! {
    static AL: RefCell<Option<Vec<(Aligner, Aligner)>>> = const { RefCell::new(None) };
}

fn build(s: &Scoring, scan16: bool) -> Aligner {
    let m = Matrix::create(b"ACGT", s.match_score, s.mismatch).unwrap();
    let mut b = Aligner::new();
    b.matrix(m)
        .gap_open(s.gap_open)
        .gap_extend(s.gap_extend)
        .global()
        .use_trace();
    if scan16 {
        b.scan().solution_width(16);
    }
    b.build()
}

fn stats(a: &Aligner, q: &[u8], r: &[u8]) -> BandedStats {
    let res = a.align(Some(q), r).unwrap();
    assert!(!res.is_saturated());
    let tb = res.get_traceback_strings(q, r).unwrap();
    let (snps, indel_events, indel_bases) = compute_alignment_stats(&tb.query, &tb.reference);
    BandedStats {
        snps,
        indel_events,
        indel_bases,
        score: res.get_score(),
    }
}

fn all_seqs(alpha: &[u8], max_len: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut cur: Vec<Vec<u8>> = vec![vec![]];
    for _ in 0..max_len {
        let mut next = Vec::new();
        for s in &cur {
            for &c in alpha {
                let mut t = s.clone();
                t.push(c);
                next.push(t);
            }
        }
        out.extend(next.iter().cloned());
        cur = next;
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let alpha = args
        .get(1)
        .map(|s| s.as_bytes().to_vec())
        .unwrap_or(b"ACGTN".to_vec());
    let max_len: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let seqs = all_seqs(&alpha, max_len);
    let n = seqs.len();
    eprintln!(
        "alphabet {} max_len {}: {} sequences, {} ordered pairs x {} presets",
        String::from_utf8_lossy(&alpha),
        max_len,
        n,
        n * n,
        PRESETS.len()
    );
    let checks = AtomicU64::new(0);
    let certified_narrow = AtomicU64::new(0);
    let failures = AtomicU64::new(0);
    let fallbacks = AtomicU64::new(0);
    (0..n).into_par_iter().for_each(|qi| {
        AL.with(|cell| {
            let mut al = cell.borrow_mut();
            let al = al.get_or_insert_with(|| {
                PRESETS
                    .iter()
                    .map(|s| (build(s, false), build(s, true)))
                    .collect()
            });
            let q = &seqs[qi];
            let mut local_checks = 0u64;
            let mut local_cert = 0u64;
            for r in &seqs {
                for (pi, s) in PRESETS.iter().enumerate() {
                    let want = stats(&al[pi].0, q, r);
                    let fail = |what: &str| {
                        failures.fetch_add(1, Ordering::Relaxed);
                        eprintln!(
                            "FAIL {what}: preset {pi} q={} r={}",
                            String::from_utf8_lossy(q),
                            String::from_utf8_lossy(r)
                        );
                    };
                    if stats(&al[pi].1, q, r) != want {
                        fail("scan16 != striped_sat");
                    }
                    // None = handed to parasail (always allowed); Some must be exact
                    match align_certified(q, r, s, 1.0) {
                        Some(got) if got != want => fail("align_certified wrong"),
                        Some(_) => {}
                        None => {
                            fallbacks.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    if verify::lower_bound(q, r, s) > want.score as i64 {
                        fail("lower bound above optimum");
                    }
                    let wmax = q.len().max(r.len()) as i64;
                    for w in 0..=wmax {
                        let res = verify::band_results(q, r, s, w);
                        local_checks += 1;
                        if res[0] != res[1] || res[0] != res[2] || res[0] != res[3] {
                            fail(&format!("implementations disagree at w={w}"));
                        }
                        match res[0] {
                            Some(got) => {
                                if got != want {
                                    fail(&format!("certified result wrong at w={w}"));
                                }
                                if w < wmax {
                                    local_cert += 1;
                                }
                            }
                            None => {
                                if w == wmax {
                                    fail("full band not certified");
                                }
                            }
                        }
                    }
                }
            }
            checks.fetch_add(local_checks, Ordering::Relaxed);
            certified_narrow.fetch_add(local_cert, Ordering::Relaxed);
        })
    });
    let f = failures.load(Ordering::Relaxed);
    println!(
        "alphabet {} max_len {}: pairs {} x presets {} ; band checks {} ; certified narrow bands {} ; align_certified fallbacks {} ; FAILURES {}",
        String::from_utf8_lossy(&alpha),
        max_len,
        n * n,
        PRESETS.len(),
        checks.load(Ordering::Relaxed),
        certified_narrow.load(Ordering::Relaxed),
        fallbacks.load(Ordering::Relaxed),
        f
    );
    if f > 0 {
        std::process::exit(1);
    }
}
