// align_bench.rs - Speed and exactness benchmark for the pairwise alignment step
//
// Input: TSV with locus, seq1, seq2, snps, indel_events, indel_bases, score as
// produced by `cgdist --save-alignments` (reference = current production path).
// Every variant must reproduce (snps, indel_events, indel_bases) exactly.
//
//   cargo run --release --example align_bench -- <pairs.tsv> [max_pairs] [timing_pairs]

use cgdist::core::alignment::compute_alignment_stats;
use parasail_rs::{Aligner, Matrix};
use rayon::prelude::*;
use std::cell::RefCell;
use std::io::{BufRead, BufReader};
use std::time::Instant;

struct Pair {
    q: Vec<u8>,
    r: Vec<u8>,
    gold: (usize, usize, usize),
}

const MATCH: i32 = 2;
const MISMATCH: i32 = -1;
const GAP_OPEN: i32 = 5;
const GAP_EXTEND: i32 = 2;

fn build(width: Option<i32>, strategy: &str) -> Aligner {
    let matrix = Matrix::create(b"ACGT", MATCH, MISMATCH).unwrap();
    let mut b = Aligner::new();
    b.matrix(matrix)
        .gap_open(GAP_OPEN)
        .gap_extend(GAP_EXTEND)
        .global()
        .use_trace();
    match strategy {
        "scan" => {
            b.scan();
        }
        "diag" => {
            b.diag();
        }
        _ => {
            b.striped();
        }
    }
    if let Some(w) = width {
        b.solution_width(w);
    }
    b.build()
}

/// Stats from a parasail CIGAR ("=", "X", "I", "D" runs); an I run directly
/// followed by a D run is one event, as in compute_alignment_stats.
fn stats_from_cigar(cigar: &str) -> (usize, usize, usize) {
    let (mut snps, mut ev, mut bases, mut n, mut in_gap) = (0, 0, 0, 0usize, false);
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as usize;
            continue;
        }
        let len = if n == 0 { 1 } else { n };
        n = 0;
        match c {
            b'I' | b'D' => {
                if !in_gap {
                    ev += 1;
                    in_gap = true;
                }
                bases += len;
            }
            b'X' => {
                in_gap = false;
                snps += len;
            }
            _ => in_gap = false,
        }
    }
    (snps, ev, bases)
}

fn current(p: &Pair) -> (usize, usize, usize) {
    // Mirrors DistanceEngine::compute_single_alignment exactly.
    let a = build(None, "striped");
    let res = a.align(Some(&p.q), &p.r).unwrap();
    let tb = res.get_traceback_strings(&p.q, &p.r).unwrap();
    compute_alignment_stats(&tb.query, &tb.reference)
}

thread_local! {
    static AL: RefCell<Vec<(String, Aligner)>> = const { RefCell::new(Vec::new()) };
}

fn with_aligner<T>(key: &str, width: Option<i32>, strat: &str, f: impl FnOnce(&Aligner) -> T) -> T {
    AL.with(|c| {
        let mut v = c.borrow_mut();
        if !v.iter().any(|(k, _)| k == key) {
            v.push((key.to_string(), build(width, strat)));
        }
        let a = &v.iter().find(|(k, _)| k == key).unwrap().1;
        f(a)
    })
}

fn variant(name: &str, p: &Pair) -> (usize, usize, usize) {
    let (width, strat, cigar) = match name {
        "reuse_sat" => (None, "striped", false),
        "reuse_16" => (Some(16), "striped", false),
        "reuse_16_scan" => (Some(16), "scan", false),
        "reuse_16_diag" => (Some(16), "diag", false),
        "reuse_16_cigar" => (Some(16), "striped", true),
        "reuse_sat_scan" => (None, "scan", false),
        "reuse_32_scan" => (Some(32), "scan", false),
        "reuse_16_scan_cigar" => (Some(16), "scan", true),
        "scan16_dp_only" => (Some(16), "scan", false),
        "reuse_sat_cigar" => (None, "striped", true),
        _ => unreachable!(),
    };
    with_aligner(name, width, strat, |a| {
        let res = a.align(Some(&p.q), &p.r).unwrap();
        if name == "scan16_dp_only" {
            return p.gold; // timing of the DP + trace table only
        }
        if cigar {
            stats_from_cigar(&res.get_cigar(&p.q, &p.r).unwrap())
        } else {
            let tb = res.get_traceback_strings(&p.q, &p.r).unwrap();
            compute_alignment_stats(&tb.query, &tb.reference)
        }
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let max: usize = args
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);
    let timing: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5000);

    let f = std::fs::File::open(path).unwrap();
    let pairs: Vec<Pair> = BufReader::new(f)
        .lines()
        .take(max)
        .map(|l| {
            let l = l.unwrap();
            let c: Vec<&str> = l.split('\t').collect();
            Pair {
                q: c[1].as_bytes().to_vec(),
                r: c[2].as_bytes().to_vec(),
                gold: (
                    c[3].parse().unwrap(),
                    c[4].parse().unwrap(),
                    c[5].parse().unwrap(),
                ),
            }
        })
        .collect();
    eprintln!("loaded {} pairs from {path}", pairs.len());
    let tp = &pairs[..timing.min(pairs.len())];

    // Single-thread timing on a fixed prefix, exactness on all pairs (parallel).
    let t = Instant::now();
    let n_ok = tp.iter().filter(|p| current(p) == p.gold).count();
    let base = t.elapsed().as_secs_f64();
    println!(
        "{:<18} {:>9.1} us/pair  1.00x  timing-set exact {}/{}",
        "current",
        base / tp.len() as f64 * 1e6,
        n_ok,
        tp.len()
    );

    let names: Vec<String> = match std::env::var("VARIANTS") {
        Ok(v) => v.split(',').map(String::from).collect(),
        Err(_) => [
            "reuse_sat",
            "reuse_sat_cigar",
            "reuse_16",
            "reuse_16_cigar",
            "reuse_16_scan",
            "reuse_16_diag",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect(),
    };
    for name in names.iter().map(String::as_str) {
        let t = Instant::now();
        for p in tp {
            std::hint::black_box(variant(name, p));
        }
        let el = t.elapsed().as_secs_f64();
        let mism: Vec<usize> = pairs
            .par_iter()
            .enumerate()
            .filter(|(_, p)| variant(name, p) != p.gold)
            .map(|(i, _)| i)
            .collect();
        println!(
            "{:<18} {:>9.1} us/pair  {:>4.2}x  mismatches {}/{}{}",
            name,
            el / tp.len() as f64 * 1e6,
            base / el,
            mism.len(),
            pairs.len(),
            if mism.is_empty() {
                String::new()
            } else {
                format!("  first at line {:?}", &mism[..mism.len().min(5)])
            }
        );
    }
}
