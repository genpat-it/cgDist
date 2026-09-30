// Where does --coding-stats spend its time? Times, on the allele pairs of
// some schema loci: banded stats only; banded + gapped strings; + string
// statistics; + coding classification.
//   cargo run --release --example coding_bench -- <schema_dir> <n_loci>
use cgdist::core::alignment::compute_alignment_stats;
use cgdist::core::banded::{align_certified, align_certified_with_strings, Scoring};
use cgdist::core::protein::{coding_counts, GeneticCode};
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut files: Vec<_> = std::fs::read_dir(&a[1])
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "fasta"))
        .collect();
    files.sort();
    let n: usize = a[2].parse().unwrap();
    let mut pairs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for f in files.iter().take(n) {
        let seqs: Vec<Vec<u8>> = bio::io::fasta::Reader::from_file(f)
            .unwrap()
            .records()
            .map(|r| r.unwrap().seq().to_vec())
            .take(120)
            .collect();
        for i in 0..seqs.len() {
            for j in i + 1..seqs.len() {
                pairs.push((seqs[i].clone(), seqs[j].clone()));
            }
        }
    }
    let s = Scoring {
        match_score: 2,
        mismatch: -1,
        gap_open: 5,
        gap_extend: 2,
    };
    let code = GeneticCode::default();
    println!("{} pairs", pairs.len());
    let t = Instant::now();
    let mut x = 0usize;
    for (p, q) in &pairs {
        if let Some(b) = align_certified(p, q, &s, 0.5) {
            x += b.snps;
        }
    }
    let t1 = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let mut strs = Vec::with_capacity(pairs.len());
    for (p, q) in &pairs {
        strs.push(align_certified_with_strings(p, q, &s, 0.5).map(|(_, st)| st));
    }
    let t2 = t.elapsed().as_secs_f64();
    let t = Instant::now();
    for st in strs.iter().flatten() {
        let q = String::from_utf8_lossy(&st.query);
        let r = String::from_utf8_lossy(&st.reference);
        x += compute_alignment_stats(&q, &r).0;
    }
    let t3 = t.elapsed().as_secs_f64();
    let t = Instant::now();
    for ((p, q), st) in pairs.iter().zip(&strs) {
        if let Some(st) = st {
            x += coding_counts(&code, p, q, &st.query, &st.reference).syn as usize;
        }
    }
    let t4 = t.elapsed().as_secs_f64();
    let us = |v: f64| v * 1e6 / pairs.len() as f64;
    println!("banded stats      {:8.2} us/pair", us(t1));
    println!("banded + strings  {:8.2} us/pair", us(t2));
    println!("string stats      {:8.2} us/pair", us(t3));
    println!("coding_counts     {:8.2} us/pair", us(t4));
    println!("(checksum {x})");
}
