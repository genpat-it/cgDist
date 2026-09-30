// Heavy loci (long, divergent alleles): how often does the certified band
// give up at max_fraction 0.5, and what do the fallbacks cost?
//   cargo run --release --example heavy_locus_bench -- <locus.fasta> <n_pairs>
use cgdist::core::banded::{align_certified_with_strings, Scoring};
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let seqs: Vec<Vec<u8>> = bio::io::fasta::Reader::from_file(&a[1])
        .unwrap()
        .records()
        .map(|r| r.unwrap().seq().to_vec())
        .collect();
    let n: usize = a[2].parse().unwrap();
    let s = Scoring {
        match_score: 2,
        mismatch: -1,
        gap_open: 5,
        gap_extend: 2,
    };
    let cfg = cgdist::core::alignment::AlignmentConfig::from_mode("dna").unwrap();
    let mut seed = 7u64;
    let mut pick = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) as usize % seqs.len()
    };
    let pairs: Vec<(usize, usize)> = (0..n)
        .map(|_| (pick(), pick()))
        .filter(|(i, j)| i != j)
        .collect();
    let (mut fail, mut t_band, mut t_par, mut t_full_band) = (0usize, 0f64, 0f64, 0f64);
    let mut same = 0usize;
    let mut checked_ok = 0usize;
    for &(i, j) in &pairs {
        let (q, r) = (&seqs[i], &seqs[j]);
        let t = Instant::now();
        let b = align_certified_with_strings(q, r, &s, 1.0);
        t_band += t.elapsed().as_secs_f64();
        if let Some((st, strs)) = &b {
            // every certified result must equal the full alignment
            let p = cgdist::core::distance::align_pair_with_strings(&cfg, q, r).unwrap();
            let full = cgdist::core::parasail_trace::NwTracer::new(
                "nw_trace_scan_32",
                parasail_rs::Matrix::create(b"ACGT", 2, -1).unwrap(),
                5,
                2,
            )
            .unwrap()
            .align(q, r)
            .unwrap();
            if strs.query != full.query.as_bytes()
                || strs.reference != full.reference.as_bytes()
                || (st.snps, st.indel_events, st.indel_bases)
                    != (p.snps, p.indel_events, p.indel_bases)
            {
                panic!("adaptive band differs from the full alignment on pair {i}/{j}");
            }
            checked_ok += 1;
            continue;
        }
        fail += 1;
        let t = Instant::now();
        let p = cgdist::core::distance::align_pair_with_strings(&cfg, q, r).unwrap();
        t_par += t.elapsed().as_secs_f64();
        let t = Instant::now();
        let fb = align_certified_with_strings(q, r, &s, 1.0);
        t_full_band += t.elapsed().as_secs_f64();
        if let Some((st, strs)) = fb {
            if (st.snps, st.indel_events, st.indel_bases) == (p.snps, p.indel_events, p.indel_bases)
                && strs.query == p.query.as_bytes()
                && strs.reference == p.reference.as_bytes()
            {
                same += 1;
            }
        }
    }
    let k = pairs.len();
    println!("certified results identical to the full parasail alignment: {checked_ok}");
    println!(
        "{k} pairs, mean length {:.0}",
        pairs.iter().map(|&(i, _)| seqs[i].len()).sum::<usize>() as f64 / k as f64
    );
    println!(
        "band (0.5) attempts: {:.2} ms/pair; gave up on {fail} ({:.1}%)",
        1e3 * t_band / k as f64,
        100.0 * fail as f64 / k as f64
    );
    if fail > 0 {
        println!(
            "fallback now (parasail path): {:.1} ms/pair",
            1e3 * t_par / fail as f64
        );
        println!("certified band at max_fraction 1.0: {:.1} ms/pair, identical to the parasail path on {same}/{fail}", 1e3 * t_full_band / fail as f64);
    }
}
