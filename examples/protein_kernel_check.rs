// protein_kernel_check.rs - cgdist's protein alignment (scan 16->32->64 bit)
// vs parasail's striped saturating kernel, on real proteins translated from
// a schema: every pair of distinct proteins of the first N loci (capped per
// locus) must give identical substitutions / InDel events / InDel residues
// and score.
//
//   cargo run --release --example protein_kernel_check -- <schema_dir> [loci] [pairs_per_locus] [matrix]

use bio::io::fasta;
use cgdist::core::alignment::compute_alignment_stats;
use cgdist::core::protein::GeneticCode;
use cgdist::core::protein_distance::{align_proteins, protein_of, ProteinSettings};
use parasail_rs::{Aligner, Matrix};
use rayon::prelude::*;
use std::collections::BTreeSet;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = &args[1];
    let n_loci: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(200);
    let per_locus: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(300);
    let matrix = args.get(4).cloned().unwrap_or_else(|| "blosum62".into());
    let settings = ProteinSettings {
        translation_table: 11,
        first_codon_as_met: true,
        matrix: matrix.clone(),
        gap_open: 11,
        gap_extend: 1,
    };
    let code = GeneticCode::default();
    let mut files: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "fasta"))
        .collect();
    files.sort();
    let mut pairs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for f in files.iter().take(n_loci) {
        let prots: BTreeSet<Vec<u8>> = fasta::Reader::from_file(f)
            .unwrap()
            .records()
            .map(|r| protein_of(&code, r.unwrap().seq()))
            .collect();
        let v: Vec<Vec<u8>> = prots.into_iter().collect();
        let mut k = 0;
        'outer: for i in 0..v.len() {
            for j in i + 1..v.len() {
                pairs.push((v[i].clone(), v[j].clone()));
                k += 1;
                if k >= per_locus {
                    break 'outer;
                }
            }
        }
    }
    // TIMING: cgdist protein alignment alone, single thread
    let t = std::time::Instant::now();
    let take = pairs.len().min(20_000);
    for (a, b) in &pairs[..take] {
        std::hint::black_box(align_proteins(&settings, a, b));
    }
    let per = t.elapsed().as_secs_f64() / take as f64 * 1e6;
    let mean_len: f64 = pairs[..take]
        .iter()
        .map(|(a, b)| (a.len() + b.len()) as f64 / 2.0)
        .sum::<f64>()
        / take as f64;
    println!("TIMING: {per:.1} us per protein pair (mean length {mean_len:.0} aa, single thread)");
    let bad = pairs
        .par_iter()
        .filter(|(a, b)| {
            let got = align_proteins(&settings, a, b).unwrap();
            // reference: parasail striped, saturating (aligner kept alive)
            let m = Matrix::from(&matrix).unwrap();
            let al = Aligner::new()
                .matrix(m)
                .gap_open(11)
                .gap_extend(1)
                .global()
                .use_trace()
                .build();
            let res = al.align(Some(a), b).unwrap();
            let tb = res.get_traceback_strings(a, b).unwrap();
            let (s, e, r) = compute_alignment_stats(&tb.query, &tb.reference);
            (
                got.aa_subs as usize,
                got.aa_indel_events as usize,
                got.aa_indel_residues as usize,
            ) != (s, e, r)
        })
        .count();
    println!(
        "{matrix}: {} protein pairs from {} loci of {dir}: mismatches vs parasail striped_sat {bad}",
        pairs.len(),
        files.len().min(n_loci)
    );
    if bad > 0 {
        std::process::exit(1);
    }
}
