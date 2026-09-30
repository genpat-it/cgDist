// Exactness of the adaptive probe path of the certified aligner, which only
// runs on long, divergent pairs: random 2-8 kb alleles with 2-25% SNPs and
// InDels of assorted lengths, three scoring presets; every certified result
// (statistics, score, gapped strings) must equal parasail's reference kernel.
//   cargo run --release --example long_divergent_fuzz -- [pairs] [seed]
use cgdist::core::alignment::compute_alignment_stats;
use cgdist::core::banded::{align_certified_with_strings, Scoring};
use cgdist::core::parasail_trace::NwTracer;
use parasail_rs::Matrix;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let pairs: usize = a.get(1).map_or(600, |x| x.parse().unwrap());
    let mut rng = Rng(a.get(2).map_or(11, |x| x.parse().unwrap()));
    let presets = [(2, -1, 5, 2), (3, -2, 8, 3), (1, 0, 3, 1)];
    let (mut ok, mut wrong, mut uncert, mut probed) = (0, 0, 0, 0);
    for k in 0..pairs {
        let (ms, mm, go, ge) = presets[k % 3];
        let len = 2000 + rng.below(6000) as usize;
        let q: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
        let snp_rate = 2 + rng.below(24); // percent
        let mut r = Vec::with_capacity(len + 100);
        let mut i = 0;
        while i < q.len() {
            let x = rng.below(1000);
            if x < snp_rate * 10 {
                r.push(b"ACGT"[((q[i] as u64 + 1 + rng.below(3)) % 4) as usize]);
                i += 1;
            } else if x < snp_rate * 10 + 3 {
                // deletion of 1-30
                i += 1 + rng.below(30) as usize;
            } else if x < snp_rate * 10 + 6 {
                // insertion of 1-30
                for _ in 0..1 + rng.below(30) {
                    r.push(b"ACGT"[rng.below(4) as usize]);
                }
            } else {
                r.push(q[i]);
                i += 1;
            }
        }
        let s = Scoring {
            match_score: ms,
            mismatch: mm,
            gap_open: go,
            gap_extend: ge,
        };
        let reference = NwTracer::new(
            "nw_trace_striped_sat",
            Matrix::create(b"ACGT", ms, mm).unwrap(),
            go,
            ge,
        )
        .unwrap()
        .align(&q, &r)
        .unwrap();
        match align_certified_with_strings(&q, &r, &s, 1.0) {
            None => uncert += 1,
            Some((b, st)) => {
                if b.score > 0 && (b.indel_events > 0 || b.snps > len / 30) {
                    probed += 1;
                }
                let want = compute_alignment_stats(&reference.query, &reference.reference);
                let same = (b.snps, b.indel_events, b.indel_bases) == want
                    && b.score == reference.score
                    && st.query == reference.query.as_bytes()
                    && st.reference == reference.reference.as_bytes();
                if same {
                    ok += 1
                } else {
                    wrong += 1;
                    eprintln!(
                        "MISMATCH pair {k} len {len} snp {snp_rate}% preset {:?}",
                        presets[k % 3]
                    );
                }
            }
        }
    }
    println!("{pairs} long divergent pairs: certified & identical {ok}, WRONG {wrong}, uncertified (parasail) {uncert}, divergent {probed}");
}
