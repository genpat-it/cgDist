// Compare every entry of a cgdist .lz4 cache with a cache store (dir,
// .cgpack or URL): SNPs, InDel events, InDel bases and allele lengths.
//   cargo run --release --example store_compare -- cache.lz4 <store>
use cgdist::core::distance::ModernCache;
use cgdist::store::remote::Source;
use cgdist::store::LocusData;
use std::collections::HashMap;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let raw = lz4_flex::decompress_size_prepended(&std::fs::read(&a[1]).unwrap()).unwrap();
    let cache: ModernCache = serde_json::from_slice(&raw).unwrap();
    let src = Source::parse(&a[2]);
    let m = src.manifest().unwrap();
    println!(
        "cache: {} entries, params {:?}",
        cache.data.len(),
        cache.metadata.alignment_config
    );
    println!("store: {}", m.params().unwrap());
    let mut by_locus: HashMap<String, Vec<(u32, u32, &cgdist::core::distance::CacheValue)>> =
        HashMap::new();
    for (k, v) in &cache.data {
        let mut it = k.rsplitn(3, ':');
        let (b, x, l) = (it.next().unwrap(), it.next().unwrap(), it.next().unwrap());
        by_locus.entry(l.to_string()).or_default().push((
            x.parse().unwrap(),
            b.parse().unwrap(),
            v,
        ));
    }
    let (mut same, mut diff, mut missing, mut lens_diff) = (0u64, 0u64, 0u64, 0u64);
    for (locus, rows) in &by_locus {
        let Some(e) = m.loci.get(locus) else {
            missing += rows.len() as u64;
            continue;
        };
        let d = LocusData::decode(&src.locus_bytes(e).unwrap()).unwrap();
        for (x, y, v) in rows {
            let key = ((*x).min(*y), (*x).max(*y));
            match d.pairs.get(&key) {
                None => missing += 1,
                Some(s)
                    if (
                        s.snps as usize,
                        s.indel_events as usize,
                        s.indel_bases as usize,
                    ) == (v.snps, v.indel_events, v.indel_bases) =>
                {
                    same += 1;
                    let (l1, l2) = (
                        d.alleles.get(&key.0).copied(),
                        d.alleles.get(&key.1).copied(),
                    );
                    let want = if x < y {
                        (v.seq1_length, v.seq2_length)
                    } else {
                        (v.seq2_length, v.seq1_length)
                    };
                    if (l1.map(|l| l as usize), l2.map(|l| l as usize)) != want {
                        lens_diff += 1;
                    }
                }
                Some(s) => {
                    diff += 1;
                    if diff <= 5 {
                        eprintln!(
                            "DIFF {locus}:{x}:{y} cache {:?} store {:?}",
                            (v.snps, v.indel_events, v.indel_bases),
                            s
                        );
                    }
                }
            }
        }
    }
    println!("identical {same}, different {diff}, missing from store {missing}, allele-length mismatches {lens_diff}");
}
