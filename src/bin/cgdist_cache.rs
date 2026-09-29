// cgdist_cache.rs - Build, convert, publish and download cgdist cache stores
//
//   cgdist-cache build  --schema DIR --out STORE      all pairs of every locus
//   cgdist-cache import --cache FILE.lz4 --out STORE  convert a cgdist cache
//   cgdist-cache pack   --store STORE --out FILE.cgpack
//   cgdist-cache pull   --from SRC --out STORE [--profiles P | --loci-list L]
//   cgdist-cache info   --store SRC
//   cgdist-cache verify --store SRC
//
// SRC is a store directory, a .cgpack file, or an http(s) URL of either.
// `build` aligns with cgdist's own engine (same certified banded alignment,
// parasail fallback and checks), so stores hold exactly what cgdist would
// compute. It is incremental and resumable: only missing pairs are aligned,
// and the manifest is saved after every locus.

use argh::FromArgs;
use cgdist::core::alignment::{AlignmentConfig, DistanceMode};
use cgdist::core::distance::{DistanceEngine, GeneticCodeMeta, ModernCache};
use cgdist::core::protein::GeneticCode;
use cgdist::data::{SequenceDatabase, SequenceInfo};
use cgdist::store::remote::{self, Source};
use cgdist::store::{AlignmentParams, LocusData, PairStats, Store};
use std::collections::{BTreeMap, HashSet};
use std::path::Path;
use std::time::Instant;

#[derive(FromArgs)]
/// cgdist-cache: build, convert, publish and download cgdist cache stores.
struct Cli {
    #[argh(subcommand)]
    cmd: Cmd,
}

#[derive(FromArgs)]
#[argh(subcommand)]
enum Cmd {
    Build(Build),
    Import(Import),
    Pack(Pack),
    Pull(Pull),
    Info(Info),
    Verify(Verify),
}

#[derive(FromArgs)]
/// align every pair of alleles of every locus of a schema into a store
#[argh(subcommand, name = "build")]
struct Build {
    /// schema directory (one <locus>.fasta per locus)
    #[argh(option)]
    schema: String,
    /// store directory (created if missing; existing pairs are kept)
    #[argh(option)]
    out: String,
    /// alignment mode: dna, dna-strict, dna-permissive (default: dna)
    #[argh(option, default = "String::from(\"dna\")")]
    alignment_mode: String,
    /// custom match score (with --mismatch-penalty, --gap-open, --gap-extend)
    #[argh(option)]
    match_score: Option<i32>,
    /// custom mismatch penalty
    #[argh(option)]
    mismatch_penalty: Option<i32>,
    /// custom gap open penalty
    #[argh(option)]
    gap_open: Option<i32>,
    /// custom gap extend penalty
    #[argh(option)]
    gap_extend: Option<i32>,
    /// only these loci (file, one per line)
    #[argh(option)]
    loci_list: Option<String>,
    /// threads (default: 1; 0 = all cores)
    #[argh(option, default = "1")]
    threads: usize,
    /// also store synonymous/nonsynonymous counts
    #[argh(switch)]
    coding_stats: bool,
    /// NCBI translation table for --coding-stats (default: 11)
    #[argh(option, default = "11")]
    translation_table: u32,
    /// do not read an alternative start codon at codon 1 as Met
    #[argh(switch)]
    no_first_codon_as_met: bool,
    /// pairs aligned per batch within a locus (default: 2000000)
    #[argh(option, default = "2_000_000")]
    batch_pairs: usize,
    /// re-check this fraction (0-1) of alignments against parasail's
    /// original kernel; any difference stops the build (default: 0)
    #[argh(option, default = "0.0")]
    verify_alignments: f64,
    /// schema name recorded in the manifest
    #[argh(option)]
    schema_name: Option<String>,
    /// schema source recorded in the manifest (e.g. chewie-ns)
    #[argh(option)]
    schema_source: Option<String>,
    /// schema version recorded in the manifest
    #[argh(option)]
    schema_version: Option<String>,
}

#[derive(FromArgs)]
/// convert a cgdist .lz4 cache into a store
#[argh(subcommand, name = "import")]
struct Import {
    /// cgdist cache file (.lz4, format v2)
    #[argh(option)]
    cache: String,
    /// store directory
    #[argh(option)]
    out: String,
}

#[derive(FromArgs)]
/// write a store as a single .cgpack file for publishing
#[argh(subcommand, name = "pack")]
struct Pack {
    /// store directory
    #[argh(option)]
    store: String,
    /// output .cgpack file
    #[argh(option)]
    out: String,
}

#[derive(FromArgs)]
/// download a published store (all loci, or only those of some profiles)
#[argh(subcommand, name = "pull")]
struct Pull {
    /// store directory, .cgpack file, or http(s) URL of either
    #[argh(option)]
    from: String,
    /// local store directory
    #[argh(option)]
    out: String,
    /// only the loci in the header of this profiles file (TSV/CSV)
    #[argh(option)]
    profiles: Option<String>,
    /// only these loci (file, one per line)
    #[argh(option)]
    loci_list: Option<String>,
    /// download threads (default: 8)
    #[argh(option, default = "8")]
    threads: usize,
}

#[derive(FromArgs)]
/// summary of a store
#[argh(subcommand, name = "info")]
struct Info {
    /// store directory, .cgpack file, or URL
    #[argh(option)]
    store: String,
}

#[derive(FromArgs)]
/// check every locus of a store (checksum, decoding, counts)
#[argh(subcommand, name = "verify")]
struct Verify {
    /// store directory, .cgpack file, or URL
    #[argh(option)]
    store: String,
}

fn threads(n: usize) {
    if n > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global();
    }
}

fn read_list(path: &str) -> Result<HashSet<String>, String> {
    Ok(std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {path}: {e}"))?
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect())
}

fn loci_of_profiles(path: &str) -> Result<HashSet<String>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let header = text
        .lines()
        .find(|l| !l.starts_with('#') && !l.trim().is_empty())
        .ok_or("empty profiles file")?;
    let sep = if header.contains('\t') { '\t' } else { ',' };
    Ok(header
        .split(sep)
        .skip(1)
        .map(|s| s.trim().to_string())
        .collect())
}

fn build(a: Build) -> Result<(), String> {
    threads(a.threads);
    let custom = [a.match_score, a.mismatch_penalty, a.gap_open, a.gap_extend];
    let config = if custom.iter().all(Option::is_some) {
        AlignmentConfig::custom(
            custom[0].unwrap(),
            custom[1].unwrap(),
            custom[2].unwrap(),
            custom[3].unwrap(),
        )
    } else if custom.iter().any(Option::is_some) {
        return Err("custom scoring needs all of --match-score, --mismatch-penalty, --gap-open, --gap-extend".into());
    } else {
        AlignmentConfig::from_mode(&a.alignment_mode)?
    };
    let code = if a.coding_stats {
        Some(GeneticCode::new(
            a.translation_table,
            !a.no_first_codon_as_met,
        )?)
    } else {
        None
    };
    let only = a.loci_list.as_deref().map(read_list).transpose()?;
    let mut store =
        Store::open_or_create(Path::new(&a.out), "crc32", AlignmentParams::from(&config))?;
    if let Some(c) = &code {
        let meta = GeneticCodeMeta::from(c);
        match store.manifest.genetic_code {
            Some(m) if m != meta => {
                return Err(format!(
                    "store {} holds coding counts for another genetic code",
                    a.out
                ))
            }
            _ => store.manifest.genetic_code = Some(meta),
        }
    }
    for (field, val) in [
        (&mut store.manifest.schema.name, &a.schema_name),
        (&mut store.manifest.schema.source, &a.schema_source),
        (&mut store.manifest.schema.version, &a.schema_version),
    ] {
        if val.is_some() {
            *field = val.clone();
        }
    }
    let _lock = store.lock()?;

    let mut files: Vec<_> = std::fs::read_dir(&a.schema)
        .map_err(|e| format!("cannot read schema {}: {e}", a.schema))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "fasta" || x == "fa"))
        .collect();
    files.sort();
    let start = Instant::now();
    let mut last_save = Instant::now();
    let (mut done_pairs, mut total_new) = (0u64, 0u64);
    let n_files = files.len();
    for (fi, path) in files.iter().enumerate() {
        let locus = path.file_stem().unwrap().to_string_lossy().to_string();
        if only.as_ref().is_some_and(|o| !o.contains(&locus)) {
            continue;
        }
        // alleles of the locus, keyed by CRC32 of the sequence as in cgdist
        let reader = bio::io::fasta::Reader::from_file(path)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let mut alleles: BTreeMap<u32, (String, Vec<u8>)> = BTreeMap::new();
        let mut collision = false;
        for rec in reader.records() {
            let rec = rec.map_err(|e| format!("{}: {e}", path.display()))?;
            let mut h = crc32fast::Hasher::new();
            h.update(rec.seq());
            let crc = h.finalize();
            match alleles.get(&crc) {
                Some((_, s)) if s.as_slice() != rec.seq() => collision = true,
                Some(_) => {}
                None => {
                    alleles.insert(crc, (rec.id().to_string(), rec.seq().to_vec()));
                }
            }
        }
        if collision {
            eprintln!("⚠️  {locus}: two different alleles share a CRC32; locus skipped");
            continue;
        }
        let mut data = store.read_locus(&locus)?.unwrap_or_default();
        for (crc, (_, seq)) in &alleles {
            data.set_allele_len(*crc, seq.len() as u32);
        }
        let crcs: Vec<u32> = alleles.keys().copied().collect();
        let n = crcs.len();
        let needs_coding = |st: &PairStats| code.is_some() && st.coding.is_none();
        // missing pairs (or pairs lacking requested coding counts)
        let mut missing: Vec<(u32, u32)> = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                match data.pairs.get(&(crcs[i], crcs[j])) {
                    Some(st) if !needs_coding(st) => {}
                    _ => missing.push((crcs[i], crcs[j])),
                }
            }
        }
        if !missing.is_empty() {
            let mut db = SequenceDatabase::new();
            for (crc, (id, seq)) in &alleles {
                db.add_sequence(
                    locus.clone(),
                    *crc,
                    SequenceInfo {
                        sequence: seq.clone(),
                        id: id.clone(),
                    },
                );
            }
            let mut engine =
                DistanceEngine::with_sequences(config.clone(), db, "crc32".to_string());
            engine.set_quiet(true);
            engine.set_verify_fraction(a.verify_alignments);
            if let Some(c) = code {
                engine.set_coding(c);
            }
            for batch in missing.chunks(a.batch_pairs.max(1)) {
                let set: HashSet<(String, u32, u32)> =
                    batch.iter().map(|&(x, y)| (locus.clone(), x, y)).collect();
                engine.precompute_alignments(&set, DistanceMode::SnpsOnly);
                for (_, d) in engine.drain_to_locus_data() {
                    data.merge(d);
                }
            }
            total_new += missing.len() as u64;
        }
        done_pairs += (n * n.saturating_sub(1) / 2) as u64;
        store.write_locus(&locus, &data)?;
        // the manifest is the resume point: save it every ~30 s (a locus
        // file written after the last save is simply rebuilt on resume)
        if last_save.elapsed().as_secs() >= 30 {
            store.save_manifest()?;
            last_save = Instant::now();
        }
        let el = start.elapsed().as_secs_f64();
        println!(
            "[{}/{}] {locus}: {n} alleles, {} new pairs ({}complete)  | {total_new} new pairs in {el:.0}s",
            fi + 1,
            n_files,
            missing.len(),
            if data.is_complete() { "" } else { "in" }
        );
    }
    store.save_manifest()?;
    println!(
        "✅ store {}: {} loci, {done_pairs} pairs covered, {total_new} newly aligned in {:.1}s",
        a.out,
        store.manifest.loci.len(),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn import(a: Import) -> Result<(), String> {
    let bytes = std::fs::read(&a.cache).map_err(|e| format!("cannot read {}: {e}", a.cache))?;
    let raw = lz4_flex::decompress_size_prepended(&bytes)
        .map_err(|e| format!("cannot decompress {}: {e}", a.cache))?;
    let cache: ModernCache = serde_json::from_slice(&raw).map_err(|e| {
        format!(
            "{} is not a cgdist v2 cache (run cgdist once to convert legacy caches): {e}",
            a.cache
        )
    })?;
    let params = AlignmentParams::from(&cache.metadata.alignment_config);
    let mut store = Store::open_or_create(Path::new(&a.out), &cache.metadata.hasher_type, params)?;
    let _lock = store.lock()?;
    let has_coding = cache.data.values().any(|v| v.syn.is_some());
    if has_coding {
        match (store.manifest.genetic_code, cache.metadata.genetic_code) {
            (Some(x), Some(y)) if x != y => {
                return Err("store and cache use different genetic codes".into())
            }
            (_, g) => store.manifest.genetic_code = g.or(store.manifest.genetic_code),
        }
    }
    let mut by_locus: BTreeMap<String, LocusData> = BTreeMap::new();
    for (k, v) in &cache.data {
        let mut it = k.rsplitn(3, ':');
        let (Some(b), Some(x), Some(locus)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let (Ok(x), Ok(b)) = (x.parse::<u32>(), b.parse::<u32>()) else {
            continue;
        };
        let d = by_locus.entry(locus.to_string()).or_default();
        d.insert_pair(
            x,
            b,
            PairStats {
                snps: v.snps as u32,
                indel_events: v.indel_events as u32,
                indel_bases: v.indel_bases as u32,
                coding: match (
                    v.syn,
                    v.nonsyn,
                    v.frame_disrupted,
                    cache.metadata.genetic_code,
                ) {
                    (Some(s), Some(n), Some(f), Some(_)) => Some((s, n, f)),
                    _ => None,
                },
            },
        );
        if let Some(l) = v.seq1_length {
            d.set_allele_len(x.min(b), l as u32);
        }
        if let Some(l) = v.seq2_length {
            d.set_allele_len(x.max(b), l as u32);
        }
    }
    let n = by_locus.len();
    let mut pairs = 0usize;
    for (locus, d) in by_locus {
        let mut cur = store.read_locus(&locus)?.unwrap_or_default();
        pairs += d.pairs.len();
        cur.merge(d);
        store.write_locus(&locus, &cur)?;
    }
    store.save_manifest()?;
    println!(
        "✅ imported {pairs} pairs in {n} loci from {} into {}",
        a.cache, a.out
    );
    Ok(())
}

fn info(a: Info) -> Result<(), String> {
    let src = Source::parse(&a.store);
    let m = src.manifest()?;
    let (alleles, pairs, bytes, complete) =
        m.loci
            .values()
            .fold((0usize, 0usize, 0u64, 0usize), |acc, e| {
                (
                    acc.0 + e.alleles,
                    acc.1 + e.pairs,
                    acc.2 + e.bytes,
                    acc.3 + usize::from(e.complete),
                )
            });
    println!("store        {}", src.describe());
    println!(
        "format       {} v{} (cgdist {})",
        m.format, m.format_version, m.cgdist_version
    );
    println!("hasher       {}", m.hasher);
    println!("alignment    {}", m.alignment);
    println!(
        "coding       {}",
        m.genetic_code
            .map(|g| format!(
                "NCBI table {}, first codon as Met: {}",
                g.table, g.first_codon_as_met
            ))
            .unwrap_or_else(|| "none".into())
    );
    println!(
        "schema       {} {} {}",
        m.schema.name.as_deref().unwrap_or("-"),
        m.schema.source.as_deref().unwrap_or("-"),
        m.schema.version.as_deref().unwrap_or("-")
    );
    println!(
        "loci         {} ({complete} complete: every pair of their alleles)",
        m.loci.len()
    );
    println!("alleles      {alleles}");
    println!("pairs        {pairs}");
    println!("size         {:.1} MB", bytes as f64 / 1e6);
    println!("modified     {}", m.last_modified);
    Ok(())
}

fn verify(a: Verify) -> Result<(), String> {
    let src = Source::parse(&a.store);
    let m = src.manifest()?;
    let mut problems = 0usize;
    for (locus, e) in &m.loci {
        match src.locus_bytes(e).and_then(|b| LocusData::decode(&b)) {
            Ok(d) => {
                if d.alleles.len() != e.alleles
                    || d.pairs.len() != e.pairs
                    || d.is_complete() != e.complete
                {
                    eprintln!("✗ {locus}: counts differ from the manifest");
                    problems += 1;
                }
            }
            Err(err) => {
                eprintln!("✗ {locus}: {err}");
                problems += 1;
            }
        }
    }
    if problems > 0 {
        return Err(format!(
            "{problems} of {} loci failed verification",
            m.loci.len()
        ));
    }
    println!("✅ {}: all {} loci verified", src.describe(), m.loci.len());
    Ok(())
}

fn main() {
    let cli: Cli = argh::from_env();
    let res = match cli.cmd {
        Cmd::Build(a) => build(a),
        Cmd::Import(a) => import(a),
        Cmd::Pack(a) => remote::pack(Path::new(&a.store), Path::new(&a.out))
            .map(|n| println!("✅ packed {n} loci into {}", a.out)),
        Cmd::Pull(a) => {
            threads(a.threads);
            let filter = match (&a.profiles, &a.loci_list) {
                (Some(p), _) => Some(loci_of_profiles(p)),
                (None, Some(l)) => Some(read_list(l)),
                (None, None) => None,
            }
            .transpose();
            filter.and_then(|f| {
                let src = Source::parse(&a.from);
                let st = remote::pull(&src, Path::new(&a.out), f.as_ref())?;
                println!(
                    "✅ pulled from {}: {} loci fetched ({:.1} MB), {} up to date, {} not in source",
                    src.describe(),
                    st.fetched,
                    st.bytes as f64 / 1e6,
                    st.up_to_date,
                    st.not_in_source
                );
                Ok(())
            })
        }
        Cmd::Info(a) => info(a),
        Cmd::Verify(a) => verify(a),
    };
    if let Err(e) = res {
        eprintln!("❌ ERROR: {e}");
        std::process::exit(1);
    }
}
