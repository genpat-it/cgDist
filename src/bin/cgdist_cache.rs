// cgdist_cache.rs - Build, convert, publish and download cgdist cache stores
//
//   cgdist-cache build  --schema DIR --out STORE      all pairs of every locus
//   cgdist-cache build  --schema DIR --out STORE --protein   protein store
//   cgdist-cache import --cache FILE.lz4 --out STORE  convert a cgdist cache
//   cgdist-cache pack   --store STORE --out FILE.cgpack
//   cgdist-cache pull   --from SRC --out STORE [--profiles P | --loci-list L]
//   cgdist-cache info   --store SRC
//   cgdist-cache verify --store SRC
//   cgdist-cache stats  --store SRC --out stats.json  distributions for reports
//
// SRC is a store directory, a .cgpack file, or an http(s) URL of either.
// `build` aligns with cgdist's own engine (same certified banded alignment,
// parasail fallback and checks), so stores hold exactly what cgdist would
// compute. It is incremental and resumable: only missing pairs are aligned,
// and the manifest is saved after every locus.

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use argh::FromArgs;
use cgdist::core::alignment::{AlignmentConfig, DistanceMode};
use cgdist::core::distance::{DistanceEngine, GeneticCodeMeta, ModernCache};
use cgdist::core::protein::GeneticCode;
use cgdist::core::protein_distance::{self, ProteinSettings};
use cgdist::data::{SequenceDatabase, SequenceInfo};
use cgdist::store::remote::{self, Source};
use cgdist::store::{AlignmentParams, LocusData, PairStats, Store, StoreParams};
use rayon::prelude::*;
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
    Stats(Stats),
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
    /// build a protein store: translate the alleles and align every pair of
    /// distinct proteins (for aa-* modes and aa_* weights)
    #[argh(switch)]
    protein: bool,
    /// protein store: substitution matrix, blosum30..blosum100, pam10..pam500
    /// or a matrix file (default: blosum62)
    #[argh(option, default = "String::from(\"blosum62\")")]
    aa_matrix: String,
    /// protein store: gap open penalty (default: 11)
    #[argh(option, default = "11")]
    aa_gap_open: i32,
    /// protein store: gap extend penalty (default: 1)
    #[argh(option, default = "1")]
    aa_gap_extend: i32,
    /// NCBI translation table for --coding-stats and --protein (default: 11)
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
/// distributions of a store (per-locus summaries and pair histograms), JSON
#[argh(subcommand, name = "stats")]
struct Stats {
    /// store directory, .cgpack file, or URL
    #[argh(option)]
    store: String,
    /// output JSON file
    #[argh(option)]
    out: String,
    /// threads (default: 0 = all cores)
    #[argh(option, default = "0")]
    threads: usize,
}

#[derive(FromArgs)]
/// check every locus of a store (checksum, decoding, counts); with
/// --schema and --realign, also recompute a fraction of the pairs from the
/// schema sequences and compare
#[argh(subcommand, name = "verify")]
struct Verify {
    /// store directory, .cgpack file, or URL
    #[argh(option)]
    store: String,
    /// schema directory the store was built from (for --realign)
    #[argh(option)]
    schema: Option<String>,
    /// fraction (0-1) of stored pairs to recompute from the schema: DNA pairs
    /// with cgdist's engine, each also checked against parasail's original
    /// kernel; protein pairs by translating and aligning again. The choice of
    /// pairs is deterministic (default: 0)
    #[argh(option, default = "0.0")]
    realign: f64,
    /// threads for --realign (default: 0 = all cores)
    #[argh(option, default = "0")]
    threads: usize,
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
    if a.protein {
        return build_protein(a);
    }
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
    let mut store = Store::open_or_create(
        Path::new(&a.out),
        "crc32",
        StoreParams::Dna(AlignmentParams::from(&config)),
    )?;
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

    let files = schema_files(&a.schema)?;
    let start = Instant::now();
    let mut last_save = Instant::now();
    let (mut done_pairs, mut total_new) = (0u64, 0u64);
    let mut clashes = 0usize;
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
        // a store allele with this CRC32 but another sequence: its pairs
        // belong to a different allele; never mix them
        let clash = alleles.iter().find(|(crc, (_, seq))| {
            data.digests
                .get(crc)
                .is_some_and(|&d| d != cgdist::store::seq_digest(seq))
        });
        if let Some((crc, _)) = clash {
            eprintln!(
                "❌ {locus}: allele {crc} of the schema has another sequence than the allele with the \
                 same CRC32 in the store (hash collision or different schema); locus skipped"
            );
            clashes += 1;
            continue;
        }
        for (crc, (_, seq)) in &alleles {
            data.set_allele_len(*crc, seq.len() as u32);
            data.set_allele_digest(*crc, cgdist::store::seq_digest(seq));
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
    if clashes > 0 {
        return Err(format!(
            "{clashes} loci skipped: schema alleles collide with different store alleles (see above)"
        ));
    }
    Ok(())
}

/// FASTA files of a schema directory, sorted.
fn schema_files(dir: &str) -> Result<Vec<std::path::PathBuf>, String> {
    let mut files: Vec<_> = std::fs::read_dir(dir)
        .map_err(|e| format!("cannot read schema {dir}: {e}"))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "fasta" || x == "fa"))
        .collect();
    files.sort();
    Ok(files)
}

fn build_protein(a: Build) -> Result<(), String> {
    if a.coding_stats {
        return Err("--coding-stats applies to DNA stores, not with --protein".into());
    }
    let settings = ProteinSettings {
        translation_table: a.translation_table,
        first_codon_as_met: !a.no_first_codon_as_met,
        matrix: a.aa_matrix.clone(),
        gap_open: a.aa_gap_open,
        gap_extend: a.aa_gap_extend,
    }
    .normalized();
    settings.validate()?;
    let code = settings.code()?;
    let only = a.loci_list.as_deref().map(read_list).transpose()?;
    let mut store = Store::open_or_create(
        Path::new(&a.out),
        "crc32",
        StoreParams::Protein(settings.clone()),
    )?;
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
    let files = schema_files(&a.schema)?;
    let start = Instant::now();
    let mut last_save = Instant::now();
    let (mut done_pairs, mut total_new, mut failed) = (0u64, 0u64, 0u64);
    let mut clashes = 0usize;
    let n_files = files.len();
    for (fi, path) in files.iter().enumerate() {
        let locus = path.file_stem().unwrap().to_string_lossy().to_string();
        if only.as_ref().is_some_and(|o| !o.contains(&locus)) {
            continue;
        }
        let reader = bio::io::fasta::Reader::from_file(path)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        // distinct proteins of the locus, keyed by protein hash as in cgdist
        let mut proteins: BTreeMap<u32, Vec<u8>> = BTreeMap::new();
        let mut collision = false;
        for rec in reader.records() {
            let rec = rec.map_err(|e| format!("{}: {e}", path.display()))?;
            let p = protein_distance::protein_of(&code, rec.seq());
            let h = protein_distance::protein_hash(&p);
            match proteins.get(&h) {
                Some(q) if *q != p => collision = true,
                Some(_) => {}
                None => {
                    proteins.insert(h, p);
                }
            }
        }
        if collision {
            eprintln!("⚠️  {locus}: two different proteins share a CRC32; locus skipped");
            continue;
        }
        let mut data = store.read_locus(&locus)?.unwrap_or_default();
        let clash = proteins.iter().find(|(h, p)| {
            data.digests
                .get(h)
                .is_some_and(|&d| d != cgdist::store::seq_digest(p))
        });
        if let Some((h, _)) = clash {
            eprintln!(
                "❌ {locus}: protein {h} of the schema differs from the protein with the same hash \
                 in the store (hash collision or different schema); locus skipped"
            );
            clashes += 1;
            continue;
        }
        for (h, p) in &proteins {
            data.set_allele_len(*h, p.len() as u32);
            data.set_allele_digest(*h, cgdist::store::seq_digest(p));
        }
        let hs: Vec<u32> = proteins.keys().copied().collect();
        let n = hs.len();
        let mut missing: Vec<(u32, u32)> = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                if !data.pairs.contains_key(&(hs[i], hs[j])) {
                    missing.push((hs[i], hs[j]));
                }
            }
        }
        for batch in missing.chunks(a.batch_pairs.max(1)) {
            // same orientation as cgdist: lower protein hash first
            let res: Vec<(u32, u32, Option<protein_distance::ProteinPair>)> = batch
                .par_iter()
                .map(|&(x, y)| {
                    (
                        x,
                        y,
                        protein_distance::align_proteins(&settings, &proteins[&x], &proteins[&y]),
                    )
                })
                .collect();
            for (x, y, r) in res {
                match r {
                    Some(p) => data.insert_pair(x, y, protein_distance::pair_stats(&p)),
                    None => failed += 1,
                }
            }
        }
        total_new += missing.len() as u64;
        done_pairs += (n * n.saturating_sub(1) / 2) as u64;
        store.write_locus(&locus, &data)?;
        if last_save.elapsed().as_secs() >= 30 {
            store.save_manifest()?;
            last_save = Instant::now();
        }
        println!(
            "[{}/{}] {locus}: {n} proteins, {} new pairs ({}complete)  | {total_new} new pairs in {:.0}s",
            fi + 1,
            n_files,
            missing.len(),
            if data.is_complete() { "" } else { "in" },
            start.elapsed().as_secs_f64()
        );
    }
    store.save_manifest()?;
    if failed > 0 {
        eprintln!("⚠️  {failed} protein pairs could not be aligned and are not in the store");
    }
    println!(
        "✅ protein store {}: {} loci, {done_pairs} protein pairs covered, {total_new} newly aligned in {:.1}s",
        a.out,
        store.manifest.loci.len(),
        start.elapsed().as_secs_f64()
    );
    if clashes > 0 {
        return Err(format!(
            "{clashes} loci skipped: schema proteins collide with different store proteins (see above)"
        ));
    }
    Ok(())
}

/// Histogram bin upper edges (inclusive): 0..=10 one by one, then wider.
const HIST_EDGES: &[u64] = &[
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    12,
    15,
    20,
    25,
    30,
    40,
    50,
    60,
    80,
    100,
    150,
    200,
    300,
    500,
    1000,
    2000,
    5000,
    u64::MAX,
];

fn bin_of(v: u64) -> usize {
    HIST_EDGES.partition_point(|&e| e < v)
}

#[derive(serde::Serialize, Default)]
struct LocusStats {
    locus: String,
    alleles: usize,
    pairs: usize,
    complete: bool,
    len_min: u32,
    len_median: u32,
    len_max: u32,
    snps_mean: f64,
    snps_median: u32,
    snps_max: u32,
    /// fraction of pairs with at least one InDel event
    indel_pair_frac: f64,
    indel_events_mean: f64,
    indel_bases_mean: f64,
    /// totals over all pairs (coding stores only)
    #[serde(skip_serializing_if = "Option::is_none")]
    syn: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    nonsyn: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    frame_disrupted: Option<u64>,
}

struct Hists {
    snps: Vec<u64>,
    indel_events: Vec<u64>,
    indel_bases: Vec<u64>,
    /// pairs by (allele-level) distance split: snps only / with indels
    nonsyn_frac: Vec<u64>,
}

impl Hists {
    fn new() -> Self {
        let n = HIST_EDGES.len();
        Self {
            snps: vec![0; n],
            indel_events: vec![0; n],
            indel_bases: vec![0; n],
            nonsyn_frac: vec![0; 11],
        }
    }
    fn add(&mut self, o: &Hists) {
        for (a, b) in [
            (&mut self.snps, &o.snps),
            (&mut self.indel_events, &o.indel_events),
            (&mut self.indel_bases, &o.indel_bases),
            (&mut self.nonsyn_frac, &o.nonsyn_frac),
        ] {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
    }
}

fn locus_stats(locus: &str, d: &LocusData, complete: bool) -> (LocusStats, Hists) {
    let mut h = Hists::new();
    let mut lens: Vec<u32> = d.alleles.values().copied().filter(|&l| l > 0).collect();
    lens.sort_unstable();
    let mut snps: Vec<u32> = Vec::with_capacity(d.pairs.len());
    let (mut ind_pairs, mut ev, mut bases) = (0u64, 0u64, 0u64);
    let (mut syn, mut non, mut fd, mut coding) = (0u64, 0u64, 0u64, false);
    for s in d.pairs.values() {
        snps.push(s.snps);
        h.snps[bin_of(s.snps as u64)] += 1;
        h.indel_events[bin_of(s.indel_events as u64)] += 1;
        h.indel_bases[bin_of(s.indel_bases as u64)] += 1;
        if s.indel_events > 0 {
            ind_pairs += 1;
        }
        ev += s.indel_events as u64;
        bases += s.indel_bases as u64;
        if let Some((a, b, c)) = s.coding {
            coding = true;
            syn += a as u64;
            non += b as u64;
            fd += c as u64;
            if a + b > 0 {
                let f = (10 * b as u64 + (a + b) as u64 / 2) / (a + b) as u64;
                h.nonsyn_frac[f as usize] += 1;
            }
        }
    }
    snps.sort_unstable();
    let n = snps.len().max(1) as f64;
    let median = |v: &[u32]| if v.is_empty() { 0 } else { v[v.len() / 2] };
    let st = LocusStats {
        locus: locus.to_string(),
        alleles: d.alleles.len(),
        pairs: d.pairs.len(),
        complete,
        len_min: lens.first().copied().unwrap_or(0),
        len_median: median(&lens),
        len_max: lens.last().copied().unwrap_or(0),
        snps_mean: snps.iter().map(|&x| x as u64).sum::<u64>() as f64 / n,
        snps_median: median(&snps),
        snps_max: snps.last().copied().unwrap_or(0),
        indel_pair_frac: ind_pairs as f64 / n,
        indel_events_mean: ev as f64 / n,
        indel_bases_mean: bases as f64 / n,
        syn: coding.then_some(syn),
        nonsyn: coding.then_some(non),
        frame_disrupted: coding.then_some(fd),
    };
    (st, h)
}

fn stats(a: Stats) -> Result<(), String> {
    threads(a.threads);
    let src = Source::parse(&a.store);
    let m = src.manifest()?;
    let entries: Vec<(&String, &cgdist::store::LocusEntry)> = m.loci.iter().collect();
    let per: Vec<(LocusStats, Hists)> = entries
        .par_iter()
        .map(|(l, e)| {
            let d = src
                .locus_bytes(e)
                .and_then(|b| LocusData::decode(&b))
                .map_err(|err| format!("locus {l}: {err}"))?;
            Ok(locus_stats(l, &d, e.complete))
        })
        .collect::<Result<_, String>>()?;
    let mut total = Hists::new();
    for (_, h) in &per {
        total.add(h);
    }
    let edges: Vec<serde_json::Value> = HIST_EDGES
        .iter()
        .map(|&e| {
            if e == u64::MAX {
                serde_json::Value::Null
            } else {
                e.into()
            }
        })
        .collect();
    let coding = per.iter().any(|(s, _)| s.syn.is_some());
    let out = serde_json::json!({
        "format": "cgdist-store-stats",
        "format_version": 1,
        "store": src.describe(),
        "kind": if m.is_protein() { "protein" } else { "dna" },
        "params": m.params()?.to_string(),
        "cgdist_version": env!("CARGO_PKG_VERSION"),
        "loci_count": per.len(),
        "histograms": {
            "edges": edges,
            "snps": total.snps,
            "indel_events": total.indel_events,
            "indel_bases": total.indel_bases,
            "nonsyn_frac_deciles": if coding { serde_json::json!(total.nonsyn_frac) } else { serde_json::Value::Null },
        },
        "loci": per.iter().map(|(s, _)| s).collect::<Vec<_>>(),
    });
    let json = serde_json::to_vec(&out).map_err(|e| e.to_string())?;
    std::fs::write(&a.out, json).map_err(|e| format!("cannot write {}: {e}", a.out))?;
    println!("✅ stats of {} loci written to {}", per.len(), a.out);
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
    let mut store = Store::open_or_create(
        Path::new(&a.out),
        &cache.metadata.hasher_type,
        StoreParams::Dna(params),
    )?;
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
    println!("content      {}", m.params()?);
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
    println!(
        "{} {alleles}",
        if m.is_protein() {
            "proteins    "
        } else {
            "alleles     "
        }
    );
    println!("pairs        {pairs}");
    println!("size         {:.1} MB", bytes as f64 / 1e6);
    println!("modified     {}", m.last_modified);
    Ok(())
}

/// Sequences of one locus of a schema keyed by DNA CRC32 (None: no file).
fn schema_locus(dir: &str, locus: &str) -> Result<Option<BTreeMap<u32, Vec<u8>>>, String> {
    let path = ["fasta", "fa"]
        .iter()
        .map(|x| Path::new(dir).join(format!("{locus}.{x}")))
        .find(|p| p.exists());
    let Some(path) = path else { return Ok(None) };
    let reader =
        bio::io::fasta::Reader::from_file(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out = BTreeMap::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| format!("{}: {e}", path.display()))?;
        let mut h = crc32fast::Hasher::new();
        h.update(rec.seq());
        out.insert(h.finalize(), rec.seq().to_vec());
    }
    Ok(Some(out))
}

#[derive(Default)]
struct Realign {
    checked: u64,
    mismatches: u64,
    /// pairs whose alleles (or proteins) are not in the schema
    skipped: u64,
}

/// Recompute the selected pairs of one locus from the schema and compare.
/// Alleles (or proteins) of a store whose sequence digest differs from the
/// schema's sequence with the same hash; alleles without digest are skipped.
fn digest_mismatches(
    m: &cgdist::store::Manifest,
    d: &LocusData,
    seqs: &BTreeMap<u32, Vec<u8>>,
) -> Result<(usize, usize), String> {
    let (mut checked, mut bad) = (0usize, 0usize);
    let own: BTreeMap<u32, u64> = match m.params()? {
        StoreParams::Dna(_) => seqs
            .iter()
            .map(|(c, s)| (*c, cgdist::store::seq_digest(s)))
            .collect(),
        StoreParams::Protein(p) => {
            let code = p.code()?;
            seqs.values()
                .map(|s| {
                    let pr = protein_distance::protein_of(&code, s);
                    (
                        protein_distance::protein_hash(&pr),
                        cgdist::store::seq_digest(&pr),
                    )
                })
                .collect()
        }
    };
    for (h, dg) in &d.digests {
        if let Some(o) = own.get(h) {
            checked += 1;
            if o != dg {
                bad += 1;
            }
        }
    }
    Ok((checked, bad))
}

fn realign_locus(
    m: &cgdist::store::Manifest,
    locus: &str,
    d: &LocusData,
    seqs: &BTreeMap<u32, Vec<u8>>,
    fraction: f64,
) -> Result<Realign, String> {
    use cgdist::core::distance::verify_selected;
    let mut r = Realign::default();
    let sel: Vec<(&(u32, u32), &PairStats)> = d
        .pairs
        .iter()
        .filter(|((a, b), _)| verify_selected(*a, *b, fraction))
        .collect();
    if sel.is_empty() {
        return Ok(r);
    }
    match m.params()? {
        StoreParams::Protein(settings) => {
            let code = settings.code()?;
            let mut prot: BTreeMap<u32, Vec<u8>> = BTreeMap::new();
            for s in seqs.values() {
                let p = protein_distance::protein_of(&code, s);
                prot.insert(protein_distance::protein_hash(&p), p);
            }
            let res: Vec<Option<bool>> = sel
                .par_iter()
                .map(|((x, y), st)| {
                    let (Some(p1), Some(p2)) = (prot.get(x), prot.get(y)) else {
                        return None;
                    };
                    let fresh = protein_distance::align_proteins(&settings, p1, p2)
                        .map(|p| protein_distance::pair_stats(&p));
                    Some(fresh.as_ref() == Some(*st))
                })
                .collect();
            for (k, v) in sel.iter().zip(res) {
                match v {
                    None => r.skipped += 1,
                    Some(true) => r.checked += 1,
                    Some(false) => {
                        r.checked += 1;
                        r.mismatches += 1;
                        eprintln!("✗ {locus}: protein pair {}-{} differs", k.0 .0, k.0 .1);
                    }
                }
            }
        }
        StoreParams::Dna(p) => {
            let config = AlignmentConfig::custom(
                p.match_score,
                p.mismatch_penalty,
                p.gap_open,
                p.gap_extend,
            );
            let mut db = SequenceDatabase::new();
            for (crc, seq) in seqs {
                db.add_sequence(
                    locus.to_string(),
                    *crc,
                    SequenceInfo {
                        sequence: seq.clone(),
                        id: crc.to_string(),
                    },
                );
            }
            let mut set: HashSet<(String, u32, u32)> = HashSet::new();
            for ((x, y), _) in &sel {
                if seqs.contains_key(x) && seqs.contains_key(y) {
                    set.insert((locus.to_string(), *x, *y));
                } else {
                    r.skipped += 1;
                }
            }
            let mut engine = DistanceEngine::with_sequences(config, db, "crc32".to_string());
            engine.set_quiet(true);
            // every recomputed pair is also checked against parasail's
            // original kernel inside the engine (a difference aborts)
            engine.set_verify_fraction(1.0);
            if let Some(g) = m.genetic_code {
                engine.set_coding(GeneticCode::new(g.table, g.first_codon_as_met)?);
            }
            engine.precompute_alignments(&set, DistanceMode::SnpsOnly);
            let fresh = engine
                .drain_to_locus_data()
                .remove(locus)
                .unwrap_or_default();
            for ((x, y), st) in &sel {
                if !set.contains(&(locus.to_string(), *x, *y)) {
                    continue;
                }
                r.checked += 1;
                let mut f = fresh.pairs.get(&(*x, *y)).copied();
                if st.coding.is_none() {
                    if let Some(f) = f.as_mut() {
                        f.coding = None;
                    }
                }
                if f.as_ref() != Some(*st) {
                    r.mismatches += 1;
                    eprintln!("✗ {locus}: pair {x}-{y} differs: store {st:?}, recomputed {f:?}");
                }
            }
        }
    }
    Ok(r)
}

fn verify(a: Verify) -> Result<(), String> {
    if !(0.0..=1.0).contains(&a.realign) {
        return Err("--realign is a fraction between 0 and 1".into());
    }
    if a.realign > 0.0 && a.schema.is_none() {
        return Err("--realign needs --schema (the sequences to recompute from)".into());
    }
    threads(a.threads);
    let src = Source::parse(&a.store);
    let m = src.manifest()?;
    let mut problems = 0usize;
    let mut total = Realign::default();
    let (mut digests_checked, mut digest_bad) = (0usize, 0usize);
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
                if let Some(dir) = &a.schema {
                    match schema_locus(dir, locus)? {
                        Some(seqs) => {
                            let (c, bad) = digest_mismatches(&m, &d, &seqs)?;
                            digests_checked += c;
                            if bad > 0 {
                                eprintln!(
                                    "✗ {locus}: {bad} store alleles have the same hash as a different \
                                     sequence of the schema"
                                );
                                digest_bad += bad;
                                problems += 1;
                            }
                            if a.realign <= 0.0 {
                                continue;
                            }
                            let r = realign_locus(&m, locus, &d, &seqs, a.realign)?;
                            if r.mismatches > 0 {
                                problems += 1;
                            }
                            total.checked += r.checked;
                            total.mismatches += r.mismatches;
                            total.skipped += r.skipped;
                        }
                        None => eprintln!("⚠️  {locus}: not in the schema, not recomputed"),
                    }
                }
            }
            Err(err) => {
                eprintln!("✗ {locus}: {err}");
                problems += 1;
            }
        }
    }
    if a.schema.is_some() {
        println!(
            "🧬 sequence digests: {digests_checked} store alleles checked against the schema, \
             {digest_bad} mismatches"
        );
    }
    if a.realign > 0.0 {
        println!(
            "🔬 recomputed {} pairs from the schema: {} differences{}",
            total.checked,
            total.mismatches,
            if total.skipped > 0 {
                format!(
                    " ({} pairs skipped: sequence not in the schema)",
                    total.skipped
                )
            } else {
                String::new()
            }
        );
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
                if st.merged > 0 {
                    println!(
                        "🔀 {} loci kept their local pairs: merged with the downloaded ones",
                        st.merged
                    );
                }
                Ok(())
            })
        }
        Cmd::Info(a) => info(a),
        Cmd::Verify(a) => verify(a),
        Cmd::Stats(a) => stats(a),
    };
    if let Err(e) = res {
        eprintln!("❌ ERROR: {e}");
        std::process::exit(1);
    }
}
