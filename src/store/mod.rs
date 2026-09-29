// store/mod.rs - Distributable per-locus alignment cache ("cache store")
//
// A cache store is a directory that can be shared, downloaded piecewise and
// updated incrementally:
//
//   <store>/manifest.json         alignment params, hasher, per-locus index
//   <store>/loci/<locus>.cgds     one compact binary file per locus
//
// Each locus file holds the allele table (CRC32 + sequence length) and the
// alignment statistics for pairs of those alleles. Pairs are addressed by
// allele index, so hashes are stored once per allele instead of twice per
// pair, and values are written column-wise as varints before LZ4 compression.
// A store is only valid for one set of alignment parameters and one hasher;
// both are checked before any entry is used.

use crate::core::alignment::AlignmentConfig;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

pub mod remote;

pub const STORE_FORMAT: &str = "cgdist-store";
pub const STORE_FORMAT_VERSION: u32 = 1;
pub const MANIFEST_FILE: &str = "manifest.json";
const LOCI_DIR: &str = "loci";
const LOCUS_EXT: &str = "cgds";
const LOCUS_MAGIC: &[u8; 4] = b"CGDS";
const LOCUS_VERSION: u8 = 1;
const LOCK_FILE: &str = ".lock";

/// Numeric alignment parameters. Descriptions are deliberately excluded:
/// two configs are interchangeable iff these four numbers match.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignmentParams {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl From<&AlignmentConfig> for AlignmentParams {
    fn from(c: &AlignmentConfig) -> Self {
        Self {
            match_score: c.match_score,
            mismatch_penalty: c.mismatch_penalty,
            gap_open: c.gap_open,
            gap_extend: c.gap_extend,
        }
    }
}

impl std::fmt::Display for AlignmentParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "match={}, mismatch={}, gap_open={}, gap_extend={}",
            self.match_score, self.mismatch_penalty, self.gap_open, self.gap_extend
        )
    }
}

/// Alignment statistics for one allele pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PairStats {
    pub snps: u32,
    pub indel_events: u32,
    pub indel_bases: u32,
    /// (synonymous, nonsynonymous, frame-disrupted) SNPs, when computed
    /// with the store's genetic code (Manifest::genetic_code)
    pub coding: Option<(u32, u32, u32)>,
}

/// Contents of one locus file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LocusData {
    /// Allele CRC32 -> nucleotide length (0 = unknown).
    pub alleles: BTreeMap<u32, u32>,
    /// (crc_lo, crc_hi) with crc_lo < crc_hi -> statistics.
    pub pairs: BTreeMap<(u32, u32), PairStats>,
}

impl LocusData {
    /// Record an allele length; a known length never gets overwritten by 0.
    pub fn set_allele_len(&mut self, crc: u32, len: u32) {
        let slot = self.alleles.entry(crc).or_insert(0);
        if len > 0 {
            *slot = len;
        }
    }

    /// Insert a pair (order-insensitive). Identical CRCs are ignored: their
    /// distance is 0 by definition and is never stored.
    pub fn insert_pair(&mut self, a: u32, b: u32, stats: PairStats) {
        if a == b {
            return;
        }
        let key = (a.min(b), a.max(b));
        self.alleles.entry(key.0).or_insert(0);
        self.alleles.entry(key.1).or_insert(0);
        self.pairs.insert(key, stats);
    }

    /// Merge `other` into `self`. On conflicting pairs `other` wins.
    pub fn merge(&mut self, other: LocusData) {
        for (crc, len) in other.alleles {
            self.set_allele_len(crc, len);
        }
        self.pairs.extend(other.pairs);
    }

    /// True when every pair of alleles in the table is present.
    pub fn is_complete(&self) -> bool {
        let n = self.alleles.len() as u64;
        self.pairs.len() as u64 == n * n.saturating_sub(1) / 2
    }

    pub fn encode(&self) -> Vec<u8> {
        let crcs: Vec<u32> = self.alleles.keys().copied().collect();
        let index: HashMap<u32, u32> = crcs
            .iter()
            .enumerate()
            .map(|(i, &c)| (c, i as u32))
            .collect();

        let mut body = Vec::with_capacity(16 + crcs.len() * 6 + self.pairs.len() * 5);
        put_varint(&mut body, crcs.len() as u64);
        let mut prev = 0u32;
        for (k, &crc) in crcs.iter().enumerate() {
            put_varint(&mut body, if k == 0 { crc } else { crc - prev } as u64);
            prev = crc;
        }
        for len in self.alleles.values() {
            put_varint(&mut body, *len as u64);
        }

        // BTreeMap order on (crc_lo, crc_hi) equals order on (i, j) because
        // the allele table is sorted by CRC.
        put_varint(&mut body, self.pairs.len() as u64);
        let idx: Vec<(u32, u32)> = self
            .pairs
            .keys()
            .map(|(a, b)| (index[a], index[b]))
            .collect();
        let mut prev_i = 0u32;
        for &(i, _) in &idx {
            put_varint(&mut body, (i - prev_i) as u64);
            prev_i = i;
        }
        let mut last: Option<(u32, u32)> = None;
        for &(i, j) in &idx {
            let d = match last {
                Some((pi, pj)) if pi == i => j - pj - 1,
                _ => j - i - 1,
            };
            put_varint(&mut body, d as u64);
            last = Some((i, j));
        }
        for s in self.pairs.values() {
            put_varint(&mut body, s.snps as u64);
        }
        for s in self.pairs.values() {
            put_varint(&mut body, s.indel_events as u64);
        }
        for s in self.pairs.values() {
            put_varint(&mut body, s.indel_bases as u64);
        }
        // optional coding columns, value + 1 (0 = not computed for the pair)
        let has_coding = self.pairs.values().any(|s| s.coding.is_some());
        body.push(u8::from(has_coding));
        if has_coding {
            for k in 0..3 {
                for s in self.pairs.values() {
                    let v = s.coding.map_or(0, |c| [c.0, c.1, c.2][k] as u64 + 1);
                    put_varint(&mut body, v);
                }
            }
        }

        let mut out = Vec::with_capacity(body.len() / 2 + 8);
        out.extend_from_slice(LOCUS_MAGIC);
        out.push(LOCUS_VERSION);
        out.extend_from_slice(&lz4_flex::compress_prepend_size(&body));
        out
    }

    pub fn decode(bytes: &[u8]) -> Result<Self, String> {
        let mut data = LocusData::default();
        decode_with(
            bytes,
            |alleles| {
                data.alleles = alleles.iter().copied().collect();
            },
            |a, b, s| {
                data.pairs.insert((a, b), s);
            },
        )?;
        Ok(data)
    }
}

/// Stream-decode a locus file: `on_alleles` receives the (crc, len) table,
/// `on_pair` every stored pair (crc_lo, crc_hi, stats). Lets callers filter
/// huge loci without materialising them.
pub fn decode_with(
    bytes: &[u8],
    mut on_alleles: impl FnMut(&[(u32, u32)]),
    mut on_pair: impl FnMut(u32, u32, PairStats),
) -> Result<(), String> {
    if bytes.len() < 5 || &bytes[..4] != LOCUS_MAGIC {
        return Err("not a cgdist locus file (bad magic)".to_string());
    }
    if bytes[4] != LOCUS_VERSION {
        return Err(format!(
            "unsupported locus file version {} (this cgdist reads {})",
            bytes[4], LOCUS_VERSION
        ));
    }
    let body = lz4_flex::decompress_size_prepended(&bytes[5..])
        .map_err(|e| format!("corrupt locus file: {e}"))?;
    let mut r = Reader { buf: &body, pos: 0 };

    let n = r.count()?;
    let mut alleles = Vec::with_capacity(n);
    let mut crc = 0u64;
    for k in 0..n {
        let d = r.varint()?;
        crc = if k == 0 { d } else { crc + d };
        if crc > u32::MAX as u64 || (k > 0 && d == 0) {
            return Err("corrupt locus file: allele table not strictly increasing".into());
        }
        alleles.push((crc as u32, 0u32));
    }
    for a in alleles.iter_mut() {
        a.1 = r.u32()?;
    }
    on_alleles(&alleles);

    let m = r.count()?;
    let mut ij = Vec::with_capacity(m);
    let mut i = 0u64;
    for _ in 0..m {
        i += r.varint()?;
        ij.push((i, 0u64));
    }
    let mut last: Option<(u64, u64)> = None;
    for e in ij.iter_mut() {
        let d = r.varint()?;
        let j = match last {
            Some((pi, pj)) if pi == e.0 => pj + 1 + d,
            _ => e.0 + 1 + d,
        };
        if e.0 >= n as u64 || j >= n as u64 {
            return Err("corrupt locus file: pair index out of range".into());
        }
        e.1 = j;
        last = Some(*e);
    }
    let mut stats = vec![PairStats::default(); m];
    for s in stats.iter_mut() {
        s.snps = r.u32()?;
    }
    for s in stats.iter_mut() {
        s.indel_events = r.u32()?;
    }
    for s in stats.iter_mut() {
        s.indel_bases = r.u32()?;
    }
    let has_coding = *body.get(r.pos).ok_or("corrupt locus file: truncated")?;
    r.pos += 1;
    match has_coding {
        0 => {}
        1 => {
            let mut cols = [vec![0u32; m], vec![0u32; m], vec![0u32; m]];
            for col in cols.iter_mut() {
                for v in col.iter_mut() {
                    *v = r.u32()?;
                }
            }
            for (i, s) in stats.iter_mut().enumerate() {
                let (a, b, c) = (cols[0][i], cols[1][i], cols[2][i]);
                s.coding = match (a, b, c) {
                    (0, 0, 0) => None,
                    (a, b, c) if a > 0 && b > 0 && c > 0 => Some((a - 1, b - 1, c - 1)),
                    _ => return Err("corrupt locus file: partial coding counts".into()),
                };
            }
        }
        _ => return Err("corrupt locus file: bad coding flag".into()),
    }
    if r.pos != body.len() {
        return Err("corrupt locus file: trailing bytes".into());
    }
    for ((i, j), s) in ij.into_iter().zip(stats) {
        on_pair(alleles[i as usize].0, alleles[j as usize].0, s);
    }
    Ok(())
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn varint(&mut self) -> Result<u64, String> {
        let mut v = 0u64;
        for shift in (0..64).step_by(7) {
            let b = *self
                .buf
                .get(self.pos)
                .ok_or("corrupt locus file: truncated")?;
            self.pos += 1;
            v |= ((b & 0x7f) as u64) << shift;
            if b < 0x80 {
                return Ok(v);
            }
        }
        Err("corrupt locus file: varint too long".into())
    }

    fn u32(&mut self) -> Result<u32, String> {
        let v = self.varint()?;
        u32::try_from(v).map_err(|_| "corrupt locus file: value overflow".to_string())
    }

    /// A count bounded by the remaining input (each item takes >= 1 byte),
    /// so a corrupt header cannot trigger a huge allocation.
    fn count(&mut self) -> Result<usize, String> {
        let v = self.varint()?;
        if v > (self.buf.len() - self.pos) as u64 {
            return Err("corrupt locus file: count exceeds data".into());
        }
        Ok(v as usize)
    }
}

/// Free-form description of the schema a store was built from.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct SchemaInfo {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LocusEntry {
    pub file: String,
    pub alleles: usize,
    pub pairs: usize,
    pub complete: bool,
    pub bytes: u64,
    pub sha256: String,
    /// byte offset of the locus blob in a pack file (packs only)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub offset: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    pub format: String,
    pub format_version: u32,
    pub hasher: String,
    pub alignment: AlignmentParams,
    #[serde(default)]
    pub schema: SchemaInfo,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    /// genetic code of the coding counts stored in the locus files, if any
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genetic_code: Option<crate::core::distance::GeneticCodeMeta>,
    pub cgdist_version: String,
    pub created: String,
    pub last_modified: String,
    pub loci: BTreeMap<String, LocusEntry>,
}

impl Manifest {
    pub fn new(hasher: &str, alignment: AlignmentParams) -> Self {
        let now = now();
        Self {
            format: STORE_FORMAT.to_string(),
            format_version: STORE_FORMAT_VERSION,
            hasher: hasher.to_string(),
            alignment,
            schema: SchemaInfo::default(),
            note: None,
            genetic_code: None,
            cgdist_version: env!("CARGO_PKG_VERSION").to_string(),
            created: now.clone(),
            last_modified: now,
            loci: BTreeMap::new(),
        }
    }

    pub fn parse(bytes: &[u8]) -> Result<Self, String> {
        let m: Manifest =
            serde_json::from_slice(bytes).map_err(|e| format!("invalid store manifest: {e}"))?;
        if m.format != STORE_FORMAT {
            return Err(format!("not a cgdist cache store (format '{}')", m.format));
        }
        if m.format_version > STORE_FORMAT_VERSION {
            return Err(format!(
                "store format v{} is newer than this cgdist supports (v{}); upgrade cgdist",
                m.format_version, STORE_FORMAT_VERSION
            ));
        }
        for (locus, e) in &m.loci {
            if e.file != locus_file_name(locus)? {
                return Err(format!(
                    "manifest entry for '{locus}' has unexpected file '{}'",
                    e.file
                ));
            }
        }
        Ok(m)
    }

    /// Refuse to mix alignment results computed under different settings.
    pub fn check_compatible(&self, hasher: &str, params: &AlignmentParams) -> Result<(), String> {
        if self.hasher != hasher {
            return Err(format!(
                "cache store uses hasher '{}', this run uses '{hasher}'",
                self.hasher
            ));
        }
        if &self.alignment != params {
            return Err(format!(
                "cache store alignment parameters differ:\n  store: {}\n  run:   {params}\n  \
                 results computed with different parameters cannot be mixed",
                self.alignment
            ));
        }
        Ok(())
    }
}

/// Relative path of a locus file. Locus names come from schema file names and
/// remote manifests, so anything that could escape the store is rejected.
pub fn locus_file_name(locus: &str) -> Result<String, String> {
    let ok = !locus.is_empty()
        && !locus.starts_with('.')
        && locus
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.' | '+'));
    if !ok {
        return Err(format!(
            "locus name '{locus}' cannot be used in a cache store"
        ));
    }
    Ok(format!("{LOCI_DIR}/{locus}.{LOCUS_EXT}"))
}

pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn now() -> String {
    chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string()
}

/// Write via a temp file + rename so readers never see a partial file.
fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let tmp = path.with_extension(format!("tmp{}", std::process::id()));
    let mut f =
        fs::File::create(&tmp).map_err(|e| format!("cannot write {}: {e}", tmp.display()))?;
    f.write_all(bytes)
        .and_then(|_| f.sync_all())
        .map_err(|e| format!("cannot write {}: {e}", tmp.display()))?;
    fs::rename(&tmp, path).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

/// A cache store on local disk.
pub struct Store {
    root: PathBuf,
    pub manifest: Manifest,
}

impl Store {
    pub fn open(root: &Path) -> Result<Self, String> {
        let path = root.join(MANIFEST_FILE);
        let bytes = fs::read(&path)
            .map_err(|e| format!("cannot read cache store manifest {}: {e}", path.display()))?;
        Ok(Self {
            root: root.to_path_buf(),
            manifest: Manifest::parse(&bytes)?,
        })
    }

    pub fn create(root: &Path, hasher: &str, params: AlignmentParams) -> Result<Self, String> {
        fs::create_dir_all(root.join(LOCI_DIR))
            .map_err(|e| format!("cannot create cache store {}: {e}", root.display()))?;
        let mut store = Self {
            root: root.to_path_buf(),
            manifest: Manifest::new(hasher, params),
        };
        store.save_manifest()?;
        Ok(store)
    }

    /// Open an existing store (checking compatibility) or create an empty one.
    pub fn open_or_create(
        root: &Path,
        hasher: &str,
        params: AlignmentParams,
    ) -> Result<Self, String> {
        if root.join(MANIFEST_FILE).exists() {
            let store = Self::open(root)?;
            store.manifest.check_compatible(hasher, &params)?;
            Ok(store)
        } else {
            Self::create(root, hasher, params)
        }
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    /// Raw bytes of a locus file, verified against the manifest checksum.
    pub fn read_locus_bytes(&self, locus: &str) -> Result<Option<Vec<u8>>, String> {
        let Some(entry) = self.manifest.loci.get(locus) else {
            return Ok(None);
        };
        let path = self.root.join(&entry.file);
        let bytes = fs::read(&path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        if sha256_hex(&bytes) != entry.sha256 {
            return Err(format!(
                "checksum mismatch for {} (file corrupt or modified)",
                path.display()
            ));
        }
        Ok(Some(bytes))
    }

    pub fn read_locus(&self, locus: &str) -> Result<Option<LocusData>, String> {
        match self.read_locus_bytes(locus)? {
            Some(b) => LocusData::decode(&b)
                .map(Some)
                .map_err(|e| format!("locus {locus}: {e}")),
            None => Ok(None),
        }
    }

    /// Write a locus file and update its manifest entry (call
    /// `save_manifest` afterwards to persist the index).
    pub fn write_locus(&mut self, locus: &str, data: &LocusData) -> Result<(), String> {
        let bytes = data.encode();
        self.put_locus_bytes(
            locus,
            &bytes,
            data.alleles.len(),
            data.pairs.len(),
            data.is_complete(),
        )
    }

    /// Store pre-encoded locus bytes (used when copying from a remote store).
    pub fn put_locus_bytes(
        &mut self,
        locus: &str,
        bytes: &[u8],
        alleles: usize,
        pairs: usize,
        complete: bool,
    ) -> Result<(), String> {
        let file = locus_file_name(locus)?;
        fs::create_dir_all(self.root.join(LOCI_DIR))
            .map_err(|e| format!("cannot create {}: {e}", self.root.display()))?;
        write_atomic(&self.root.join(&file), bytes)?;
        self.manifest.loci.insert(
            locus.to_string(),
            LocusEntry {
                file,
                alleles,
                pairs,
                complete,
                bytes: bytes.len() as u64,
                sha256: sha256_hex(bytes),
                offset: None,
            },
        );
        Ok(())
    }

    pub fn save_manifest(&mut self) -> Result<(), String> {
        self.manifest.last_modified = now();
        self.manifest.cgdist_version = env!("CARGO_PKG_VERSION").to_string();
        let json = serde_json::to_vec_pretty(&self.manifest)
            .map_err(|e| format!("cannot serialise manifest: {e}"))?;
        write_atomic(&self.root.join(MANIFEST_FILE), &json)
    }

    /// Take an exclusive write lock (a lock file) for the store.
    pub fn lock(&self) -> Result<StoreLock, String> {
        let path = self.root.join(LOCK_FILE);
        match fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
        {
            Ok(mut f) => {
                let _ = writeln!(f, "{}", std::process::id());
                Ok(StoreLock { path })
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => Err(format!(
                "cache store {} is locked by another cgdist process \
                 (remove {} if no other process is running)",
                self.root.display(),
                path.display()
            )),
            Err(e) => Err(format!("cannot lock {}: {e}", self.root.display())),
        }
    }

    /// Check every locus file against the manifest; returns problems found.
    pub fn verify(&self) -> Vec<String> {
        let mut problems = Vec::new();
        for (locus, entry) in &self.manifest.loci {
            match self.read_locus(locus) {
                Ok(Some(d)) => {
                    if d.alleles.len() != entry.alleles || d.pairs.len() != entry.pairs {
                        problems.push(format!("{locus}: counts differ from manifest"));
                    }
                    if d.is_complete() != entry.complete {
                        problems.push(format!("{locus}: 'complete' flag is wrong"));
                    }
                }
                Ok(None) => {}
                Err(e) => problems.push(e),
            }
        }
        problems
    }
}

/// Removes the lock file when dropped.
pub struct StoreLock {
    path: PathBuf,
}

impl Drop for StoreLock {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> LocusData {
        let mut d = LocusData::default();
        d.set_allele_len(10, 900);
        d.set_allele_len(4_000_000_000, 903);
        d.set_allele_len(77, 0);
        d.insert_pair(
            4_000_000_000,
            10,
            PairStats {
                snps: 3,
                indel_events: 1,
                indel_bases: 3,
                coding: None,
            },
        );
        d.insert_pair(
            10,
            77,
            PairStats {
                snps: 200,
                indel_events: 0,
                indel_bases: 0,
                coding: None,
            },
        );
        d.insert_pair(
            77,
            4_000_000_000,
            PairStats {
                snps: 0,
                indel_events: 2,
                indel_bases: 70_000,
                coding: None,
            },
        );
        d
    }

    #[test]
    fn roundtrip() {
        let d = sample();
        assert!(d.is_complete());
        let back = LocusData::decode(&d.encode()).unwrap();
        assert_eq!(back, d);
    }

    #[test]
    fn roundtrip_with_coding() {
        let mut d = sample();
        let k = *d.pairs.keys().next().unwrap();
        d.pairs.get_mut(&k).unwrap().coding = Some((2, 1, 0));
        let back = LocusData::decode(&d.encode()).unwrap();
        assert_eq!(back, d);
        assert_eq!(back.pairs[&k].coding, Some((2, 1, 0)));
        assert!(back.pairs.values().filter(|p| p.coding.is_none()).count() == 2);
    }

    #[test]
    fn roundtrip_empty_and_sparse() {
        let e = LocusData::default();
        assert_eq!(LocusData::decode(&e.encode()).unwrap(), e);

        let mut s = LocusData::default();
        for c in 0..50u32 {
            s.set_allele_len(c * 1000 + 7, 100 + c);
        }
        s.insert_pair(
            7,
            49_007,
            PairStats {
                snps: 1,
                ..Default::default()
            },
        );
        s.insert_pair(
            20_007,
            30_007,
            PairStats {
                snps: 2,
                ..Default::default()
            },
        );
        assert!(!s.is_complete());
        assert_eq!(LocusData::decode(&s.encode()).unwrap(), s);
    }

    #[test]
    fn identical_pair_ignored_and_merge_keeps_lengths() {
        let mut a = LocusData::default();
        a.insert_pair(5, 5, PairStats::default());
        assert!(a.pairs.is_empty());
        a.set_allele_len(5, 300);
        let mut b = LocusData::default();
        b.insert_pair(
            5,
            6,
            PairStats {
                snps: 1,
                ..Default::default()
            },
        );
        a.merge(b);
        assert_eq!(a.alleles[&5], 300);
        assert_eq!(a.alleles[&6], 0);
        assert_eq!(a.pairs.len(), 1);
    }

    #[test]
    fn corrupt_input_rejected() {
        let bytes = sample().encode();
        assert!(LocusData::decode(&bytes[..3]).is_err());
        let mut bad = bytes.clone();
        bad[0] = b'X';
        assert!(LocusData::decode(&bad).is_err());
        let mut v = bytes.clone();
        v[4] = 99;
        assert!(LocusData::decode(&v).is_err());
        assert!(LocusData::decode(&bytes[..bytes.len() - 2]).is_err());
    }

    #[test]
    fn locus_names_are_sanitised() {
        assert!(locus_file_name("cgMLST-00079395").is_ok());
        assert!(locus_file_name("Lm_0001.cds").is_ok());
        for bad in ["", "../x", "a/b", ".hidden", "a b", "a\\b"] {
            assert!(locus_file_name(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn store_write_read_verify_and_compat() {
        let dir = std::env::temp_dir().join(format!("cgdist_store_test_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let params = AlignmentParams::from(&AlignmentConfig::default());
        let mut st = Store::create(&dir, "crc32", params).unwrap();
        st.write_locus("L1", &sample()).unwrap();
        st.save_manifest().unwrap();

        let st2 = Store::open(&dir).unwrap();
        assert_eq!(st2.read_locus("L1").unwrap().unwrap(), sample());
        assert!(st2.read_locus("missing").unwrap().is_none());
        assert!(st2.verify().is_empty());
        assert!(st2.manifest.check_compatible("crc32", &params).is_ok());
        let other = AlignmentParams {
            gap_open: 8,
            ..params
        };
        assert!(st2.manifest.check_compatible("crc32", &other).is_err());
        assert!(st2.manifest.check_compatible("sha256", &params).is_err());

        // a preset with the same numbers but a different description is compatible
        let dna = AlignmentConfig::from_mode("dna").unwrap();
        assert!(st2
            .manifest
            .check_compatible("crc32", &AlignmentParams::from(&dna))
            .is_ok());

        // tampering is detected
        let f = dir.join(locus_file_name("L1").unwrap());
        let mut b = fs::read(&f).unwrap();
        let last = b.len() - 1;
        b[last] ^= 1;
        fs::write(&f, b).unwrap();
        assert!(st2.read_locus("L1").is_err());

        // lock is exclusive and released on drop
        let l = st2.lock().unwrap();
        assert!(st2.lock().is_err());
        drop(l);
        assert!(st2.lock().is_ok());
        let _ = fs::remove_dir_all(&dir);
    }
}
