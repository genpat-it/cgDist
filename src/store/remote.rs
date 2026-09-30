// remote.rs - Publish and mirror cache stores
//
// A store can be published in two ways, both static files:
//
// * as a directory (manifest.json + loci/<locus>.cgds), e.g. on a web server
//   or object store;
// * as a single *pack* file, for archives with a file-count limit (Zenodo):
//
//     b"CGDPACK1" | manifest length (u64 LE) | manifest JSON | locus blobs
//
//   The pack manifest is the store manifest with each locus' byte offset
//   (relative to the first blob). Over HTTP only the header, the manifest
//   and the needed loci are downloaded (Range requests); a server without
//   Range support is read once in full.
//
// Every locus downloaded is verified against its manifest sha256. Loci
// already present with the same checksum are skipped, so a repeated pull is
// an incremental update.

use super::{sha256_hex, LocusData, LocusEntry, Manifest, Store, MANIFEST_FILE};
use rayon::prelude::*;
use std::collections::HashSet;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

const PACK_MAGIC: &[u8; 8] = b"CGDPACK1";
const MAX_DOWNLOAD: u64 = 64 << 30;

/// Where a store is published.
pub enum Source {
    Dir(PathBuf),
    PackFile(PathBuf),
    HttpDir(String),
    HttpPack {
        url: String,
        /// whole pack, when the server ignored the Range header
        full: OnceLock<Vec<u8>>,
        /// byte offset of the first locus blob
        data_start: OnceLock<u64>,
    },
}

impl Source {
    pub fn parse(s: &str) -> Self {
        let is_pack = s.ends_with(".cgpack");
        if s.starts_with("http://") || s.starts_with("https://") {
            if is_pack {
                Source::HttpPack {
                    url: s.to_string(),
                    full: OnceLock::new(),
                    data_start: OnceLock::new(),
                }
            } else {
                Source::HttpDir(s.trim_end_matches('/').to_string())
            }
        } else {
            let p = PathBuf::from(s.strip_prefix("file://").unwrap_or(s));
            if is_pack || p.is_file() {
                Source::PackFile(p)
            } else {
                Source::Dir(p)
            }
        }
    }

    pub fn describe(&self) -> String {
        match self {
            Source::Dir(p) => format!("directory {}", p.display()),
            Source::PackFile(p) => format!("pack {}", p.display()),
            Source::HttpDir(u) => format!("URL {u}"),
            Source::HttpPack { url, .. } => format!("pack URL {url}"),
        }
    }

    /// Store manifest (with locus offsets for packs).
    pub fn manifest(&self) -> Result<Manifest, String> {
        match self {
            Source::Dir(root) => {
                let p = root.join(MANIFEST_FILE);
                Manifest::parse(
                    &std::fs::read(&p).map_err(|e| format!("cannot read {}: {e}", p.display()))?,
                )
            }
            Source::HttpDir(base) => {
                Manifest::parse(&http_get(&format!("{base}/{MANIFEST_FILE}"), None)?.0)
            }
            Source::PackFile(path) => {
                let mut f = std::fs::File::open(path)
                    .map_err(|e| format!("cannot open {}: {e}", path.display()))?;
                let mut head = [0u8; 16];
                f.read_exact(&mut head)
                    .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
                let len = pack_manifest_len(&head)?;
                let mut m = vec![0u8; len as usize];
                f.read_exact(&mut m)
                    .map_err(|e| format!("truncated pack {}: {e}", path.display()))?;
                Manifest::parse(&m)
            }
            Source::HttpPack {
                url,
                full,
                data_start,
            } => {
                let (head, partial) = http_get(url, Some((0, 15)))?;
                if !partial {
                    let _ = full.set(head);
                    let all = full.get().unwrap();
                    let len = pack_manifest_len(all.get(..16).ok_or("pack too short")?)?;
                    let _ = data_start.set(16 + len);
                    return Manifest::parse(
                        all.get(16..16 + len as usize).ok_or("truncated pack")?,
                    );
                }
                let len = pack_manifest_len(&head)?;
                let (m, _) = http_get(url, Some((16, 16 + len - 1)))?;
                let _ = data_start.set(16 + len);
                Manifest::parse(&m)
            }
        }
    }

    /// One locus file, verified against its manifest entry.
    pub fn locus_bytes(&self, entry: &LocusEntry) -> Result<Vec<u8>, String> {
        let bytes = match self {
            Source::Dir(root) => {
                let p = root.join(&entry.file);
                std::fs::read(&p).map_err(|e| format!("cannot read {}: {e}", p.display()))?
            }
            Source::HttpDir(base) => http_get(&format!("{base}/{}", entry.file), None)?.0,
            Source::PackFile(path) => {
                let off = entry.offset.ok_or("pack manifest entry without offset")?;
                let mut f = std::fs::File::open(path)
                    .map_err(|e| format!("cannot open {}: {e}", path.display()))?;
                let mut head = [0u8; 16];
                f.read_exact(&mut head).map_err(|e| e.to_string())?;
                let start = 16 + pack_manifest_len(&head)? + off;
                f.seek(SeekFrom::Start(start)).map_err(|e| e.to_string())?;
                let mut b = vec![0u8; entry.bytes as usize];
                f.read_exact(&mut b)
                    .map_err(|e| format!("truncated pack {}: {e}", path.display()))?;
                b
            }
            Source::HttpPack {
                url,
                full,
                data_start,
            } => {
                let off = entry.offset.ok_or("pack manifest entry without offset")?;
                let start = *data_start.get().ok_or("pack manifest not loaded")? + off;
                if let Some(all) = full.get() {
                    all.get(start as usize..(start + entry.bytes) as usize)
                        .ok_or("truncated pack")?
                        .to_vec()
                } else {
                    let (b, partial) = http_get(url, Some((start, start + entry.bytes - 1)))?;
                    if !partial {
                        return Err(format!("{url}: server stopped honouring Range requests"));
                    }
                    b
                }
            }
        };
        if sha256_hex(&bytes) != entry.sha256 {
            return Err(format!(
                "checksum mismatch for {} from {}",
                entry.file,
                self.describe()
            ));
        }
        Ok(bytes)
    }
}

fn pack_manifest_len(head: &[u8]) -> Result<u64, String> {
    if head.len() < 16 || &head[..8] != PACK_MAGIC {
        return Err("not a cgdist cache pack (bad magic)".to_string());
    }
    let len = u64::from_le_bytes(head[8..16].try_into().unwrap());
    if len > 1 << 32 {
        return Err("cache pack manifest too large".to_string());
    }
    Ok(len)
}

/// GET a URL (optionally a byte range, inclusive). Returns the body and
/// whether the server answered with a partial (206) response.
fn http_get(url: &str, range: Option<(u64, u64)>) -> Result<(Vec<u8>, bool), String> {
    let mut last_err = String::new();
    for attempt in 0..3 {
        if attempt > 0 {
            std::thread::sleep(std::time::Duration::from_millis(500 << attempt));
        }
        let mut req =
            ureq::get(url).set("User-Agent", concat!("cgdist/", env!("CARGO_PKG_VERSION")));
        if let Some((a, b)) = range {
            req = req.set("Range", &format!("bytes={a}-{b}"));
        }
        match req.call() {
            Ok(resp) => {
                let partial = resp.status() == 206;
                let mut buf = Vec::new();
                resp.into_reader()
                    .take(MAX_DOWNLOAD)
                    .read_to_end(&mut buf)
                    .map_err(|e| format!("download of {url} failed: {e}"))?;
                return Ok((buf, partial));
            }
            // a 4xx will not get better by retrying
            Err(ureq::Error::Status(code, _)) if (400..500).contains(&code) => {
                return Err(format!("download of {url} failed: HTTP {code}"));
            }
            Err(e) => last_err = format!("download of {url} failed: {e}"),
        }
    }
    Err(last_err)
}

/// Write a store as a single pack file.
pub fn pack(store_dir: &Path, out: &Path) -> Result<usize, String> {
    let store = Store::open(store_dir)?;
    let mut manifest = store.manifest.clone();
    let mut blobs: Vec<Vec<u8>> = Vec::with_capacity(manifest.loci.len());
    let mut offset = 0u64;
    for (locus, entry) in manifest.loci.iter_mut() {
        let b = store
            .read_locus_bytes(locus)?
            .ok_or_else(|| format!("locus {locus} missing from store"))?;
        entry.offset = Some(offset);
        offset += b.len() as u64;
        blobs.push(b);
    }
    let json =
        serde_json::to_vec(&manifest).map_err(|e| format!("cannot serialise manifest: {e}"))?;
    let tmp = out.with_extension("cgpack.tmp");
    let mut f = std::io::BufWriter::new(
        std::fs::File::create(&tmp).map_err(|e| format!("cannot write {}: {e}", tmp.display()))?,
    );
    let w = |f: &mut std::io::BufWriter<std::fs::File>, b: &[u8]| {
        f.write_all(b)
            .map_err(|e| format!("cannot write {}: {e}", tmp.display()))
    };
    w(&mut f, PACK_MAGIC)?;
    w(&mut f, &(json.len() as u64).to_le_bytes())?;
    w(&mut f, &json)?;
    for b in &blobs {
        w(&mut f, b)?;
    }
    f.flush().map_err(|e| e.to_string())?;
    drop(f);
    std::fs::rename(&tmp, out).map_err(|e| format!("cannot write {}: {e}", out.display()))?;
    Ok(blobs.len())
}

#[derive(Debug, Default)]
pub struct PullStats {
    pub fetched: usize,
    pub up_to_date: usize,
    /// loci merged with local content (pairs the source does not have)
    pub merged: usize,
    pub not_in_source: usize,
    pub bytes: u64,
}

/// Mirror `loci` (all loci when `None`) of `source` into the store at `dest`.
/// The destination adopts the source's hasher, alignment parameters, genetic
/// code and schema description; an existing destination must match them.
///
/// A locus the destination already holds with other content (e.g. pairs of
/// new alleles added locally) is merged, never replaced: the result has every
/// pair of both. A pair present in both with different statistics is an
/// error (the two stores are not from the same computation).
pub fn pull(
    source: &Source,
    dest: &Path,
    loci: Option<&HashSet<String>>,
) -> Result<PullStats, String> {
    let remote = source.manifest()?;
    let mut store = Store::open_or_create(dest, &remote.hasher, remote.params()?)?;
    if store.manifest.genetic_code.is_some() && store.manifest.genetic_code != remote.genetic_code {
        return Err(format!(
            "{} holds coding counts for another genetic code than {}",
            dest.display(),
            source.describe()
        ));
    }
    store.manifest.genetic_code = remote.genetic_code;
    let _lock = store.lock()?;
    if store.manifest.schema == Default::default() {
        store.manifest.schema = remote.schema.clone();
    }

    let mut stats = PullStats::default();
    let mut todo: Vec<(&String, &LocusEntry)> = Vec::new();
    match loci {
        Some(wanted) => {
            for l in wanted {
                match remote.loci.get_key_value(l) {
                    Some(kv) => todo.push(kv),
                    None => stats.not_in_source += 1,
                }
            }
        }
        None => todo.extend(remote.loci.iter()),
    }
    todo.retain(|(l, e)| match store.manifest.loci.get(*l) {
        Some(local) if local.sha256 == e.sha256 => {
            stats.up_to_date += 1;
            false
        }
        _ => true,
    });

    // Download in parallel, write sequentially (manifest updates are
    // serial). Chunks bound memory; the manifest is saved after each chunk so
    // an interrupted pull keeps what it already fetched.
    for chunk in todo.chunks(64) {
        let got: Vec<(String, Vec<u8>, &LocusEntry)> = chunk
            .par_iter()
            .map(|(l, e)| source.locus_bytes(e).map(|b| ((*l).clone(), b, *e)))
            .collect::<Result<_, _>>()?;
        for (locus, bytes, e) in got {
            stats.bytes += bytes.len() as u64;
            stats.fetched += 1;
            if store.manifest.loci.contains_key(&locus) {
                let local = store.read_locus(&locus)?.unwrap_or_default();
                let remote_data =
                    LocusData::decode(&bytes).map_err(|err| format!("locus {locus}: {err}"))?;
                match merge_locus(&locus, remote_data, local)? {
                    None => {
                        store.put_locus_bytes(&locus, &bytes, e.alleles, e.pairs, e.complete)?
                    }
                    Some(merged) => {
                        store.write_locus(&locus, &merged)?;
                        stats.merged += 1;
                    }
                }
            } else {
                store.put_locus_bytes(&locus, &bytes, e.alleles, e.pairs, e.complete)?;
            }
        }
        store.save_manifest()?;
    }
    store.save_manifest()?;
    Ok(stats)
}

/// Union of a downloaded locus and the local one. None when the local locus
/// adds nothing (the downloaded bytes can be stored as they are).
fn merge_locus(
    locus: &str,
    remote: LocusData,
    local: LocusData,
) -> Result<Option<LocusData>, String> {
    let mut merged = remote.clone();
    let mut added = false;
    for (k, l) in &local.pairs {
        match merged.pairs.get_mut(k) {
            Some(r) => {
                if (r.snps, r.indel_events, r.indel_bases)
                    != (l.snps, l.indel_events, l.indel_bases)
                    || (r.coding.is_some() && l.coding.is_some() && r.coding != l.coding)
                {
                    return Err(format!(
                        "locus {locus}: pair {}-{} differs between the local store and the source \
                         ({l:?} vs {r:?}); refusing to merge stores from different computations",
                        k.0, k.1
                    ));
                }
                if r.coding.is_none() && l.coding.is_some() {
                    r.coding = l.coding;
                    added = true;
                }
            }
            None => {
                merged.pairs.insert(*k, *l);
                added = true;
            }
        }
    }
    for (crc, d) in &local.digests {
        match merged.digests.get(crc) {
            Some(r) if r != d => {
                return Err(format!(
                    "locus {locus}: allele {crc} has different sequences in the local store and \
                     the source (hash collision); refusing to merge"
                ))
            }
            Some(_) => {}
            None => {
                merged.set_allele_digest(*crc, *d);
                added = true;
            }
        }
    }
    for (crc, len) in &local.alleles {
        let known = merged.alleles.get(crc).copied();
        if known.is_none() || (known == Some(0) && *len > 0) {
            added = true;
        }
        merged.set_allele_len(*crc, *len);
    }
    Ok(added.then_some(merged))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::alignment::AlignmentConfig;
    use crate::store::{AlignmentParams, LocusData, PairStats};

    #[test]
    fn pack_and_pull_from_file() {
        let base = std::env::temp_dir().join(format!("cgdist_pack_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        let (src, pack_path, dst) = (base.join("src"), base.join("s.cgpack"), base.join("dst"));
        let params = AlignmentParams::from(&AlignmentConfig::default());
        let mut st = Store::create(&src, "crc32", crate::store::StoreParams::Dna(params)).unwrap();
        for (i, locus) in ["L1", "L2", "L3"].iter().enumerate() {
            let mut d = LocusData::default();
            d.set_allele_len(10, 900);
            d.set_allele_len(20 + i as u32, 903);
            d.insert_pair(
                10,
                20 + i as u32,
                PairStats {
                    snps: i as u32,
                    indel_events: 1,
                    indel_bases: 3,
                    coding: None,
                },
            );
            st.write_locus(locus, &d).unwrap();
        }
        st.save_manifest().unwrap();
        assert_eq!(pack(&src, &pack_path).unwrap(), 3);

        let source = Source::parse(pack_path.to_str().unwrap());
        let wanted: HashSet<String> = ["L1", "L3", "L9"].iter().map(|s| s.to_string()).collect();
        let stats = pull(&source, &dst, Some(&wanted)).unwrap();
        assert_eq!((stats.fetched, stats.not_in_source), (2, 1));
        let got = Store::open(&dst).unwrap();
        assert_eq!(
            got.read_locus("L3").unwrap(),
            Store::open(&src).unwrap().read_locus("L3").unwrap()
        );
        assert!(got.read_locus("L2").unwrap().is_none());
        // second pull: nothing to fetch
        let again = pull(&source, &dst, Some(&wanted)).unwrap();
        assert_eq!((again.fetched, again.up_to_date), (0, 2));

        // a corrupted pack is detected
        let mut bytes = std::fs::read(&pack_path).unwrap();
        let last = bytes.len() - 1;
        bytes[last] ^= 0xff;
        std::fs::write(&pack_path, bytes).unwrap();
        let m = Source::parse(pack_path.to_str().unwrap())
            .manifest()
            .unwrap();
        let bad = m
            .loci
            .values()
            .filter(|e| {
                Source::parse(pack_path.to_str().unwrap())
                    .locus_bytes(e)
                    .is_err()
            })
            .count();
        assert_eq!(bad, 1);
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn pull_merges_local_pairs_instead_of_replacing_them() {
        let st = |snps| PairStats {
            snps,
            ..Default::default()
        };
        let mut remote = LocusData::default();
        remote.set_allele_len(1, 900);
        remote.set_allele_len(2, 900);
        remote.insert_pair(1, 2, st(3));
        // local = remote + a new allele (3) with its pairs
        let mut local = remote.clone();
        local.set_allele_len(3, 903);
        local.insert_pair(1, 3, st(1));
        local.insert_pair(2, 3, st(4));
        let m = merge_locus("L", remote.clone(), local).unwrap().unwrap();
        assert_eq!(m.pairs.len(), 3);
        assert_eq!(m.alleles[&3], 903);
        assert!(m.is_complete());
        // nothing local to add: the downloaded bytes are kept as they are
        assert!(merge_locus("L", remote.clone(), remote.clone())
            .unwrap()
            .is_none());
        // same hash, different sequence digest: refuse
        let mut r2 = remote.clone();
        r2.set_allele_digest(1, crate::store::seq_digest(b"AAAA"));
        let mut l2 = remote.clone();
        l2.set_allele_digest(1, crate::store::seq_digest(b"CCCC"));
        assert!(merge_locus("L", r2, l2).is_err());
        // same pair, different values: refuse
        let mut bad = remote.clone();
        bad.insert_pair(1, 2, st(4));
        assert!(merge_locus("L", remote, bad).is_err());
    }
}
