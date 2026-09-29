// protein_distance.rs - Amino-acid level distances between alleles
//
// Alleles are translated on the fly from the DNA schema (configurable genetic
// code). Each allele gets a protein hash (CRC32 of the amino-acid sequence,
// terminal stop removed); alleles with the same protein are identical at
// protein level and need no alignment. Pairs of distinct proteins are
// aligned globally with parasail (configurable substitution matrix and gap
// penalties) and their amino-acid substitutions / InDel events / InDel
// residues are counted exactly as the DNA statistics are
// (compute_alignment_stats on the gapped protein strings).
//
// Results live in a separate protein cache (--protein-cache-file): its
// values depend on the genetic code, matrix and gap penalties, which are
// recorded in its metadata and must match to be reused.

use crate::core::alignment::compute_alignment_stats;
use crate::core::protein::GeneticCode;
use crate::data::SequenceDatabase;
use parasail_rs::{Aligner, Matrix};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap, HashSet};

/// Settings that protein-level results depend on.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProteinSettings {
    pub translation_table: u32,
    pub first_codon_as_met: bool,
    /// parasail matrix name (blosum62, pam250, ...) or a matrix file path
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl ProteinSettings {
    pub fn code(&self) -> Result<GeneticCode, String> {
        GeneticCode::new(self.translation_table, self.first_codon_as_met)
    }

    fn load_matrix(&self) -> Result<Matrix, String> {
        if std::path::Path::new(&self.matrix).is_file() {
            Matrix::from_file(&self.matrix)
                .map_err(|e| format!("cannot read matrix file {}: {e:?}", self.matrix))
        } else {
            Matrix::from(&self.matrix.to_lowercase()).map_err(|_| {
                format!(
                    "unknown substitution matrix '{}' (use blosum30..blosum100, pam10..pam500, or a matrix file)",
                    self.matrix
                )
            })
        }
    }

    /// Check that the settings are usable (matrix exists, penalties valid).
    pub fn validate(&self) -> Result<(), String> {
        self.code()?;
        self.load_matrix()?;
        if self.gap_open < 0 || self.gap_extend < 0 {
            return Err("--aa-gap-open and --aa-gap-extend must be non-negative".to_string());
        }
        Ok(())
    }
}

/// Amino-acid statistics of one pair of distinct proteins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProteinPair {
    pub aa_subs: u32,
    pub aa_indel_events: u32,
    pub aa_indel_residues: u32,
    pub aa_length1: u32,
    pub aa_length2: u32,
}

/// Protein of a CDS for comparison: translated, terminal stop removed.
pub fn protein_of(code: &GeneticCode, cds: &[u8]) -> Vec<u8> {
    let mut p = code.translate(cds);
    if p.last() == Some(&b'*') {
        p.pop();
    }
    p
}

fn protein_hash(p: &[u8]) -> u32 {
    let mut h = crc32fast::Hasher::new();
    h.update(p);
    h.finalize()
}

thread_local! {
    static AA_ALIGNERS: RefCell<Vec<(ProteinSettings, i32, Aligner)>> = const { RefCell::new(Vec::new()) };
}

/// Align two proteins (global, traceback) and count substitutions and
/// InDels. Scan kernel at 16-bit, escalated to 32/64-bit on saturation, as
/// for DNA. None if parasail cannot align them.
pub fn align_proteins(s: &ProteinSettings, p1: &[u8], p2: &[u8]) -> Option<ProteinPair> {
    if p1.is_empty() || p2.is_empty() {
        // parasail needs non-empty input: an empty protein vs a protein of
        // length L is one InDel event of L residues
        let l = p1.len().max(p2.len()) as u32;
        return Some(ProteinPair {
            aa_subs: 0,
            aa_indel_events: u32::from(l > 0),
            aa_indel_residues: l,
            aa_length1: p1.len() as u32,
            aa_length2: p2.len() as u32,
        });
    }
    let run = |width: i32| -> Option<(Option<ProteinPair>, bool)> {
        AA_ALIGNERS.with(|cell| {
            let mut v = cell.borrow_mut();
            let idx = match v.iter().position(|(k, w, _)| k == s && *w == width) {
                Some(i) => i,
                None => {
                    let m = s.load_matrix().ok()?;
                    let a = Aligner::new()
                        .matrix(m)
                        .gap_open(s.gap_open)
                        .gap_extend(s.gap_extend)
                        .global()
                        .use_trace()
                        .scan()
                        .solution_width(width)
                        .build();
                    v.push((s.clone(), width, a));
                    v.len() - 1
                }
            };
            let res = v[idx].2.align(Some(p1), p2).ok()?;
            if res.is_saturated() {
                return Some((None, true));
            }
            let tb = res.get_traceback_strings(p1, p2).ok()?;
            let (subs, ev, residues) = compute_alignment_stats(&tb.query, &tb.reference);
            Some((
                Some(ProteinPair {
                    aa_subs: subs as u32,
                    aa_indel_events: ev as u32,
                    aa_indel_residues: residues as u32,
                    aa_length1: p1.len() as u32,
                    aa_length2: p2.len() as u32,
                }),
                false,
            ))
        })
    };
    for width in [16, 32, 64] {
        match run(width)? {
            (Some(r), _) => return Some(r),
            (None, true) => continue,
            (None, false) => return None,
        }
    }
    panic!("protein alignment saturated even at 64-bit precision");
}

#[derive(Serialize, Deserialize)]
struct ProteinCacheFile {
    metadata: ProteinCacheMetadata,
    /// "locus:dna_crc" -> protein hash
    proteins: BTreeMap<String, u32>,
    /// "locus:protein_hash_lo:protein_hash_hi" -> statistics
    data: BTreeMap<String, ProteinPair>,
}

#[derive(Serialize, Deserialize)]
struct ProteinCacheMetadata {
    format: String,
    format_version: u32,
    cgdist_version: String,
    settings: ProteinSettings,
    created: String,
}

const PROTEIN_CACHE_FORMAT: &str = "cgdist-protein-cache";

/// Protein-level view of the alleles of a run, and cached pair results.
pub struct ProteinStore {
    settings: ProteinSettings,
    code: GeneticCode,
    /// (locus, dna crc) -> protein hash
    map: HashMap<(String, u32), u32>,
    /// (locus, lo, hi) -> statistics, lo < hi protein hashes
    pairs: HashMap<(String, u32, u32), ProteinPair>,
    has_new: bool,
    /// allele pairs (locus, dna lo, dna hi) whose protein is unknown (no
    /// sequence in the schema)
    pub unknown_alleles: usize,
}

impl ProteinStore {
    pub fn new(settings: ProteinSettings) -> Result<Self, String> {
        settings.validate()?;
        let code = settings.code()?;
        Ok(Self {
            settings,
            code,
            map: HashMap::new(),
            pairs: HashMap::new(),
            has_new: false,
            unknown_alleles: 0,
        })
    }

    pub fn settings(&self) -> &ProteinSettings {
        &self.settings
    }

    /// Load a protein cache; its settings must equal the run's.
    pub fn load(&mut self, path: &str) -> Result<usize, String> {
        let bytes =
            std::fs::read(path).map_err(|e| format!("cannot read protein cache {path}: {e}"))?;
        let raw = lz4_flex::decompress_size_prepended(&bytes)
            .map_err(|e| format!("cannot decompress protein cache {path}: {e}"))?;
        let f: ProteinCacheFile = serde_json::from_slice(&raw)
            .map_err(|e| format!("not a cgdist protein cache ({path}): {e}"))?;
        if f.metadata.format != PROTEIN_CACHE_FORMAT {
            return Err(format!("{path} is not a cgdist protein cache"));
        }
        if f.metadata.settings != self.settings {
            return Err(format!(
                "protein cache {path} was built with different settings:\n  cache: {:?}\n  run:   {:?}\n  \
                 use the same --translation-table/--aa-matrix/--aa-gap-open/--aa-gap-extend, or another file",
                f.metadata.settings, self.settings
            ));
        }
        for (k, v) in f.proteins {
            if let Some((locus, crc)) = k.rsplit_once(':') {
                if let Ok(crc) = crc.parse() {
                    self.map.insert((locus.to_string(), crc), v);
                }
            }
        }
        for (k, v) in f.data {
            let mut it = k.rsplitn(3, ':');
            if let (Some(hi), Some(lo), Some(locus)) = (it.next(), it.next(), it.next()) {
                if let (Ok(lo), Ok(hi)) = (lo.parse(), hi.parse()) {
                    self.pairs.insert((locus.to_string(), lo, hi), v);
                }
            }
        }
        Ok(self.pairs.len())
    }

    pub fn save(&self, path: &str) -> Result<(), String> {
        let f = ProteinCacheFile {
            metadata: ProteinCacheMetadata {
                format: PROTEIN_CACHE_FORMAT.to_string(),
                format_version: 1,
                cgdist_version: env!("CARGO_PKG_VERSION").to_string(),
                settings: self.settings.clone(),
                created: chrono::Utc::now()
                    .format("%Y-%m-%d %H:%M:%S UTC")
                    .to_string(),
            },
            proteins: self
                .map
                .iter()
                .map(|((l, c), p)| (format!("{l}:{c}"), *p))
                .collect(),
            data: self
                .pairs
                .iter()
                .map(|((l, a, b), v)| (format!("{l}:{a}:{b}"), *v))
                .collect(),
        };
        let json =
            serde_json::to_vec(&f).map_err(|e| format!("cannot serialise protein cache: {e}"))?;
        std::fs::write(path, lz4_flex::compress_prepend_size(&json))
            .map_err(|e| format!("cannot write protein cache {path}: {e}"))
    }

    pub fn has_new_entries(&self) -> bool {
        self.has_new
    }

    /// Translate the alleles of the run and align every missing pair of
    /// distinct proteins. `dna_pairs` are the (locus, crc lo, crc hi) allele
    /// pairs of the run.
    pub fn precompute(
        &mut self,
        seq_db: &SequenceDatabase,
        dna_pairs: &HashSet<(String, u32, u32)>,
    ) -> Result<(), String> {
        // 1. proteins of all alleles involved
        let mut alleles: HashSet<(&str, u32)> = HashSet::new();
        for (l, a, b) in dna_pairs {
            alleles.insert((l.as_str(), *a));
            alleles.insert((l.as_str(), *b));
        }
        let mut protein_seq: HashMap<(String, u32), Vec<u8>> = HashMap::new();
        for (locus, crc) in alleles {
            let Some(si) = seq_db.get_sequence(locus, crc) else {
                continue;
            };
            let p = protein_of(&self.code, &si.sequence);
            let h = protein_hash(&p);
            // distinct proteins of one locus must not share a hash
            match protein_seq.get(&(locus.to_string(), h)) {
                Some(other) if *other != p => {
                    return Err(format!(
                        "two different proteins of locus {locus} have the same CRC32 ({h}); \
                         protein-level modes cannot be used for this locus"
                    ))
                }
                _ => {
                    protein_seq.insert((locus.to_string(), h), p);
                }
            }
            self.map.insert((locus.to_string(), crc), h);
        }
        // 2. distinct protein pairs not in the cache
        let mut todo: HashSet<(String, u32, u32)> = HashSet::new();
        self.unknown_alleles = 0;
        for (l, a, b) in dna_pairs {
            match (
                self.map.get(&(l.clone(), *a)),
                self.map.get(&(l.clone(), *b)),
            ) {
                (Some(&pa), Some(&pb)) if pa != pb => {
                    let key = (l.clone(), pa.min(pb), pa.max(pb));
                    if !self.pairs.contains_key(&key) {
                        todo.insert(key);
                    }
                }
                (Some(_), Some(_)) => {}
                _ => self.unknown_alleles += 1,
            }
        }
        println!(
            "🧬 Protein level: {} allele pairs, {} distinct protein pairs to align (table {}, {}, gaps {}/{})",
            dna_pairs.len(),
            todo.len(),
            self.settings.translation_table,
            self.settings.matrix,
            self.settings.gap_open,
            self.settings.gap_extend
        );
        let todo: Vec<(String, u32, u32)> = todo.into_iter().collect();
        let settings = &self.settings;
        let results: Vec<((String, u32, u32), Option<ProteinPair>)> = todo
            .into_par_iter()
            .map(|k| {
                let p1 = &protein_seq[&(k.0.clone(), k.1)];
                let p2 = &protein_seq[&(k.0.clone(), k.2)];
                let r = align_proteins(settings, p1, p2);
                (k, r)
            })
            .collect();
        for (k, r) in results {
            if let Some(r) = r {
                self.pairs.insert(k, r);
                self.has_new = true;
            }
        }
        if self.unknown_alleles > 0 {
            eprintln!(
                "⚠️  WARNING: {} allele pairs have an allele without sequence in the schema; \
                 they count as 0 at protein level",
                self.unknown_alleles
            );
        }
        Ok(())
    }

    /// Protein-level view of a DNA allele pair: None if a protein is
    /// unknown; Some(None) if both alleles encode the same protein; else the
    /// pair statistics.
    pub fn lookup(&self, locus: &str, crc1: u32, crc2: u32) -> Option<Option<&ProteinPair>> {
        let pa = *self.map.get(&(locus.to_string(), crc1))?;
        let pb = *self.map.get(&(locus.to_string(), crc2))?;
        if pa == pb {
            return Some(None);
        }
        self.pairs
            .get(&(locus.to_string(), pa.min(pb), pa.max(pb)))
            .map(Some)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings() -> ProteinSettings {
        ProteinSettings {
            translation_table: 11,
            first_codon_as_met: true,
            matrix: "blosum62".into(),
            gap_open: 11,
            gap_extend: 1,
        }
    }

    #[test]
    fn protein_of_strips_terminal_stop() {
        let code = GeneticCode::default();
        assert_eq!(protein_of(&code, b"GTGGGATAA"), b"MG");
        assert_eq!(protein_of(&code, b"ATGGGA"), b"MG");
    }

    #[test]
    fn protein_alignment_counts() {
        let s = settings();
        let r = align_proteins(&s, b"MKTAYIAKQR", b"MKTAYIAKQR").unwrap();
        assert_eq!(
            (r.aa_subs, r.aa_indel_events, r.aa_indel_residues),
            (0, 0, 0)
        );
        let r = align_proteins(&s, b"MKTAYIAKQR", b"MKTAWIAKQR").unwrap();
        assert_eq!((r.aa_subs, r.aa_indel_events), (1, 0));
        let r = align_proteins(
            &s,
            b"MKTAYIAKQRQISFVKSHFSRQ",
            b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
        )
        .unwrap();
        assert_eq!(
            (r.aa_subs, r.aa_indel_events, r.aa_indel_residues),
            (0, 1, 11)
        );
        assert!(ProteinStore::new(ProteinSettings {
            matrix: "nope".into(),
            ..settings()
        })
        .is_err());
    }
}
