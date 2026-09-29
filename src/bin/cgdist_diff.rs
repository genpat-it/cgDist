// cgdist_diff.rs - Show where two alleles of a locus differ
//
// Given a schema and two allele hashes (CRC32, as in hashed profiles and in
// the cache keys), aligns the two alleles exactly as cgdist does and lists
// every SNP, insertion and deletion with its position and bases, plus the
// CIGAR. With --cache-file it also checks the counts against the cache entry.

use argh::FromArgs;
use bio::io::fasta;
use cgdist::core::alignment::{cigar_from_aligned, AlignmentConfig};
use cgdist::core::distance::{align_pair_with_strings, ModernCache};
use cgdist::core::protein::{self, GeneticCode, SnpAnnotation, SnpEffect};
use cgdist::core::protein_distance::{align_proteins_with_strings, protein_of, ProteinSettings};
use std::path::Path;

#[derive(FromArgs)]
/// cgdist-diff: list the SNPs and InDels between two alleles of one locus.
///
/// Alleles are ordered as in the cache: allele 1 is the one with the smaller
/// hash. Positions are 1-based on each allele.
struct Args {
    /// schema directory (one <locus>.fasta per locus)
    #[argh(option)]
    schema: String,

    /// locus name (FASTA file stem)
    #[argh(option)]
    locus: String,

    /// hash (CRC32) of the first allele
    #[argh(option)]
    hash1: u32,

    /// hash (CRC32) of the second allele
    #[argh(option)]
    hash2: u32,

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

    /// cgdist cache (.lz4) to compare the counts with
    #[argh(option)]
    cache_file: Option<String>,

    /// print the differences as TSV instead of text
    #[argh(switch)]
    tsv: bool,

    /// NCBI translation table for the protein annotation (default: 11,
    /// Bacterial, Archaeal and Plant Plastid)
    #[argh(option, default = "11")]
    translation_table: u32,

    /// do not read an alternative start codon (e.g. GTG, TTG) as Met when it
    /// is the first codon
    #[argh(switch)]
    no_first_codon_as_met: bool,

    /// substitution matrix for the protein alignment (default: blosum62)
    #[argh(option, default = "String::from(\"blosum62\")")]
    aa_matrix: String,

    /// protein gap open penalty (default: 11)
    #[argh(option, default = "11")]
    aa_gap_open: i32,

    /// protein gap extend penalty (default: 1)
    #[argh(option, default = "1")]
    aa_gap_extend: i32,

    /// also list the amino-acid substitutions and InDels of the protein
    /// alignment (TSV rows AA_SUB / AA_INS / AA_DEL)
    #[argh(switch)]
    protein_diffs: bool,
}

/// One difference between the two alleles.
struct Diff {
    kind: &'static str, // SNP, INS (bases only in allele 2), DEL (bases only in allele 1)
    event: usize,       // InDel event number (0 for SNPs), as counted by cgdist
    a1: (usize, usize), // allele 1 range, 1-based; for INS: (p, p) = after position p
    a2: (usize, usize), // allele 2 range, 1-based; for DEL: (p, p) = after position p
    bases1: String,
    bases2: String,
}

/// Walk the gapped strings. A maximal run of gap columns (of either kind)
/// is one InDel event, as in compute_alignment_stats.
fn differences(query: &[u8], reference: &[u8]) -> Vec<Diff> {
    let mut out: Vec<Diff> = Vec::new();
    let (mut p1, mut p2) = (0usize, 0usize);
    let mut event = 0usize;
    let mut in_gap = false;
    for (&a, &b) in query.iter().zip(reference) {
        if a == b'-' || b == b'-' {
            if !in_gap {
                event += 1;
                in_gap = true;
            }
            let kind = if a == b'-' { "INS" } else { "DEL" };
            if a == b'-' {
                p2 += 1;
            } else {
                p1 += 1;
            }
            match out.last_mut() {
                Some(d) if d.kind == kind && d.event == event => {
                    if kind == "INS" {
                        d.a2.1 = p2;
                        d.bases2.push(b as char);
                    } else {
                        d.a1.1 = p1;
                        d.bases1.push(a as char);
                    }
                }
                _ => out.push(if kind == "INS" {
                    Diff {
                        kind,
                        event,
                        a1: (p1, p1),
                        a2: (p2, p2),
                        bases1: String::new(),
                        bases2: (b as char).to_string(),
                    }
                } else {
                    Diff {
                        kind,
                        event,
                        a1: (p1, p1),
                        a2: (p2, p2),
                        bases1: (a as char).to_string(),
                        bases2: String::new(),
                    }
                }),
            }
        } else {
            in_gap = false;
            p1 += 1;
            p2 += 1;
            // byte-wise, as cgdist counts SNPs
            if a != b {
                out.push(Diff {
                    kind: "SNP",
                    event: 0,
                    a1: (p1, p1),
                    a2: (p2, p2),
                    bases1: (a as char).to_string(),
                    bases2: (b as char).to_string(),
                });
            }
        }
    }
    out
}

fn load_alleles(
    schema: &str,
    locus: &str,
    want: [u32; 2],
) -> Result<[(String, Vec<u8>); 2], String> {
    let dir = Path::new(schema);
    let path = ["fasta", "fa"]
        .iter()
        .map(|ext| dir.join(format!("{locus}.{ext}")))
        .find(|p| p.exists())
        .ok_or_else(|| format!("no {locus}.fasta (or .fa) in {schema}"))?;
    let reader = fasta::Reader::from_file(&path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut found: [Option<(String, Vec<u8>)>; 2] = [None, None];
    for rec in reader.records() {
        let rec = rec.map_err(|e| format!("invalid FASTA record in {}: {e}", path.display()))?;
        // CRC32 of the raw sequence bytes, as cgdist and chewBBACA compute it
        let mut h = crc32fast::Hasher::new();
        h.update(rec.seq());
        let crc = h.finalize();
        for (k, w) in want.iter().enumerate() {
            if crc == *w && found[k].is_none() {
                found[k] = Some((rec.id().to_string(), rec.seq().to_vec()));
            }
        }
    }
    let missing: Vec<String> = (0..2)
        .filter(|&k| found[k].is_none())
        .map(|k| want[k].to_string())
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "allele hash(es) {} not found in {} (profiles called with another schema?)",
            missing.join(", "),
            path.display()
        ));
    }
    let [a, b] = found;
    Ok([a.unwrap(), b.unwrap()])
}

fn cache_entry(
    path: &str,
    key: &str,
) -> Result<Option<(usize, usize, usize, AlignmentConfig)>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let raw = lz4_flex::decompress_size_prepended(&bytes)
        .map_err(|e| format!("cannot decompress {path}: {e}"))?;
    let cache: ModernCache =
        serde_json::from_slice(&raw).map_err(|e| format!("not a cgdist cache ({path}): {e}"))?;
    Ok(cache.data.get(key).map(|v| {
        (
            v.snps,
            v.indel_events,
            v.indel_bases,
            cache.metadata.alignment_config.clone(),
        )
    }))
}

fn run(args: Args) -> Result<i32, String> {
    let custom = [
        args.match_score,
        args.mismatch_penalty,
        args.gap_open,
        args.gap_extend,
    ];
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
        AlignmentConfig::from_mode(&args.alignment_mode)?
    };

    // cache order: allele 1 = smaller hash (the alignment query)
    let (lo, hi) = (args.hash1.min(args.hash2), args.hash1.max(args.hash2));
    if lo == hi {
        return Err("the two hashes are identical: the alleles do not differ".into());
    }
    let [(id1, s1), (id2, s2)] = load_alleles(&args.schema, &args.locus, [lo, hi])?;
    let aln =
        align_pair_with_strings(&config, &s1, &s2).ok_or("the alleles could not be aligned")?;
    let diffs = differences(aln.query.as_bytes(), aln.reference.as_bytes());
    let cigar = cigar_from_aligned(aln.query.as_bytes(), aln.reference.as_bytes());

    // coding effects; positions in `Diff` are 1-based
    let code = GeneticCode::new(args.translation_table, !args.no_first_codon_as_met)?;
    let map = protein::map_allele1_to_allele2(aln.query.as_bytes(), aln.reference.as_bytes());
    let snp_ann: Vec<Option<SnpAnnotation>> = diffs
        .iter()
        .map(|d| {
            (d.kind == "SNP")
                .then(|| protein::classify_snp(&code, &s1, &s2, &map, d.a1.0 - 1, d.a2.0 - 1))
        })
        .collect();
    let effect_of = |i: usize| -> String {
        let d = &diffs[i];
        match &snp_ann[i] {
            Some(a) => a.effect.label().to_string(),
            None => protein::indel_effect(d.bases1.len().max(d.bases2.len()), d.kind == "INS")
                .to_string(),
        }
    };
    let (prot1, prot2) = (code.translate(&s1), code.translate(&s2));

    // protein alignment, exactly as the aa-* distance modes compute it
    let psettings = ProteinSettings {
        translation_table: args.translation_table,
        first_codon_as_met: !args.no_first_codon_as_met,
        matrix: args.aa_matrix.clone(),
        gap_open: args.aa_gap_open,
        gap_extend: args.aa_gap_extend,
    };
    psettings.validate()?;
    let (pp1, pp2) = (protein_of(&code, &s1), protein_of(&code, &s2));
    let (paln, pq, pr) = align_proteins_with_strings(&psettings, &pp1, &pp2)
        .ok_or("the proteins could not be aligned")?;
    let pdiffs = if pp1 == pp2 {
        Vec::new()
    } else {
        differences(pq.as_bytes(), pr.as_bytes())
    };

    if args.tsv {
        println!("type\tevent\tallele1_start\tallele1_end\tallele2_start\tallele2_end\tallele1_bases\tallele2_bases\teffect\tcodon1\tcodon2\tprotein_change");
        for (i, d) in diffs.iter().enumerate() {
            let (c1, c2, pc) = match &snp_ann[i] {
                Some(a) => (
                    format!("{}:{}", a.codon1, a.codon_seq1),
                    format!("{}:{}", a.codon2, a.codon_seq2),
                    a.protein_change(),
                ),
                None => ("-".into(), "-".into(), "-".into()),
            };
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                d.kind,
                d.event,
                d.a1.0,
                d.a1.1,
                d.a2.0,
                d.a2.1,
                if d.bases1.is_empty() { "-" } else { &d.bases1 },
                if d.bases2.is_empty() { "-" } else { &d.bases2 },
                effect_of(i),
                c1,
                c2,
                pc
            );
        }
        if args.protein_diffs {
            for d in &pdiffs {
                let kind = match d.kind {
                    "SNP" => "AA_SUB",
                    "INS" => "AA_INS",
                    _ => "AA_DEL",
                };
                let change = if d.kind == "SNP" {
                    format!(
                        "p.{}{}{}",
                        protein::aa_name(d.bases1.as_bytes()[0]),
                        d.a1.0,
                        protein::aa_name(d.bases2.as_bytes()[0])
                    )
                } else {
                    "-".to_string()
                };
                println!(
                    "{kind}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t-\t-\t-\t{change}",
                    d.event,
                    d.a1.0,
                    d.a1.1,
                    d.a2.0,
                    d.a2.1,
                    if d.bases1.is_empty() { "-" } else { &d.bases1 },
                    if d.bases2.is_empty() { "-" } else { &d.bases2 }
                );
            }
        }
    } else {
        println!("locus      {}", args.locus);
        println!("allele 1   hash {lo}  {id1}  {} bp", s1.len());
        println!("allele 2   hash {hi}  {id2}  {} bp", s2.len());
        println!(
            "scoring    match={} mismatch={} gap_open={} gap_extend={}",
            config.match_score, config.mismatch_penalty, config.gap_open, config.gap_extend
        );
        println!(
            "counts     snps={}  indel_events={}  indel_bases={}  score={}",
            aln.snps, aln.indel_events, aln.indel_bases, aln.score
        );
        let ins: usize = diffs
            .iter()
            .filter(|d| d.kind == "INS")
            .map(|d| d.bases2.len())
            .sum();
        let del: usize = diffs
            .iter()
            .filter(|d| d.kind == "DEL")
            .map(|d| d.bases1.len())
            .sum();
        println!(
            "           bases only in allele 2 (INS): {ins}   bases only in allele 1 (DEL): {del}"
        );
        let count = |e: SnpEffect| snp_ann.iter().flatten().filter(|a| a.effect == e).count();
        let syn = snp_ann
            .iter()
            .flatten()
            .filter(|a| a.effect.is_synonymous())
            .count();
        let nonsyn = snp_ann
            .iter()
            .flatten()
            .filter(|a| a.effect.is_nonsynonymous())
            .count();
        println!(
            "code       NCBI table {} ({}), first codon as Met: {}",
            code.table_id(),
            code.table_name(),
            if code.first_codon_as_met { "yes" } else { "no" }
        );
        println!(
            "coding     synonymous={syn} (synonymous_variant={} stop_retained_variant={} start_retained_variant={})",
            count(SnpEffect::Synonymous),
            count(SnpEffect::StopRetained),
            count(SnpEffect::StartRetained)
        );
        println!(
            "           nonsynonymous={nonsyn} (missense_variant={} stop_gained={} stop_lost={} start_lost={})",
            count(SnpEffect::Missense),
            count(SnpEffect::StopGained),
            count(SnpEffect::StopLost),
            count(SnpEffect::StartLost)
        );
        println!(
            "           frame-disrupted (coding_sequence_variant)={}",
            count(SnpEffect::FrameDisrupted)
        );
        println!(
            "protein    allele1 {} aa  allele2 {} aa  {}",
            prot1.len(),
            prot2.len(),
            if prot1 == prot2 {
                "identical proteins"
            } else {
                "proteins differ"
            }
        );
        println!(
            "protein    alignment ({}, gaps {}/{}): aa_subs={}  aa_indel_events={}  aa_indel_residues={}",
            psettings.matrix,
            psettings.gap_open,
            psettings.gap_extend,
            paln.aa_subs,
            paln.aa_indel_events,
            paln.aa_indel_residues
        );
        println!("cigar      {cigar}");
        println!(
            "           (INS/DEL are relative to allele 1: INS = bases only in allele 2, \
             DEL = bases only in allele 1; in the CIGAR, allele 1 is the query, so \
             D = base only in allele 2, I = base only in allele 1)"
        );
        if !diffs.is_empty() {
            println!();
        }
        for (i, d) in diffs.iter().enumerate() {
            match d.kind {
                "SNP" => {
                    let a = snp_ann[i].as_ref().unwrap();
                    println!(
                        "SNP        allele1 {:>6} {}  ->  allele2 {:>6} {}   codon {} {}>{}  {:<16} {}",
                        d.a1.0,
                        d.bases1,
                        d.a2.0,
                        d.bases2,
                        a.codon1,
                        a.codon_seq1,
                        a.codon_seq2,
                        a.effect.label(),
                        a.protein_change()
                    )
                }
                "INS" => println!(
                    "INS #{:<4}  after allele1 {:>6}  allele2 {}-{}  +{} ({} bp, {})",
                    d.event,
                    d.a1.0,
                    d.a2.0,
                    d.a2.1,
                    d.bases2,
                    d.bases2.len(),
                    effect_of(i)
                ),
                _ => println!(
                    "DEL #{:<4}  allele1 {}-{}  after allele2 {:>6}  -{} ({} bp, {})",
                    d.event,
                    d.a1.0,
                    d.a1.1,
                    d.a2.0,
                    d.bases1,
                    d.bases1.len(),
                    effect_of(i)
                ),
            }
        }
    }

    if args.protein_diffs && !args.tsv && !pdiffs.is_empty() {
        println!();
        for d in &pdiffs {
            match d.kind {
                "SNP" => println!(
                    "AA_SUB     protein1 {:>5} {}  ->  protein2 {:>5} {}   p.{}{}{}",
                    d.a1.0,
                    d.bases1,
                    d.a2.0,
                    d.bases2,
                    protein::aa_name(d.bases1.as_bytes()[0]),
                    d.a1.0,
                    protein::aa_name(d.bases2.as_bytes()[0])
                ),
                "INS" => println!(
                    "AA_INS #{:<3} after protein1 {:>5}  protein2 {}-{}  +{} ({} aa)",
                    d.event,
                    d.a1.0,
                    d.a2.0,
                    d.a2.1,
                    d.bases2,
                    d.bases2.len()
                ),
                _ => println!(
                    "AA_DEL #{:<3} protein1 {}-{}  after protein2 {:>5}  -{} ({} aa)",
                    d.event,
                    d.a1.0,
                    d.a1.1,
                    d.a2.0,
                    d.bases1,
                    d.bases1.len()
                ),
            }
        }
    }

    if let Some(cache) = &args.cache_file {
        let key = format!("{}:{lo}:{hi}", args.locus);
        match cache_entry(cache, &key)? {
            None => eprintln!("cache      no entry {key} in {cache}"),
            Some((snps, ev, bases, cfg)) => {
                let same_params = (
                    cfg.match_score,
                    cfg.mismatch_penalty,
                    cfg.gap_open,
                    cfg.gap_extend,
                ) == (
                    config.match_score,
                    config.mismatch_penalty,
                    config.gap_open,
                    config.gap_extend,
                );
                let ok = (snps, ev, bases) == (aln.snps, aln.indel_events, aln.indel_bases);
                eprintln!(
                    "cache      snps={snps}  indel_events={ev}  indel_bases={bases}  -> {}{}",
                    if ok { "MATCH" } else { "MISMATCH" },
                    if same_params {
                        ""
                    } else {
                        " (cache built with different alignment parameters)"
                    }
                );
                if !ok && same_params {
                    return Ok(2);
                }
            }
        }
    }
    Ok(0)
}

fn main() {
    let args: Args = argh::from_env();
    match run(args) {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            eprintln!("❌ ERROR: {e}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::differences;

    #[test]
    fn lists_snps_and_indels_with_positions() {
        // allele1: ATGCA-TTG, allele2: ATCCATT-G
        let d = differences(b"ATGCA-TTG", b"ATCCATT-G");
        assert_eq!(d.len(), 3);
        assert_eq!(
            (
                d[0].kind,
                d[0].a1,
                d[0].a2,
                d[0].bases1.as_str(),
                d[0].bases2.as_str()
            ),
            ("SNP", (3, 3), (3, 3), "G", "C")
        );
        // insertion in allele 2 after allele-1 position 5, at allele-2 position 6
        assert_eq!(
            (
                d[1].kind,
                d[1].event,
                d[1].a1,
                d[1].a2,
                d[1].bases2.as_str()
            ),
            ("INS", 1, (5, 5), (6, 6), "T")
        );
        // deletion: allele-1 position 7 missing in allele 2 (after allele-2 position 7)
        assert_eq!(
            (
                d[2].kind,
                d[2].event,
                d[2].a1,
                d[2].a2,
                d[2].bases1.as_str()
            ),
            ("DEL", 2, (7, 7), (7, 7), "T")
        );
    }

    #[test]
    fn adjacent_gap_kinds_are_one_event() {
        // an insertion directly followed by a deletion is one InDel event
        let d = differences(b"A-CG", b"AT-G");
        assert_eq!(d.iter().map(|x| x.event).collect::<Vec<_>>(), vec![1, 1]);
    }
}
