// protein.rs - Coding effect of allele differences (translation table 11)
//
// cgMLST alleles (chewBBACA) are complete coding sequences in frame: length
// a multiple of 3, a start codon, a stop codon at the end. This module
// translates them with the bacterial/archaeal/plant-plastid code (NCBI
// translation table 11) and classifies each difference found by an alignment:
// SNPs by the change of the whole codon they fall in, InDels as in-frame or
// frameshift. It only annotates; no distance or cache value depends on it.

/// Start codons of translation table 11. At the first codon they are
/// translated as Met.
const STARTS_11: [&[u8; 3]; 7] = [b"ATG", b"GTG", b"TTG", b"CTG", b"ATT", b"ATC", b"ATA"];

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Standard codon table (identical amino acids in table 11), indexed by
/// TCAG order: first base * 16 + second * 4 + third.
const CODE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// Amino acid of a codon ('X' if it has a non-ACGT base). `first` applies
/// the initiation rule: any table-11 start codon is Met.
pub fn translate_codon(codon: &[u8], first: bool) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    if first {
        let up = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        if STARTS_11.iter().any(|s| **s == up) {
            return b'M';
        }
    }
    match (
        base_index(codon[0]),
        base_index(codon[1]),
        base_index(codon[2]),
    ) {
        (Some(a), Some(b), Some(c)) => CODE[a * 16 + b * 4 + c],
        _ => b'X',
    }
}

/// Protein of a coding sequence (table 11, first codon as initiator).
/// A trailing incomplete codon is ignored.
pub fn translate(cds: &[u8]) -> Vec<u8> {
    cds.chunks_exact(3)
        .enumerate()
        .map(|(i, c)| translate_codon(c, i == 0))
        .collect()
}

pub fn aa_name(aa: u8) -> &'static str {
    match aa {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'Q' => "Gln",
        b'E' => "Glu",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        b'*' => "Ter",
        _ => "Xaa",
    }
}

/// Effect of a SNP on the protein.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnpEffect {
    /// same amino acid
    Synonymous,
    /// different amino acid
    Missense,
    /// amino acid -> stop (premature stop codon)
    Nonsense,
    /// stop -> amino acid
    StopLost,
    /// the initiator codon is no longer a start codon
    StartLost,
    /// the codon is shifted or split by an InDel: no codon-to-codon comparison
    FrameDisrupted,
}

impl SnpEffect {
    pub fn label(self) -> &'static str {
        match self {
            SnpEffect::Synonymous => "synonymous",
            SnpEffect::Missense => "missense",
            SnpEffect::Nonsense => "nonsense",
            SnpEffect::StopLost => "stop_lost",
            SnpEffect::StartLost => "start_lost",
            SnpEffect::FrameDisrupted => "frame_disrupted",
        }
    }

    /// Changes the encoded protein (everything except synonymous and
    /// frame-disrupted, which is not classifiable at codon level).
    pub fn is_nonsynonymous(self) -> bool {
        matches!(
            self,
            SnpEffect::Missense | SnpEffect::Nonsense | SnpEffect::StopLost | SnpEffect::StartLost
        )
    }
}

/// Codon-level annotation of one SNP.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnpAnnotation {
    pub effect: SnpEffect,
    /// 1-based codon (= amino-acid) number in allele 1 and allele 2
    pub codon1: usize,
    pub codon2: usize,
    pub codon_seq1: String,
    pub codon_seq2: String,
    pub aa1: u8,
    pub aa2: u8,
}

impl SnpAnnotation {
    /// HGVS-like protein notation, numbered on allele 1: p.Gly8Ser,
    /// p.Gly8= (synonymous), p.Trp20Ter (nonsense).
    pub fn protein_change(&self) -> String {
        if self.effect == SnpEffect::FrameDisrupted {
            return "p.?".to_string();
        }
        let a = aa_name(self.aa1);
        if self.aa1 == self.aa2 {
            format!("p.{a}{}=", self.codon1)
        } else {
            format!("p.{a}{}{}", self.codon1, aa_name(self.aa2))
        }
    }
}

/// Column-by-column mapping of an alignment: for each base of allele 1
/// (0-based), the aligned base of allele 2, if any.
pub fn map_allele1_to_allele2(query_aln: &[u8], ref_aln: &[u8]) -> Vec<Option<usize>> {
    let mut map = Vec::new();
    let (mut p1, mut p2) = (0usize, 0usize);
    for (&a, &b) in query_aln.iter().zip(ref_aln) {
        match (a == b'-', b == b'-') {
            (false, false) => {
                map.push(Some(p2));
                p1 += 1;
                p2 += 1;
            }
            (false, true) => {
                map.push(None);
                p1 += 1;
            }
            (true, false) => p2 += 1,
            (true, true) => {}
        }
    }
    debug_assert_eq!(map.len(), p1);
    map
}

/// Classify the SNP at 0-based positions (p1 in allele 1, p2 in allele 2).
/// The whole codon of allele 1 must be aligned, gap-free and in the same
/// phase, to one codon of allele 2; otherwise the SNP is FrameDisrupted.
/// Several SNPs in one codon all receive that codon's effect.
pub fn classify_snp(
    seq1: &[u8],
    seq2: &[u8],
    map: &[Option<usize>],
    p1: usize,
    p2: usize,
) -> SnpAnnotation {
    let (c1, c2) = (p1 / 3, p2 / 3);
    let aligned = p1 % 3 == p2 % 3
        && (0..3).all(|k| map.get(c1 * 3 + k).copied().flatten() == Some(c2 * 3 + k))
        && c1 * 3 + 3 <= seq1.len()
        && c2 * 3 + 3 <= seq2.len();
    let codon = |s: &[u8], c: usize| -> String {
        s.get(c * 3..(c * 3 + 3).min(s.len()))
            .map(|x| String::from_utf8_lossy(x).into_owned())
            .unwrap_or_default()
    };
    let (cs1, cs2) = (codon(seq1, c1), codon(seq2, c2));
    let (aa1, aa2) = (
        translate_codon(cs1.as_bytes(), c1 == 0),
        translate_codon(cs2.as_bytes(), c2 == 0),
    );
    let effect = if !aligned {
        SnpEffect::FrameDisrupted
    } else if c1 == 0 && aa1 == b'M' && aa2 != b'M' {
        SnpEffect::StartLost
    } else if aa1 == aa2 {
        SnpEffect::Synonymous
    } else if aa2 == b'*' {
        SnpEffect::Nonsense
    } else if aa1 == b'*' {
        SnpEffect::StopLost
    } else {
        SnpEffect::Missense
    };
    SnpAnnotation {
        effect,
        codon1: c1 + 1,
        codon2: c2 + 1,
        codon_seq1: cs1,
        codon_seq2: cs2,
        aa1,
        aa2,
    }
}

/// Effect of an InDel run of `len` bases.
pub fn indel_effect(len: usize) -> &'static str {
    if len.is_multiple_of(3) {
        "in_frame"
    } else {
        "frameshift"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translation_table_11() {
        assert_eq!(translate(b"ATGGGATAA"), b"MG*");
        // alternative start codons are Met only as initiators
        assert_eq!(translate(b"GTGGTGTAA"), b"MV*");
        assert_eq!(translate(b"TTGTTGTGA"), b"ML*");
        assert_eq!(translate(b"ATTATTTAG"), b"MI*");
        assert_eq!(translate_codon(b"NNN", false), b'X');
        assert_eq!(aa_name(b'*'), "Ter");
    }

    fn classify(a: &[u8], b: &[u8], p: usize) -> SnpAnnotation {
        let map = map_allele1_to_allele2(a, b);
        classify_snp(a, b, &map, p, p)
    }

    #[test]
    fn snp_effects() {
        // GGA (Gly) -> GGG (Gly)
        let s = classify(b"ATGGGATAA", b"ATGGGGTAA", 5);
        assert_eq!(
            (s.effect, s.protein_change().as_str()),
            (SnpEffect::Synonymous, "p.Gly2=")
        );
        // GGA (Gly) -> AGA (Arg)
        let s = classify(b"ATGGGATAA", b"ATGAGATAA", 3);
        assert_eq!(
            (s.effect, s.protein_change().as_str()),
            (SnpEffect::Missense, "p.Gly2Arg")
        );
        // TGG (Trp) -> TGA (stop)
        assert_eq!(
            classify(b"ATGTGGTAA", b"ATGTGATAA", 5).effect,
            SnpEffect::Nonsense
        );
        // TAA (stop) -> CAA (Gln)
        assert_eq!(
            classify(b"ATGGGATAA", b"ATGGGACAA", 6).effect,
            SnpEffect::StopLost
        );
        // ATG -> ACG: no longer a start codon
        assert_eq!(
            classify(b"ATGGGATAA", b"ACGGGATAA", 1).effect,
            SnpEffect::StartLost
        );
        // ATG -> GTG: still a (table 11) start codon, Met
        assert_eq!(
            classify(b"ATGGGATAA", b"GTGGGATAA", 0).effect,
            SnpEffect::Synonymous
        );
    }

    #[test]
    fn shifted_codons_are_frame_disrupted() {
        // allele 2 has one extra base before the SNP: codons are shifted
        let a = b"ATGGGA-CCATAA";
        let b = b"ATGGGATCGATAA";
        let seq1: Vec<u8> = a.iter().copied().filter(|&c| c != b'-').collect();
        let seq2: Vec<u8> = b.to_vec();
        let map = map_allele1_to_allele2(a, b);
        // SNP at alignment column 8: allele1 base 7 (C), allele2 base 8 (G)
        let s = classify_snp(&seq1, &seq2, &map, 7, 8);
        assert_eq!(s.effect, SnpEffect::FrameDisrupted);
        assert_eq!(indel_effect(1), "frameshift");
        assert_eq!(indel_effect(6), "in_frame");
    }
}
