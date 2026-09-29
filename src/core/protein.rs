// protein.rs - Coding effect of allele differences (any NCBI genetic code)
//
// cgMLST alleles (chewBBACA) are complete coding sequences in frame: length
// a multiple of 3, a start codon, a stop codon at the end. This module
// translates them with a configurable NCBI translation table (default 11,
// bacterial) and classifies each difference found by an alignment:
// SNPs by the change of the whole codon they fall in, InDels as in-frame or
// frameshift. It only annotates; no distance or cache value depends on it.

use crate::core::codon_tables::{self, CodonTableDef};

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// A genetic code: an NCBI translation table plus the initiation rule.
#[derive(Clone, Copy)]
pub struct GeneticCode {
    table: &'static CodonTableDef,
    /// Translate the first codon as Met when it is a start codon of the
    /// table (as translation does at initiation: GTG/TTG start -> Met).
    pub first_codon_as_met: bool,
}

impl std::fmt::Debug for GeneticCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GeneticCode({} {})", self.table.id, self.table.name)
    }
}

impl Default for GeneticCode {
    /// NCBI table 11 (Bacterial, Archaeal and Plant Plastid), first codon Met.
    fn default() -> Self {
        Self::new(11, true).expect("table 11 exists")
    }
}

impl GeneticCode {
    pub fn new(table_id: u32, first_codon_as_met: bool) -> Result<Self, String> {
        let table = codon_tables::table(table_id).ok_or_else(|| {
            let ids: Vec<String> = codon_tables::TABLES
                .iter()
                .map(|t| t.id.to_string())
                .collect();
            format!(
                "unknown translation table {table_id} (available: {})",
                ids.join(", ")
            )
        })?;
        Ok(Self {
            table,
            first_codon_as_met,
        })
    }

    pub fn table_id(&self) -> u32 {
        self.table.id
    }

    pub fn table_name(&self) -> &'static str {
        self.table.name
    }

    fn index(codon: &[u8]) -> Option<usize> {
        if codon.len() != 3 {
            return None;
        }
        Some(base_index(codon[0])? * 16 + base_index(codon[1])? * 4 + base_index(codon[2])?)
    }

    /// Whether a codon is an initiation codon of this table.
    pub fn is_start(&self, codon: &[u8]) -> bool {
        Self::index(codon).is_some_and(|i| self.table.starts[i] == b'M')
    }

    /// Amino acid of a codon ('*' stop, 'X' if it has a non-ACGT base).
    /// `first` applies the initiation rule to the first codon of a CDS.
    pub fn translate_codon(&self, codon: &[u8], first: bool) -> u8 {
        match Self::index(codon) {
            None => b'X',
            Some(i) => {
                if first && self.first_codon_as_met && self.table.starts[i] == b'M' {
                    b'M'
                } else {
                    self.table.aa[i]
                }
            }
        }
    }

    /// Protein of a coding sequence; a trailing incomplete codon is ignored.
    pub fn translate(&self, cds: &[u8]) -> Vec<u8> {
        cds.chunks_exact(3)
            .enumerate()
            .map(|(i, c)| self.translate_codon(c, i == 0))
            .collect()
    }
}

/// Amino acid of a codon with the default code (table 11, first codon Met).
pub fn translate_codon(codon: &[u8], first: bool) -> u8 {
    GeneticCode::default().translate_codon(codon, first)
}

/// Protein of a CDS with the default code (table 11, first codon Met).
pub fn translate(cds: &[u8]) -> Vec<u8> {
    GeneticCode::default().translate(cds)
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

/// Effect of a SNP on the protein, as Sequence Ontology consequence terms
/// (the vocabulary of snpEff, Ensembl VEP and bcftools csq).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnpEffect {
    /// same amino acid (SO:0001819 synonymous_variant)
    Synonymous,
    /// stop codon changed to another stop (SO:0001567 stop_retained_variant)
    StopRetained,
    /// initiator codon changed to another start codon (SO:0002019
    /// start_retained_variant)
    StartRetained,
    /// different amino acid (SO:0001583 missense_variant)
    Missense,
    /// amino acid -> stop (SO:0001587 stop_gained)
    StopGained,
    /// stop -> amino acid (SO:0001578 stop_lost)
    StopLost,
    /// initiator codon no longer a start codon (SO:0002012 start_lost)
    StartLost,
    /// the codon is shifted or split by an InDel, so no codon-to-codon
    /// comparison is possible (SO:0001580 coding_sequence_variant)
    FrameDisrupted,
}

impl SnpEffect {
    /// Sequence Ontology term.
    pub fn label(self) -> &'static str {
        match self {
            SnpEffect::Synonymous => "synonymous_variant",
            SnpEffect::StopRetained => "stop_retained_variant",
            SnpEffect::StartRetained => "start_retained_variant",
            SnpEffect::Missense => "missense_variant",
            SnpEffect::StopGained => "stop_gained",
            SnpEffect::StopLost => "stop_lost",
            SnpEffect::StartLost => "start_lost",
            SnpEffect::FrameDisrupted => "coding_sequence_variant",
        }
    }

    /// The encoded amino acid is unchanged.
    pub fn is_synonymous(self) -> bool {
        matches!(
            self,
            SnpEffect::Synonymous | SnpEffect::StopRetained | SnpEffect::StartRetained
        )
    }

    /// The encoded protein changes (missense, stop gained/lost, start lost).
    /// Frame-disrupted SNPs are in neither class.
    pub fn is_nonsynonymous(self) -> bool {
        matches!(
            self,
            SnpEffect::Missense
                | SnpEffect::StopGained
                | SnpEffect::StopLost
                | SnpEffect::StartLost
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
    code: &GeneticCode,
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
        code.translate_codon(cs1.as_bytes(), c1 == 0),
        code.translate_codon(cs2.as_bytes(), c2 == 0),
    );
    let initiator = c1 == 0 && c2 == 0 && code.is_start(cs1.as_bytes());
    let effect = if !aligned {
        SnpEffect::FrameDisrupted
    } else if initiator && !code.is_start(cs2.as_bytes()) {
        SnpEffect::StartLost
    } else if initiator {
        SnpEffect::StartRetained
    } else if aa1 == aa2 && aa1 == b'*' {
        SnpEffect::StopRetained
    } else if aa1 == aa2 {
        SnpEffect::Synonymous
    } else if aa2 == b'*' {
        SnpEffect::StopGained
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

/// Sequence Ontology term of an InDel run of `len` bases; `insertion` is
/// relative to allele 1 (bases present only in allele 2).
pub fn indel_effect(len: usize, insertion: bool) -> &'static str {
    match (len.is_multiple_of(3), insertion) {
        (false, _) => "frameshift_variant",
        (true, true) => "inframe_insertion",
        (true, false) => "inframe_deletion",
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
        let no_met = GeneticCode::new(11, false).unwrap();
        assert_eq!(no_met.translate(b"GTGGTGTAA"), b"VV*");
    }

    #[test]
    fn other_tables() {
        // table 4 (Mycoplasma): TGA = Trp
        assert_eq!(
            GeneticCode::new(4, true).unwrap().translate(b"ATGTGATAA"),
            b"MW*"
        );
        // table 1: TGA = stop
        assert_eq!(
            GeneticCode::new(1, true).unwrap().translate(b"ATGTGATAA"),
            b"M**"
        );
        assert!(GeneticCode::new(7, true).is_err());
    }

    fn classify(a: &[u8], b: &[u8], p: usize) -> SnpAnnotation {
        let map = map_allele1_to_allele2(a, b);
        classify_snp(&GeneticCode::default(), a, b, &map, p, p)
    }

    #[test]
    fn snp_effects() {
        let s = classify(b"ATGGGATAA", b"ATGGGGTAA", 5);
        assert_eq!(
            (s.effect, s.protein_change().as_str()),
            (SnpEffect::Synonymous, "p.Gly2=")
        );
        let s = classify(b"ATGGGATAA", b"ATGAGATAA", 3);
        assert_eq!(
            (s.effect, s.protein_change().as_str()),
            (SnpEffect::Missense, "p.Gly2Arg")
        );
        assert_eq!(
            classify(b"ATGTGGTAA", b"ATGTGATAA", 5).effect,
            SnpEffect::StopGained
        );
        assert_eq!(
            classify(b"ATGGGATAA", b"ATGGGACAA", 6).effect,
            SnpEffect::StopLost
        );
        assert_eq!(
            classify(b"ATGGGATAA", b"ATGGGATAG", 8).effect,
            SnpEffect::StopRetained
        );
        assert_eq!(
            classify(b"ATGGGATAA", b"ACGGGATAA", 1).effect,
            SnpEffect::StartLost
        );
        // ATG -> GTG: still a (table 11) start codon
        assert_eq!(
            classify(b"ATGGGATAA", b"GTGGGATAA", 0).effect,
            SnpEffect::StartRetained
        );
        assert!(SnpEffect::StartRetained.is_synonymous());
        assert!(SnpEffect::StopGained.is_nonsynonymous());
        assert!(
            !SnpEffect::FrameDisrupted.is_synonymous()
                && !SnpEffect::FrameDisrupted.is_nonsynonymous()
        );
    }

    #[test]
    fn shifted_codons_are_frame_disrupted() {
        let a = b"ATGGGA-CCATAA";
        let b = b"ATGGGATCGATAA";
        let seq1: Vec<u8> = a.iter().copied().filter(|&c| c != b'-').collect();
        let map = map_allele1_to_allele2(a, b);
        let s = classify_snp(&GeneticCode::default(), &seq1, b, &map, 7, 8);
        assert_eq!(s.effect, SnpEffect::FrameDisrupted);
        assert_eq!(s.effect.label(), "coding_sequence_variant");
        assert_eq!(indel_effect(1, true), "frameshift_variant");
        assert_eq!(indel_effect(6, true), "inframe_insertion");
        assert_eq!(indel_effect(3, false), "inframe_deletion");
    }
}
