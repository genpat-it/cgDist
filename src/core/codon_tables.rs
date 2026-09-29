// codon_tables.rs - NCBI genetic codes (translation tables)
//
// Generated from Biopython 1.86 Bio.Data.CodonTable.unambiguous_dna_by_id
// (source: NCBI gc.prt). Codons are in TCAG order: index = 16*b1 + 4*b2 + b3
// with T=0, C=1, A=2, G=3. `aa` holds the amino acid ('*' = stop) and
// `starts` marks initiation codons with 'M'. In tables 27, 28 and 31 some
// codons are stop or sense depending on context; they are listed as stop.

/// One NCBI translation table.
pub struct CodonTableDef {
    pub id: u32,
    pub name: &'static str,
    pub aa: &'static [u8; 64],
    pub starts: &'static [u8; 64],
}

pub const TABLES: &[CodonTableDef] = &[
    CodonTableDef {
        id: 1,
        name: "Standard",
        aa: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M---------------M---------------M----------------------------",
    },
    CodonTableDef {
        id: 2,
        name: "Vertebrate Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        starts: b"--------------------------------MMMM---------------M------------",
    },
    CodonTableDef {
        id: 3,
        name: "Yeast Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------------------------------MM---------------M------------",
    },
    CodonTableDef {
        id: 4,
        name: "Mold Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--MM---------------M------------MMMM---------------M------------",
    },
    CodonTableDef {
        id: 5,
        name: "Invertebrate Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        starts: b"---M----------------------------MMMM---------------M------------",
    },
    CodonTableDef {
        id: 6,
        name: "Ciliate Nuclear",
        aa: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 9,
        name: "Echinoderm Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M---------------M------------",
    },
    CodonTableDef {
        id: 10,
        name: "Euplotid Nuclear",
        aa: b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 11,
        name: "Bacterial",
        aa: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M---------------M------------MMMM---------------M------------",
    },
    CodonTableDef {
        id: 12,
        name: "Alternative Yeast Nuclear",
        aa: b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-------------------M---------------M----------------------------",
    },
    CodonTableDef {
        id: 13,
        name: "Ascidian Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        starts: b"---M------------------------------MM---------------M------------",
    },
    CodonTableDef {
        id: 14,
        name: "Alternative Flatworm Mitochondrial",
        aa: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 15,
        name: "Blepharisma Macronuclear",
        aa: b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 16,
        name: "Chlorophycean Mitochondrial",
        aa: b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 21,
        name: "Trematode Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M---------------M------------",
    },
    CodonTableDef {
        id: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        aa: b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 23,
        name: "Thraustochytrium Mitochondrial",
        aa: b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------------------------M--M---------------M------------",
    },
    CodonTableDef {
        id: 24,
        name: "Pterobranchia Mitochondrial",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M---------------M---------------M---------------M------------",
    },
    CodonTableDef {
        id: 25,
        name: "Candidate Division SR1",
        aa: b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M-------------------------------M---------------M------------",
    },
    CodonTableDef {
        id: 26,
        name: "Pachysolen tannophilus Nuclear",
        aa: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-------------------M---------------M----------------------------",
    },
    CodonTableDef {
        id: 27,
        name: "Karyorelict Nuclear",
        aa: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 28,
        name: "Condylostoma Nuclear",
        aa: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 29,
        name: "Mesodinium Nuclear",
        aa: b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 30,
        name: "Peritrich Nuclear",
        aa: b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 31,
        name: "Blastocrithidia Nuclear",
        aa: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"-----------------------------------M----------------------------",
    },
    CodonTableDef {
        id: 32,
        name: "Balanophoraceae Plastid",
        aa: b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M---------------M------------MMMM---------------M------------",
    },
    CodonTableDef {
        id: 33,
        name: "Cephalodiscidae Mitochondrial",
        aa: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M---------------M---------------M---------------M------------",
    },
];

pub fn table(id: u32) -> Option<&'static CodonTableDef> {
    TABLES.iter().find(|t| t.id == id)
}
