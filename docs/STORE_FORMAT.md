# cgDist cache store format

A **cache store** holds precomputed alignment statistics for pairs of alleles
of a cgMLST/wgMLST schema, one file per locus, plus a JSON manifest. A store
can be a directory or a single `.cgpack` file; both are readable over HTTP.

This document is normative for store format version 1 with locus files of
version 1 and 2. `scripts/cgds_reader.py` is a reference reader written from
it (Python standard library only), and `scripts/test_store_format.py` checks
that it decodes every locus exactly as cgdist does.

All integers are unsigned. "varint" is LEB128: 7 bits per byte, least
significant group first, high bit set on every byte but the last; at most
10 bytes; a value that must fit in u32 and does not is an error.

## 1. Directory store

```
<store>/manifest.json
<store>/loci/<locus>.cgds
<store>/.lock                 (present only while a process writes the store)
```

Locus names contain only ASCII letters, digits, `_`, `-`, `.`, `+`, and do not
start with `.`; the file of locus `L` is `loci/L.cgds`.

## 2. Manifest (`manifest.json`)

```json
{
  "format": "cgdist-store",
  "format_version": 1,
  "hasher": "crc32",
  "alignment": {"match_score": 2, "mismatch_penalty": -1, "gap_open": 5, "gap_extend": 2},
  "genetic_code": {"table": 11, "first_codon_as_met": true},
  "schema": {"name": "...", "source": "chewie-ns", "version": "..."},
  "note": "...",
  "cgdist_version": "0.1.4",
  "created": "2026-09-29T14:43:04Z",
  "last_modified": "2026-09-29T15:22:41Z",
  "loci": {
    "cgMLST-00085850": {
      "file": "loci/cgMLST-00085850.cgds", "alleles": 2, "pairs": 1,
      "complete": true, "bytes": 35, "sha256": "378e…", "offset": 0
    }
  }
}
```

* `format` must be `cgdist-store`; a reader refuses a `format_version` above
  the one it knows.
* Exactly one of `alignment` (**DNA store**) and `protein` (**protein store**)
  is present. A protein store has
  `"protein": {"translation_table": 11, "first_codon_as_met": true, "matrix": "blosum62", "gap_open": 11, "gap_extend": 1}`
  and no `genetic_code`.
* `alignment` holds the four scoring numbers of the global alignment (a gap of
  L bases costs `gap_open + (L-1)·gap_extend`). Two stores are interchangeable
  only if these numbers, the hasher and (for protein stores) every protein
  setting are equal; readers must refuse a store whose parameters differ from
  the run's.
* `schema` describes where the alleles come from and does not affect the
  results. Besides `name`, `source` and `version` it may carry any other field,
  e.g. `url`, `species_id`, `schema_id`, `nr_loci`, `nr_alleles`, `citation`;
  tools keep these fields when they copy a manifest (`pull`, `pack`).
* `genetic_code`, when present, is the code of the synonymous/nonsynonymous
  counts stored in the locus files.
* Per locus: `alleles` and `pairs` are the counts in the locus file,
  `complete` is true when every pair of its alleles is present
  (`pairs == alleles·(alleles-1)/2`), `bytes` and `sha256` describe the locus
  file, and `offset` (packs only) is its position in the pack.
* `schema`, `note`, `cgdist_version`, `created`, `last_modified` are
  informative.

## 3. Locus file (`.cgds`)

```
offset 0  4 bytes  magic "CGDS"
offset 4  1 byte   version: 1 or 2
offset 5  ...      body, compressed: u32 LE uncompressed size, then one LZ4 block
```

The body (after decompression) is, in order:

| field | encoding |
|---|---|
| `n` | varint: number of alleles |
| allele hashes | `n` varints: the first hash, then the difference to the previous hash; strictly increasing, so every difference after the first is > 0 |
| allele lengths | `n` varints: sequence length of each allele in table order, 0 = unknown |
| `m` | varint: number of pairs |
| pair first indices | `m` varints: delta-encoded allele index `i` of each pair (non-decreasing) |
| pair second indices | `m` varints: for the first pair of each `i`, `j - i - 1`; for the next pairs with the same `i`, `j - j_prev - 1` |
| SNPs | `m` varints |
| InDel events | `m` varints |
| InDel bases | `m` varints |
| coding flag | 1 byte: 0 = no coding columns, 1 = three columns follow |
| synonymous, nonsynonymous, frame-disrupted | if flag = 1: three columns of `m` varints each, value + 1 (0 = not computed for this pair; either all three are 0 or all three are > 0) |
| digests (version 2 only) | `n` × u64 little-endian, in allele table order: sequence digest of each allele, 0 = unknown |

Nothing may follow. Pairs are listed in increasing `(i, j)` order with
`i < j < n`; the pair refers to the alleles with hashes `hash[i] < hash[j]`.
Identical alleles have distance 0 and are never stored.

* In a **DNA store** the hash is the CRC32 of the allele sequence (as in
  chewBBACA hashed profiles), the length is in bases, and the three
  statistics are the SNPs, InDel events and InDel bases of the global
  alignment of the two alleles (allele with the lower hash as query).
* In a **protein store** the hash is the CRC32 of the translated protein
  (terminal stop codon removed), the length is in amino acids, and the three
  statistics are amino-acid substitutions, InDel events and InDel residues.
  There are no coding columns.
* **Digest** of a sequence: the first 8 bytes of its SHA-256, read as u64
  little-endian; the value 0 is replaced by 1. It lets a reader tell apart two
  different sequences with the same CRC32: an allele of a run whose sequence
  digest differs from the stored digest (or, for alleles without digest,
  whose length differs from a known stored length) must not take its pairs
  from the store.

Writers emit version 1 when no digest is known and version 2 otherwise.

## 4. Pack file (`.cgpack`)

```
offset 0   8 bytes  magic "CGDPACK1"
offset 8   u64 LE   length L of the manifest
offset 16  L bytes  manifest (JSON, as above, with "offset" in every locus entry)
offset 16+L         locus files concatenated; locus X starts at 16 + L + offset(X)
                    and is bytes(X) long
```

A reader can fetch the first 16 bytes, then the manifest, then each needed
locus with an HTTP Range request. Every locus file must match its `sha256`.

## 5. Merging and updating

A pair present in two stores built with the same parameters has the same
statistics; a difference, or two different digests for the same hash, means
the stores come from different computations or schemas, and tools must refuse
to merge them. `cgdist-cache pull` merges a downloaded locus with local
content instead of replacing it.

## 6. Compatibility rules

* Readers accept locus versions they know (1, 2) and refuse others.
* New optional manifest fields may be added; readers ignore unknown fields.
* Any change to the locus body layout gets a new locus version.
