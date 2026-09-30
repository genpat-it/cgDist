#!/usr/bin/env python3
"""Reference reader for cgdist cache stores (format described in
docs/STORE_FORMAT.md). Python standard library only, LZ4 included: it exists
to show that the format can be read without cgdist, and to test conformance.

    cgds_reader.py info   <store dir | .cgpack>
    cgds_reader.py verify <store dir | .cgpack>        sha256 + decoding + counts
    cgds_reader.py dump   <store dir | .cgpack> <locus>   TSV of the locus pairs
    cgds_reader.py summary <store dir | .cgpack>       per-locus JSON summaries
    cgds_reader.py digests <store dir | .cgpack> <schema dir>   (DNA stores)
                   check every stored sequence digest against the schema

Exit code 1 on any problem.
"""

import hashlib
import json
import struct
import sys
from pathlib import Path

LOCUS_MAGIC = b"CGDS"
LOCUS_VERSIONS = (1, 2)  # 2 adds a sequence digest per allele
PACK_MAGIC = b"CGDPACK1"


# ------------------------------------------------------------------ LZ4 block

def lz4_block_decompress(src, size):
    """Decompress one raw LZ4 block to exactly `size` bytes."""
    out = bytearray()
    i, n = 0, len(src)
    while i < n:
        token = src[i]
        i += 1
        lit = token >> 4
        if lit == 15:
            while True:
                b = src[i]
                i += 1
                lit += b
                if b != 255:
                    break
        out += src[i:i + lit]
        i += lit
        if i >= n:  # the last sequence has literals only
            break
        offset = src[i] | (src[i + 1] << 8)
        i += 2
        if offset == 0 or offset > len(out):
            raise ValueError("corrupt LZ4 block: bad match offset")
        mlen = token & 15
        if mlen == 15:
            while True:
                b = src[i]
                i += 1
                mlen += b
                if b != 255:
                    break
        mlen += 4
        start = len(out) - offset
        for k in range(mlen):  # byte by byte: matches may overlap
            out.append(out[start + k])
    if len(out) != size:
        raise ValueError(f"corrupt LZ4 block: {len(out)} bytes, expected {size}")
    return bytes(out)


def lz4_size_prepended(data):
    """lz4_flex::compress_prepend_size: u32 LE uncompressed size + block."""
    (size,) = struct.unpack_from("<I", data, 0)
    return lz4_block_decompress(data[4:], size)


# ------------------------------------------------------------------ locus file

class Reader:
    def __init__(self, buf):
        self.buf, self.pos = buf, 0

    def varint(self):
        v, shift = 0, 0
        while True:
            if self.pos >= len(self.buf):
                raise ValueError("corrupt locus file: truncated")
            b = self.buf[self.pos]
            self.pos += 1
            v |= (b & 0x7F) << shift
            if b < 0x80:
                return v
            shift += 7
            if shift >= 64:
                raise ValueError("corrupt locus file: varint too long")

    def u32(self):
        v = self.varint()
        if v > 0xFFFFFFFF:
            raise ValueError("corrupt locus file: value overflow")
        return v


def decode_locus(blob):
    """Decode a .cgds blob -> (alleles, pairs).

    alleles: list of (hash, length, digest) sorted by hash (length or digest
    0 = unknown; digest = first 8 bytes of SHA-256 of the sequence, LE);
    pairs: list of (hash_lo, hash_hi, snps, indel_events, indel_bases,
    coding) with coding = (syn, nonsyn, frame_disrupted) or None.
    In a protein store "hash" is the protein hash, "length" the amino-acid
    length and the three counts are substitutions, InDel events, InDel
    residues."""
    if len(blob) < 5 or blob[:4] != LOCUS_MAGIC:
        raise ValueError("not a cgdist locus file (bad magic)")
    version = blob[4]
    if version not in LOCUS_VERSIONS:
        raise ValueError(f"unsupported locus file version {version}")
    r = Reader(lz4_size_prepended(blob[5:]))
    n = r.varint()
    hashes, h = [], 0
    for k in range(n):
        d = r.varint()
        if k > 0 and d == 0:
            raise ValueError("corrupt locus file: allele table not strictly increasing")
        h = d if k == 0 else h + d
        hashes.append(h)
    lengths = [r.u32() for _ in range(n)]
    m = r.varint()
    first = []
    i = 0
    for _ in range(m):
        i += r.varint()
        first.append(i)
    second, last = [], None
    for i in first:
        d = r.varint()
        j = (last[1] + 1 + d) if (last is not None and last[0] == i) else (i + 1 + d)
        if i >= n or j >= n:
            raise ValueError("corrupt locus file: pair index out of range")
        second.append(j)
        last = (i, j)
    snps = [r.u32() for _ in range(m)]
    ev = [r.u32() for _ in range(m)]
    bases = [r.u32() for _ in range(m)]
    flag = r.buf[r.pos]
    r.pos += 1
    coding = [None] * m
    if flag == 1:
        cols = [[r.u32() for _ in range(m)] for _ in range(3)]
        for k in range(m):
            a, b, c = cols[0][k], cols[1][k], cols[2][k]
            if (a, b, c) == (0, 0, 0):
                continue
            if a == 0 or b == 0 or c == 0:
                raise ValueError("corrupt locus file: partial coding counts")
            coding[k] = (a - 1, b - 1, c - 1)
    elif flag != 0:
        raise ValueError("corrupt locus file: bad coding flag")
    digests = [0] * n
    if version == 2:
        need = 8 * n
        if r.pos + need > len(r.buf):
            raise ValueError("corrupt locus file: truncated digests")
        digests = list(struct.unpack_from(f"<{n}Q", r.buf, r.pos))
        r.pos += need
    if r.pos != len(r.buf):
        raise ValueError("corrupt locus file: trailing bytes")
    alleles = list(zip(hashes, lengths, digests))
    pairs = [(hashes[i], hashes[j], snps[k], ev[k], bases[k], coding[k])
             for k, (i, j) in enumerate(zip(first, second))]
    return alleles, pairs


# ------------------------------------------------------------------ stores

class Store:
    """A store directory or a .cgpack file."""

    def __init__(self, path):
        self.path = Path(path)
        if self.path.is_dir():
            self.pack = None
            self.manifest = json.loads((self.path / "manifest.json").read_text())
        else:
            self.pack = open(self.path, "rb")
            head = self.pack.read(16)
            if head[:8] != PACK_MAGIC:
                raise ValueError("not a cgdist cache pack (bad magic)")
            (mlen,) = struct.unpack("<Q", head[8:16])
            self.manifest = json.loads(self.pack.read(mlen))
            self.data_start = 16 + mlen
        if self.manifest.get("format") != "cgdist-store":
            raise ValueError("not a cgdist cache store")
        if self.manifest.get("format_version", 0) > 1:
            raise ValueError("store format newer than this reader")
        if ("alignment" in self.manifest) == ("protein" in self.manifest):
            raise ValueError("manifest must have exactly one of 'alignment' and 'protein'")

    def blob(self, locus):
        e = self.manifest["loci"][locus]
        if self.pack is None:
            data = (self.path / e["file"]).read_bytes()
        else:
            self.pack.seek(self.data_start + e["offset"])
            data = self.pack.read(e["bytes"])
        if hashlib.sha256(data).hexdigest() != e["sha256"]:
            raise ValueError(f"{locus}: checksum mismatch")
        return data

    def locus(self, locus):
        return decode_locus(self.blob(locus))


def seq_digest(seq):
    """Digest of an allele/protein sequence as stored in version-2 files."""
    return max(1, struct.unpack("<Q", hashlib.sha256(seq).digest()[:8])[0])


def summary(alleles, pairs):
    """Per-locus numbers, defined as in `cgdist-cache stats`."""
    lens = sorted(a[1] for a in alleles if a[1] > 0)
    snps = sorted(p[2] for p in pairs)
    n = len(pairs) or 1
    med = lambda v: v[len(v) // 2] if v else 0
    coding = [p[5] for p in pairs if p[5] is not None]
    out = {
        "alleles": len(alleles), "pairs": len(pairs),
        "len_min": lens[0] if lens else 0, "len_median": med(lens), "len_max": lens[-1] if lens else 0,
        "snps_mean": sum(snps) / n, "snps_median": med(snps), "snps_max": snps[-1] if snps else 0,
        "indel_pair_frac": sum(1 for p in pairs if p[3] > 0) / n,
        "indel_events_mean": sum(p[3] for p in pairs) / n,
        "indel_bases_mean": sum(p[4] for p in pairs) / n,
    }
    if coding:
        out.update(syn=sum(c[0] for c in coding), nonsyn=sum(c[1] for c in coding),
                   frame_disrupted=sum(c[2] for c in coding))
    return out


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__)
    cmd, st = argv[1], Store(argv[2])
    m = st.manifest
    if cmd == "info":
        kind = "protein" if "protein" in m else "dna"
        print(json.dumps({k: m.get(k) for k in ("format", "format_version", "hasher", "alignment",
                                                   "protein", "genetic_code", "schema", "cgdist_version")},
                         indent=2))
        print(f"{kind} store: {len(m['loci'])} loci, {sum(e['pairs'] for e in m['loci'].values())} pairs")
    elif cmd == "verify":
        bad = 0
        for locus, e in m["loci"].items():
            try:
                alleles, pairs = st.locus(locus)
                n = len(alleles)
                complete = len(pairs) == n * (n - 1) // 2
                if (len(alleles), len(pairs), complete) != (e["alleles"], e["pairs"], e["complete"]):
                    raise ValueError("counts differ from the manifest")
            except Exception as err:  # noqa: BLE001 - report every locus
                print(f"✗ {locus}: {err}")
                bad += 1
        print(f"{len(m['loci']) - bad} of {len(m['loci'])} loci verified")
        return 1 if bad else 0
    elif cmd == "dump":
        alleles, pairs = st.locus(argv[3])
        print("hash1\thash2\tsnps\tindel_events\tindel_bases\tsyn\tnonsyn\tframe_disrupted")
        for a, b, s, e, ib, c in pairs:
            print(f"{a}\t{b}\t{s}\t{e}\t{ib}\t" + ("\t".join(map(str, c)) if c else "\t\t"))
    elif cmd == "summary":
        out = {}
        for locus in m["loci"]:
            out[locus] = summary(*st.locus(locus))
        json.dump(out, sys.stdout)
    elif cmd == "digests":
        import zlib
        schema = Path(argv[3])
        checked = bad = 0
        for locus in m["loci"]:
            f = schema / f"{locus}.fasta"
            if not f.exists():
                continue
            seqs, cur = {}, []
            for line in f.read_text().splitlines() + [">"]:
                if line.startswith(">"):
                    if cur:
                        sq = "".join(cur).encode()
                        seqs[zlib.crc32(sq) & 0xFFFFFFFF] = sq
                    cur = []
                elif line.strip():
                    cur.append(line.strip())
            for h, _, d in st.locus(locus)[0]:
                if d and h in seqs:
                    checked += 1
                    if seq_digest(seqs[h]) != d:
                        bad += 1
                        print(f"✗ {locus}: allele {h} has another sequence than in the store")
        print(f"{checked} digests checked, {bad} mismatches")
        return 1 if bad else 0
    else:
        sys.exit(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
