#!/usr/bin/env python3
"""Conformance test of the store format: the independent reference reader
(cgds_reader.py) must decode every locus of each store exactly as cgdist
does. For every locus, the Python summary is compared with the output of
`cgdist-cache stats` (Rust): counts, lengths, SNP/InDel statistics and
coding totals must be equal (means to 1e-12).

    test_store_format.py <path to cgdist-cache> <store dir | .cgpack> ...
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cgds_reader  # noqa: E402

INTS = ["alleles", "pairs", "len_min", "len_median", "len_max", "snps_median", "snps_max",
        "syn", "nonsyn", "frame_disrupted"]
FLOATS = ["snps_mean", "indel_pair_frac", "indel_events_mean", "indel_bases_mean"]


def check(exe, path):
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "stats.json"
        subprocess.run([exe, "stats", "--store", path, "--out", str(out), "--threads", "8"],
                       check=True, stdout=subprocess.DEVNULL)
        rust = {x["locus"]: x for x in json.loads(out.read_text())["loci"]}
    st = cgds_reader.Store(path)
    problems, pairs = 0, 0
    for locus in st.manifest["loci"]:
        py = cgds_reader.summary(*st.locus(locus))
        rs = rust[locus]
        pairs += py["pairs"]
        for k in INTS:
            if py.get(k) != rs.get(k):
                problems += 1
                print(f"  ✗ {locus} {k}: python {py.get(k)} rust {rs.get(k)}")
        for k in FLOATS:
            if abs(py[k] - rs[k]) > 1e-12 * max(1.0, abs(rs[k])):
                problems += 1
                print(f"  ✗ {locus} {k}: python {py[k]} rust {rs[k]}")
    kind = "protein" if "protein" in st.manifest else "dna"
    print(f"{'OK ' if not problems else 'FAIL'} {path}: {kind}, {len(rust)} loci, {pairs:,} pairs, {problems} differences")
    return problems


def main():
    exe, stores = sys.argv[1], sys.argv[2:]
    bad = sum(check(exe, s) for s in stores)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
