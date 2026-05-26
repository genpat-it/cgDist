#!/usr/bin/env python3
"""
Functional tests for cgDist's data-handling, filtering, and output-format
arguments — verifying *behaviour*, not just that the flag is accepted:

  --min-loci, --missing-char (default and custom), --sample-threshold,
  --locus-threshold, --include-samples / --exclude-samples,
  --include-loci / --exclude-loci, and --format (tsv/csv/phylip/nexus).

Self-contained: locates the cgdist binary, builds tiny profiles on the fly
from the committed schema_crc32 alleles, and asserts the expected outcome.
Uses --mode hamming (no sequence alignment needed) so the focus stays on the
filtering/missing logic.
"""

import os
import re
import sys
import shutil
import subprocess
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCHEMA = "schema_crc32"

# Known-valid allele hashes from schema_crc32 (Sample_Ref alleles).
H1, H2, H3 = "334687599", "2310869957", "474191502"

results = []  # (name, ok, detail)


def find_cgdist():
    candidates = []
    env_bin = os.environ.get("CGDIST_BIN")
    if env_bin:
        candidates.append(env_bin)
    candidates.append(os.path.join(SCRIPT_DIR, "..", "target", "release", "cgdist"))
    on_path = shutil.which("cgdist")
    if on_path:
        candidates.append(on_path)
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    return None


CGDIST = find_cgdist()


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))


def write_profile(path, rows, missing="-"):
    """rows: list of (sample, a1, a2, a3) where '' means missing."""
    with open(path, "w") as f:
        f.write("sample\tlocus1\tlocus2\tlocus3\n")
        for s, a1, a2, a3 in rows:
            cells = [x if x != "" else missing for x in (a1, a2, a3)]
            f.write(s + "\t" + "\t".join(cells) + "\n")


def run(profile, out, extra, sep="\t"):
    """Run cgdist (hamming mode) and return (stdout, returncode, matrix dict, samples)."""
    cmd = [CGDIST, "--schema", SCHEMA, "--profiles", profile, "--output", out,
           "--mode", "hamming", "--hasher-type", "crc32"] + extra
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    matrix, samples = None, None
    if os.path.exists(out):
        lines = [l.rstrip("\n") for l in open(out) if not l.startswith("#") and l.strip()]
        if lines:
            samples = lines[0].split(sep)[1:]
            matrix = {}
            for row in lines[1:]:
                cells = row.split(sep)
                for s2, v in zip(samples, cells[1:]):
                    matrix[(cells[0], s2)] = v
    return p.stdout + p.stderr, p.returncode, matrix, samples


def after_counts(stdout):
    """Extract post-filter sample/loci counts from the 'A → B' filter summary."""
    s = re.search(r"Samples:\s*\d+\s*\S+\s*(\d+)", stdout)
    l = re.search(r"Loci:\s*\d+\s*\S+\s*(\d+)", stdout)
    return (int(s.group(1)) if s else None, int(l.group(1)) if l else None)


def main():
    os.chdir(SCRIPT_DIR)
    if CGDIST is None:
        print("ERROR: cgdist binary not found "
              "($CGDIST_BIN / ../target/release / PATH).")
        sys.exit(2)
    print(f"Using cgdist binary: {CGDIST}\n")

    d = tempfile.mkdtemp(prefix="cgdist_filt_")
    try:
        # ---- missing data (default '-') + --min-loci ----
        # S_full shares 2 loci with S_miss1, only 1 (locus3) with S_miss2.
        miss = os.path.join(d, "miss.tsv")
        write_profile(miss, [
            ("S_full",  H1, H2, H3),
            ("S_miss1", H1, "", H3),   # locus2 missing
            ("S_miss2", "", "", H3),   # only locus3 present
        ])

        _, rc, m, _ = run(miss, os.path.join(d, "m0.tsv"), [])
        check("missing '-' recognised; no min-loci -> no NA",
              rc == 0 and m is not None and "NA" not in m.values())

        _, rc, m, _ = run(miss, os.path.join(d, "m2.tsv"), ["--min-loci", "2"])
        ok = (rc == 0 and m is not None
              and m.get(("S_full", "S_miss2")) == "NA"     # 1 shared locus -> dropped
              and m.get(("S_miss1", "S_miss2")) == "NA"    # 1 shared locus -> dropped
              and m.get(("S_full", "S_miss1")) not in (None, "NA"))  # 2 shared -> kept
        check("--min-loci 2 -> under-covered pairs become NA", ok,
              "" if ok else f"matrix={m}")

        # ---- custom --missing-char '?' reproduces the same behaviour ----
        miss_q = os.path.join(d, "miss_q.tsv")
        write_profile(miss_q, [
            ("S_full",  H1, H2, H3),
            ("S_miss2", "", "", H3),
        ], missing="?")
        _, rc, m, _ = run(miss_q, os.path.join(d, "mq.tsv"),
                          ["--missing-char", "?", "--min-loci", "2"])
        check("--missing-char '?' honoured (same NA outcome)",
              rc == 0 and m is not None and m.get(("S_full", "S_miss2")) == "NA")

        # ---- --sample-threshold drops a mostly-missing sample ----
        qual = os.path.join(d, "qual.tsv")
        write_profile(qual, [
            ("S_full", H1, H2, H3),   # 100% complete
            ("S_ok",   H1, H2, ""),   # 67%
            ("S_bad",  H1, "", ""),   # 33%
        ])
        _, rc, m, samples = run(qual, os.path.join(d, "st.tsv"),
                                ["--sample-threshold", "0.5"])
        check("--sample-threshold 0.5 drops the 33%-complete sample",
              rc == 0 and samples is not None
              and "S_bad" not in samples and "S_full" in samples,
              "" if samples is None else f"samples={samples}")

        # ---- --locus-threshold drops a mostly-missing locus ----
        lq = os.path.join(d, "lq.tsv")
        write_profile(lq, [
            ("S1", H1, H2, H3),
            ("S2", H1, "", H3),   # locus2 missing
            ("S3", H1, "", H3),   # locus2 missing -> locus2 present in 1/3
        ])
        out, rc, _, _ = run(lq, os.path.join(d, "lt.tsv"),
                            ["--locus-threshold", "0.5"])
        _, loci_after = after_counts(out)
        check("--locus-threshold 0.5 drops the 33%-present locus (3 -> 2)",
              rc == 0 and loci_after == 2, f"loci_after={loci_after}")

        # ---- include/exclude samples (regex) ----
        sel = os.path.join(d, "sel.tsv")
        write_profile(sel, [
            ("Keep_A", H1, H2, H3),
            ("Keep_B", H2, H2, H3),
            ("Drop_C", H3, H2, H3),
        ])
        _, rc, _, samples = run(sel, os.path.join(d, "inc.tsv"),
                                ["--include-samples", "Keep_.*"])
        check("--include-samples regex keeps only matches",
              rc == 0 and samples is not None
              and set(samples) == {"Keep_A", "Keep_B"}, f"samples={samples}")

        _, rc, _, samples = run(sel, os.path.join(d, "exc.tsv"),
                                ["--exclude-samples", "Drop_.*"])
        check("--exclude-samples regex removes matches",
              rc == 0 and samples is not None and "Drop_C" not in samples,
              f"samples={samples}")

        # ---- include/exclude loci ----
        out, rc, _, _ = run(sel, os.path.join(d, "il.tsv"),
                            ["--include-loci", "locus1"])
        _, loci_after = after_counts(out)
        check("--include-loci locus1 keeps a single locus (3 -> 1)",
              rc == 0 and loci_after == 1, f"loci_after={loci_after}")

        out, rc, _, _ = run(sel, os.path.join(d, "el.tsv"),
                            ["--exclude-loci", "locus[23]"])
        _, loci_after = after_counts(out)
        check("--exclude-loci locus[23] leaves one locus (3 -> 1)",
              rc == 0 and loci_after == 1, f"loci_after={loci_after}")

        # ---- output formats ----
        # CSV: data rows are comma-separated.
        csvp = os.path.join(d, "f.csv")
        p = subprocess.run([CGDIST, "--schema", SCHEMA, "--profiles", sel,
                            "--output", csvp, "--mode", "hamming",
                            "--hasher-type", "crc32", "--format", "csv"],
                           capture_output=True, text=True)
        csv_ok = False
        if p.returncode == 0 and os.path.exists(csvp):
            data = [l for l in open(csvp) if not l.startswith("#") and l.strip()]
            csv_ok = len(data) >= 2 and "," in data[1] and "\t" not in data[1]
        check("--format csv produces comma-separated rows", csv_ok)

        # PHYLIP: first data line is the sample count.
        phyp = os.path.join(d, "f.phy")
        p = subprocess.run([CGDIST, "--schema", SCHEMA, "--profiles", sel,
                            "--output", phyp, "--mode", "hamming",
                            "--hasher-type", "crc32", "--format", "phylip"],
                           capture_output=True, text=True)
        phy_ok = False
        if p.returncode == 0 and os.path.exists(phyp):
            data = [l.strip() for l in open(phyp) if not l.startswith("#") and l.strip()]
            phy_ok = bool(data) and data[0].isdigit() and int(data[0]) == 3
        check("--format phylip starts with the sample count", phy_ok)

        # NEXUS: contains the #NEXUS header.
        nexp = os.path.join(d, "f.nex")
        p = subprocess.run([CGDIST, "--schema", SCHEMA, "--profiles", sel,
                            "--output", nexp, "--mode", "hamming",
                            "--hasher-type", "crc32", "--format", "nexus"],
                           capture_output=True, text=True)
        nex_ok = (p.returncode == 0 and os.path.exists(nexp)
                  and "#NEXUS" in open(nexp).read().upper())
        check("--format nexus emits a #NEXUS block", nex_ok)

        # ---- --hamming-fallback toggles SNPs-mode behaviour for InDel-only pairs ----
        # Uses the committed validation dataset (Sample_Dels_Only / Sample_Ins_Only
        # differ from Sample_Ref only by InDels; Sample_SNPs_Only by real SNPs).
        def snps_matrix(out, extra):
            cmd = [CGDIST, "--schema", SCHEMA, "--profiles",
                   "profiles/test_profiles_crc32.tsv", "--output", out,
                   "--mode", "snps", "--hasher-type", "crc32",
                   "--force-recompute"] + extra
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            m = {}
            if os.path.exists(out):
                rows = [l.rstrip("\n") for l in open(out)
                        if not l.startswith("#") and l.strip()]
                hdr = rows[0].split("\t")[1:]
                for row in rows[1:]:
                    c = row.split("\t")
                    for s2, v in zip(hdr, c[1:]):
                        m[(c[0], s2)] = v
            return p.returncode, m

        rc_off, off = snps_matrix(os.path.join(d, "snps_off.tsv"), [])
        rc_on, on = snps_matrix(os.path.join(d, "snps_on.tsv"), ["--hamming-fallback"])
        ref = "Sample_Ref"
        ok = (rc_off == 0 and rc_on == 0
              # without fallback, InDel-only pairs contribute 0 SNPs
              and off.get((ref, "Sample_Dels_Only")) == "0"
              and off.get((ref, "Sample_Ins_Only")) == "0"
              # with fallback, each differing locus contributes +1 (3 loci -> 3)
              and on.get((ref, "Sample_Dels_Only")) == "3"
              and on.get((ref, "Sample_Ins_Only")) == "3"
              # pairs with real SNPs are unaffected by the flag
              and off.get((ref, "Sample_SNPs_Only")) == "4"
              and on.get((ref, "Sample_SNPs_Only")) == "4")
        check("--hamming-fallback: InDel-only 0 (off) -> 3 (on); real SNPs stay 4", ok,
              f"off={off.get((ref,'Sample_Dels_Only'))} on={on.get((ref,'Sample_Dels_Only'))} "
              f"snps_off={off.get((ref,'Sample_SNPs_Only'))}")
    finally:
        shutil.rmtree(d, ignore_errors=True)

    # Report
    width = max(len(n) for n, _, _ in results)
    npass = sum(1 for _, ok, _ in results if ok)
    print("=" * (width + 16))
    for name, ok, detail in results:
        line = f"  [{'PASS' if ok else 'FAIL'}] {name.ljust(width)}"
        if not ok and detail:
            line += f"  <- {detail}"
        print(line)
    print("=" * (width + 16))
    print(f"{npass}/{len(results)} functional checks passed")
    if npass == len(results):
        print("ALL FILTER / MISSING / FORMAT CHECKS PASSED")
        return True
    return False


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
