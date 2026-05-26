#!/usr/bin/env python3
"""
Comprehensive CLI smoke-test for cgDist.

Exercises every cgdist argument / combination on the committed synthetic
fixtures (schema_crc32 + profiles/test_profiles_crc32.tsv) and the
recombination_candidate_analyzer workflow, asserting each invocation runs
successfully (and produces a non-empty output where one is expected).

This is a *smoke* test: it checks that every argument still works after code
changes — it does NOT check numerical correctness (run_validation.py does).

Self-contained: locates binaries via $CGDIST_BIN / ../target/release / PATH,
writes all outputs to a temp directory, and runs from any working directory.
"""

import os
import sys
import shutil
import tempfile
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCHEMA = "schema_crc32"
PROFILES = "profiles/test_profiles_crc32.tsv"


def find_binary(name):
    """Locate a cgDist binary: $CGDIST_BIN (or its dir), ../target/release, PATH."""
    candidates = []
    env_bin = os.environ.get("CGDIST_BIN")
    if env_bin:
        if os.path.basename(env_bin) == name:
            candidates.append(env_bin)
        candidates.append(os.path.join(os.path.dirname(os.path.abspath(env_bin)), name))
    candidates.append(os.path.join(SCRIPT_DIR, "..", "target", "release", name))
    on_path = shutil.which(name)
    if on_path:
        candidates.append(on_path)
    for cand in candidates:
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)
    return None


CGDIST = find_binary("cgdist")
ANALYZER = find_binary("recombination_candidate_analyzer")
INSPECTOR = find_binary("inspector")

results = []  # (name, ok, detail)


def run_case(name, argv, expect_file=None, expect_rc=0, timeout=60):
    """Run one argv (list), record PASS/FAIL. Optionally require a non-empty file."""
    try:
        proc = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        results.append((name, False, "timeout"))
        return False
    except Exception as e:  # noqa: BLE001
        results.append((name, False, f"exception: {e}"))
        return False

    if proc.returncode != expect_rc:
        detail = (proc.stderr or proc.stdout or "").strip().splitlines()
        detail = detail[-1] if detail else f"rc={proc.returncode}"
        results.append((name, False, f"rc={proc.returncode}: {detail}"))
        return False

    if expect_file is not None:
        if not os.path.exists(expect_file) or os.path.getsize(expect_file) == 0:
            results.append((name, False, "expected output missing/empty"))
            return False

    results.append((name, True, "ok"))
    return True


def cg(out_dir, name, extra, output="out.tsv", check_output=True, **kw):
    """Helper: a standard cgdist run with --schema/--profiles + extra args."""
    out_path = os.path.join(out_dir, f"{name.replace(' ', '_').replace('/', '_')}_{output}")
    argv = [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES,
            "--output", out_path] + extra
    run_case(name, argv, expect_file=out_path if check_output else None, **kw)
    return out_path


def main():
    os.chdir(SCRIPT_DIR)
    missing = [n for n, p in [("cgdist", CGDIST),
                              ("recombination_candidate_analyzer", ANALYZER)] if p is None]
    if missing:
        print(f"ERROR: missing binary/binaries: {', '.join(missing)}")
        print("Build them (cd .. && cargo build --release) or 'cargo install cgdist'.")
        sys.exit(2)

    print("cgDist CLI smoke-test")
    print(f"  cgdist:   {CGDIST}")
    print(f"  analyzer: {ANALYZER}\n")

    d = tempfile.mkdtemp(prefix="cgdist_cli_")
    try:
        # --version
        run_case("--version", [CGDIST, "--version"])

        # Distance modes (incl. legacy aliases)
        for m in ["hamming", "snps", "snps-indel-contiguous", "snps-indel-bases",
                  "snps-indel-events", "snps+indel-events"]:
            cg(d, f"mode={m}", ["--mode", m, "--hasher-type", "crc32"])

        # Hamming fallback flags
        cg(d, "snps +hamming-fallback", ["--mode", "snps", "--hamming-fallback", "--hasher-type", "crc32"])
        cg(d, "snps +no-hamming-fallback (no-op)", ["--mode", "snps", "--no-hamming-fallback", "--hasher-type", "crc32"])

        # Output formats
        for fmt in ["tsv", "csv", "phylip", "nexus"]:
            cg(d, f"format={fmt}", ["--format", fmt, "--mode", "snps", "--hasher-type", "crc32"])

        # Hashers
        for h in ["crc32", "sha256", "md5", "sequence"]:
            cg(d, f"hasher={h}", ["--mode", "snps", "--hasher-type", h])
        # The hamming hasher is allelic-only: no --schema, default mode.
        hh_out = os.path.join(d, "hamming_hasher.tsv")
        run_case("hasher=hamming (allelic, no schema)",
                 [CGDIST, "--profiles", PROFILES, "--output", hh_out, "--hasher-type", "hamming"],
                 expect_file=hh_out)

        # Alignment modes
        for am in ["dna", "dna-strict", "dna-permissive"]:
            cg(d, f"alignment-mode={am}", ["--mode", "snps-indel-bases", "--alignment-mode", am, "--hasher-type", "crc32"])
        # Custom alignment scoring (enables custom mode)
        cg(d, "custom alignment scores",
           ["--mode", "snps-indel-bases", "--alignment-mode", "custom",
            "--match-score", "2", "--mismatch-penalty", "-1",
            "--gap-open", "6", "--gap-extend", "1", "--hasher-type", "crc32"])
        # Friendly validation: wrong-sign custom scores must error cleanly, not abort.
        run_case("custom scores wrong sign -> clean error",
                 [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES,
                  "--output", os.path.join(d, "bad.tsv"), "--alignment-mode", "custom",
                  "--match-score", "2", "--mismatch-penalty", "4", "--hasher-type", "crc32"],
                 expect_rc=1)

        # Filters (regex)
        cg(d, "include-samples regex", ["--include-samples", "Sample_S.*", "--hasher-type", "crc32"])
        cg(d, "exclude-samples regex", ["--exclude-samples", "Sample_Large.*", "--hasher-type", "crc32"])
        cg(d, "include-loci regex", ["--include-loci", "locus[12]", "--hasher-type", "crc32"])
        cg(d, "exclude-loci regex", ["--exclude-loci", "locus3", "--hasher-type", "crc32"])

        # Filters (list files)
        loci_list = os.path.join(d, "loci.txt")
        with open(loci_list, "w") as fh:
            fh.write("locus1\nlocus2\n")
        samples_list = os.path.join(d, "samples.txt")
        with open(samples_list, "w") as fh:
            fh.write("Sample_Ref\nSample_SNPs_Only\nSample_Dels_Only\n")
        cg(d, "include-loci-list", ["--include-loci-list", loci_list, "--hasher-type", "crc32"])
        cg(d, "exclude-loci-list", ["--exclude-loci-list", loci_list, "--hasher-type", "crc32"], check_output=False)
        cg(d, "include-samples-list", ["--include-samples-list", samples_list, "--hasher-type", "crc32"])
        cg(d, "exclude-samples-list", ["--exclude-samples-list", samples_list, "--hasher-type", "crc32"])

        # Quality / threshold filters
        cg(d, "min-loci", ["--min-loci", "1", "--hasher-type", "crc32"])
        cg(d, "sample-threshold", ["--sample-threshold", "0.5", "--hasher-type", "crc32"])
        cg(d, "locus-threshold", ["--locus-threshold", "0.5", "--hasher-type", "crc32"])
        cg(d, "missing-char", ["--missing-char", "?", "--hasher-type", "crc32"])

        # Threads
        for t in ["1", "2", "0"]:
            cg(d, f"threads={t}", ["--threads", t, "--mode", "snps", "--hasher-type", "crc32"])

        # save-alignments
        aln = os.path.join(d, "alignments.tsv")
        cg(d, "save-alignments", ["--mode", "snps-indel-bases", "--save-alignments", aln,
                                  "--hasher-type", "crc32", "--force-recompute"])
        run_case("save-alignments output exists", ["test", "-s", aln])

        # stats-only / dry-run (no --output)
        run_case("stats-only", [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES,
                                "--stats-only", "--hasher-type", "crc32"])
        run_case("dry-run", [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES,
                             "--output", os.path.join(d, "dry.tsv"),
                             "--dry-run", "--hasher-type", "crc32"])

        # generate-config + use it
        cfg = os.path.join(d, "cgdist-config.toml")
        run_case("generate-config", [CGDIST, "--generate-config", "--config", cfg])
        if os.path.exists(cfg):
            cg(d, "use --config", ["--config", cfg, "--hasher-type", "crc32"], check_output=False)

        # Cache lifecycle: create, cache-only, reuse, force-recompute, enrich, inspect
        cache = os.path.join(d, "cache.lz4")
        cg(d, "cache create", ["--mode", "snps-indel-bases", "--hasher-type", "crc32",
                               "--cache-file", cache, "--cache-note", "smoke", "--force-recompute"])
        cache_only = os.path.join(d, "cache_only.lz4")
        run_case("cache-only (no matrix)",
                 [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES, "--mode", "snps-indel-bases",
                  "--hasher-type", "crc32", "--cache-file", cache_only, "--cache-only"])
        run_case("cache-only produced cache", ["test", "-s", cache_only])
        cg(d, "cache reuse", ["--mode", "snps-indel-bases", "--hasher-type", "crc32", "--cache-file", cache])
        enriched = os.path.join(d, "enriched.bin")
        run_case("enrich-lengths + enrich-output",
                 [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES, "--output", os.path.join(d, "enr.tsv"),
                  "--mode", "snps-indel-bases", "--hasher-type", "crc32",
                  "--cache-file", cache, "--enrich-lengths", "--enrich-output", enriched])
        run_case("enriched cache produced", ["test", "-s", enriched])
        run_case("inspector (cgdist --inspector)",
                 [CGDIST, "--inspector", cache])
        # NOTE: the standalone `inspector` binary is currently stale (its loader
        # does not handle the v2 cache format that cgdist writes); use
        # `cgdist --inspector` instead. Tracked as a separate issue.

        # benchmark (short)
        run_case("benchmark (1s)", [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES,
                                    "--benchmark", "--benchmark-duration", "1", "--hasher-type", "crc32"],
                 timeout=30)

        # Analyzer workflow (needs the enriched cache + a distance matrix)
        dist = os.path.join(d, "analyzer_dist.tsv")
        run_case("matrix for analyzer",
                 [CGDIST, "--schema", SCHEMA, "--profiles", PROFILES, "--output", dist,
                  "--mode", "snps-indel-bases", "--hasher-type", "crc32", "--cache-file", cache],
                 expect_file=dist)
        corrected = os.path.join(d, "corrected.tsv")
        log = os.path.join(d, "candidate_recombination_loci.tsv")
        run_case("analyzer: candidate flagging",
                 [ANALYZER, "--enriched-cache", enriched, "--profiles", PROFILES,
                  "--distance-matrix", dist, "--output-matrix", corrected,
                  "--candidate-recombination-log", log, "--threshold", "3.0"],
                 expect_file=log)
        run_case("analyzer: legacy --recombination-log alias",
                 [ANALYZER, "--enriched-cache", enriched, "--profiles", PROFILES,
                  "--distance-matrix", dist, "--output-matrix", os.path.join(d, "corr2.tsv"),
                  "--recombination-log", os.path.join(d, "log2.tsv"), "--threshold", "3.0"],
                 expect_file=os.path.join(d, "log2.tsv"))
    finally:
        shutil.rmtree(d, ignore_errors=True)

    # Report
    width = max(len(n) for n, _, _ in results)
    print("=" * (width + 20))
    npass = sum(1 for _, ok, _ in results if ok)
    for name, ok, detail in results:
        mark = "PASS" if ok else "FAIL"
        line = f"  [{mark}] {name.ljust(width)}"
        if not ok:
            line += f"  <- {detail}"
        print(line)
    print("=" * (width + 20))
    print(f"{npass}/{len(results)} CLI cases passed")
    if npass == len(results):
        print("ALL CLI ARGUMENT CASES PASSED")
        return True
    print(f"{len(results) - npass} CASE(S) FAILED")
    return False


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
