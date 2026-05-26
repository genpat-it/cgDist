#!/usr/bin/env python3
"""
Feature tests for cgDist:
  1. cache-only mode (pre-compute alignments without producing a matrix)
  2. recombination-candidate flagging via the recombination_candidate_analyzer
     binary, reading an enriched cache — the workflow described in the paper
     (Supplementary §S6): cgdist --enrich-lengths  ->  analyzer.

Self-contained: locates the binaries via $CGDIST_BIN / ../target/release / PATH,
generates everything from the committed synthetic fixtures, and runs from any
working directory.
"""

import os
import sys
import shutil
import subprocess

import pandas as pd

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


def run_command(cmd, description, expect_success=True):
    """Run a shell command; report and return success boolean."""
    print(f"\n🧪 {description}")
    print(f"   $ {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        print(f"   ⏰ TIMEOUT")
        return False
    ok = (result.returncode == 0) == expect_success
    if not ok:
        print(f"   ❌ unexpected exit {result.returncode}")
        print(f"   STDERR: {result.stderr.strip()[:300]}")
    return ok


def test_cache_only_mode():
    print("\n" + "=" * 70)
    print("CACHE-ONLY MODE")
    print("=" * 70)
    for f in ["results/test_cache_only_new.lz4", "results/test_from_cache_new.tsv"]:
        if os.path.exists(f):
            os.remove(f)

    passed = 0
    total = 3

    cmd = (f"{CGDIST} --schema {SCHEMA} --profiles {PROFILES} "
           "--mode snps-indel-bases --hasher-type crc32 "
           "--cache-file results/test_cache_only_new.lz4 --cache-only "
           "--cache-note 'cache-only test'")
    if run_command(cmd, "Cache-only execution") and os.path.exists("results/test_cache_only_new.lz4"):
        print("   ✅ cache created")
        passed += 1

    cmd = (f"{CGDIST} --schema {SCHEMA} --profiles {PROFILES} "
           "--output results/test_from_cache_new.tsv "
           "--mode snps-indel-bases --hasher-type crc32 "
           "--cache-file results/test_cache_only_new.lz4")
    if run_command(cmd, "Matrix from cached data") and os.path.exists("results/test_from_cache_new.tsv"):
        print("   ✅ matrix built from cache")
        passed += 1

    cmd = (f"{CGDIST} --schema {SCHEMA} --profiles {PROFILES} "
           "--mode snps-indel-bases --hasher-type crc32 --cache-only")
    if run_command(cmd, "Cache-only without --cache-file must fail", expect_success=False):
        print("   ✅ correctly rejected (cache-only requires --cache-file)")
        passed += 1

    return passed, total


def test_recombination_candidate_flagging():
    print("\n" + "=" * 70)
    print("RECOMBINATION-CANDIDATE FLAGGING (enriched cache -> analyzer)")
    print("=" * 70)

    cache = "results/recomb_cache.lz4"
    enriched = "results/recomb_enriched.bin"
    dist = "results/recomb_distance.tsv"
    corrected = "results/recomb_corrected.tsv"
    log = "results/candidate_recombination_loci.tsv"
    for f in [cache, enriched, dist, corrected, log]:
        if os.path.exists(f):
            os.remove(f)

    passed = 0
    total = 4

    # Step 1: build cache + distance matrix
    step1 = (f"{CGDIST} --schema {SCHEMA} --profiles {PROFILES} --output {dist} "
             f"--mode snps-indel-bases --hasher-type crc32 --cache-file {cache} --force-recompute")
    # Step 2: reload the cache and enrich it with sequence lengths (-> .bin)
    step2 = (f"{CGDIST} --schema {SCHEMA} --profiles {PROFILES} --output {dist} "
             f"--mode snps-indel-bases --hasher-type crc32 --cache-file {cache} "
             f"--enrich-lengths --enrich-output {enriched}")
    if (run_command(step1, "Build cache + distance matrix")
            and run_command(step2, "Enrich cache with sequence lengths")
            and os.path.exists(enriched)):
        print("   ✅ enriched cache created")
        passed += 1

    # Step 3: run the analyzer (canonical args)
    step3 = (f"{ANALYZER} --enriched-cache {enriched} --profiles {PROFILES} "
             f"--distance-matrix {dist} --output-matrix {corrected} "
             f"--candidate-recombination-log {log} --threshold 3.0")
    if run_command(step3, "Candidate flagging (threshold 3%)") and os.path.exists(log):
        print("   ✅ candidate flagging log created")
        passed += 1

    # Step 4: validate the log columns and values
    if os.path.exists(log):
        try:
            df = pd.read_csv(log, sep="\t")
            expected = [
                "sample_i", "sample_j", "locus", "cgdist_snps", "cgdist_indel_events",
                "cgdist_indel_bases", "cgdist_total", "hamming_distance",
                "avg_seq_length", "mutation_density_%", "recombination_excess_%",
            ]
            if all(c in df.columns for c in expected):
                print("   ✅ log columns correct")
                passed += 1
            else:
                print(f"   ❌ columns mismatch. found: {list(df.columns)}")
            if len(df) > 0 and df["mutation_density_%"].between(0, 100).all():
                print(f"   ✅ {len(df)} candidate loci flagged, densities in [0,100]")
                passed += 1
            else:
                print("   ❌ no candidates or out-of-range densities")
        except Exception as e:  # noqa: BLE001
            print(f"   ❌ error reading log: {e}")

    return passed, total


def main():
    print("🧪 cgDist cache-only & recombination-candidate feature tests")
    print("=" * 70)
    os.chdir(SCRIPT_DIR)
    missing = [n for n, p in [("cgdist", CGDIST),
                              ("recombination_candidate_analyzer", ANALYZER)] if p is None]
    if missing:
        print(f"ERROR: missing binary/binaries: {', '.join(missing)}")
        print("Build them (cd .. && cargo build --release) or 'cargo install cgdist'.")
        sys.exit(2)
    print(f"cgdist:   {CGDIST}")
    print(f"analyzer: {ANALYZER}")
    os.makedirs("results", exist_ok=True)

    cache_passed, cache_total = test_cache_only_mode()
    recomb_passed, recomb_total = test_recombination_candidate_flagging()

    total_passed = cache_passed + recomb_passed
    total = cache_total + recomb_total

    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    print(f"Cache-only mode:                  {cache_passed}/{cache_total}")
    print(f"Recombination-candidate flagging: {recomb_passed}/{recomb_total}")
    print(f"\nOverall: {total_passed}/{total}")
    if total_passed == total:
        print("🎉 ALL FEATURE TESTS PASSED!")
        return True
    print(f"⚠️  {total - total_passed} test(s) failed")
    return False


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
