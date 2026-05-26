# cgDist Validation Test Suite

This directory contains a comprehensive validation test suite for cgDist that demonstrates the algorithm's correctness in calculating SNPs, InDel events, and InDel bases from controlled sequences.

## Purpose

This validation suite serves to:
1. **Verify algorithmic correctness**: Ensure cgDist accurately counts mutations at the nucleotide level
2. **Validate Parasail integration**: Confirm sequence alignment produces expected results
3. **Check mathematical invariants**: Ensure cgDist ≥ Hamming distance relationship is maintained
4. **Provide scientific confidence**: Give the research community controlled test cases with known answers

## Test Design

### Controlled Schema
The test uses 3 loci with carefully designed allelic variants:

- **Locus 1**: Tests various SNP patterns and InDel types
- **Locus 2**: Tests deletions and insertions of different lengths  
- **Locus 3**: Tests complex patterns and large InDels

### Sample Profiles
10 test samples with specific mutation patterns:
- `Sample_Ref`: All reference alleles (baseline)
- `Sample_Identical`: Identical to reference (distance = 0)
- `Sample_SNPs_Only`: Only SNPs (4 total: 1+1+2 across loci)
- `Sample_Dels_Only`: Only deletions (4 total bases deleted)
- `Sample_Ins_Only`: Only insertions (7 total bases inserted)
- `Sample_Mixed*`: Complex combinations
- `Sample_Large_*`: Large InDels for stress testing

### Distance Modes Tested
- **Hamming**: Traditional allelic differences
- **SNPs**: Single nucleotide polymorphisms only
- **SNPs + InDel events**: SNPs plus number of InDel events
- **SNPs + InDel bases**: SNPs plus total InDel bases

## Key Validation Results ✅

### Test Cases
1. **Identical sequences**: Distance = 0 in all modes
2. **SNPs only**: Hamming=3, SNPs=4, Events=4, Bases=4
3. **Deletions only**: Hamming=3, SNPs=3 (fallback), Events=3, Bases=4
4. **Insertions only**: Hamming=3, SNPs=3 (fallback), Events=3, Bases=7
5. **Mathematical invariant**: cgDist ≥ Hamming for all pairs

### Key Insights
- **Hamming fallback**: When InDels exist but SNPs=0, algorithm correctly applies +1 per locus
- **Event vs Base counting**: Algorithm correctly distinguishes between number of InDel events and total bases affected
- **Parasail alignment**: Global alignment produces expected mutation counts
- **Cache efficiency**: Unified cache stores all statistics for rapid mode switching

## Files in this Directory

### Core Test Files
- `requirements.txt`: Python dependencies for the validation scripts (`pip install -r requirements.txt`)
- `setup_validation_test.py`: Generates CRC32-based schema and profiles
- `run_validation.py`: Main validation script with correct expected values (self-contained)
- `validate_cache.py`: Validates cache integrity and consistency (self-contained)
- `test_cache_only_and_recombination.py`: Tests cache-only mode and the recombination-candidate analyzer workflow
- `test_cli_arguments.py`: Comprehensive smoke-test exercising every cgdist argument/combination plus the analyzer
- `test_filters_and_missing.py`: Functional tests for the data-handling args — `--min-loci`, `--missing-char`, `--sample-threshold`, `--locus-threshold`, `--include`/`--exclude` samples & loci, and output formats (verifies behaviour, not just that the flag runs)
- `schema_crc32/`: FASTA files with controlled sequences
- `profiles/test_profiles_crc32.tsv`: Sample-to-allele mappings
- `results/`: Output distance matrices from cgDist

### Documentation
- `VALIDATION_SUMMARY.md`: Summary of validation results
- `README.md`: This documentation

## Running the Validation

### Quickest path: one command

`run_validation.py` is self-contained. It locates the `cgdist` binary
(checking `$CGDIST_BIN`, then `../target/release/cgdist`, then your
`PATH`), regenerates the four distance matrices, and validates them —
so after installing or building cgDist you only need:

```bash
pip install -r requirements.txt   # one-time: installs pandas
python3 run_validation.py
```

This works whether you installed via `cargo install cgdist` (binary on
your `PATH`) or built locally with `cargo build --release` (binary in
`../target/release/`). If the binary lives elsewhere, point to it with
`CGDIST_BIN=/path/to/cgdist python3 run_validation.py`. Expected final
line: `🎉 ALL VALIDATION TESTS PASSED!`.

The manual, step-by-step instructions below are equivalent and useful
if you want to inspect each intermediate matrix.

There are two ways to run the validation suite from a fresh clone:

- **Quick smoke test (recommended for first-time users)**: the input fixture
  `profiles/test_profiles_crc32.tsv` and the schema FASTA files are
  committed to the repository, so you can skip the setup step and go
  straight to running `cgdist`. Use this path if you just want to verify
  that the tool installs and runs correctly.
- **Regenerate from scratch**: run `setup_validation_test.py` to
  regenerate the input fixture, the FASTA schema, and the
  `EXPECTED_RESULTS_CRC32.md` documentation. Use this path if you want
  to verify that the test generator itself is reproducible, or if you
  modify the test scenarios. The regenerated fixture is byte-identical
  to the committed one (CRC32 hashes are deterministic).

### Step 1: Build cgDist
```bash
cd ..
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### Step 2: Generate Test Data (optional — fixture already committed)
```bash
cd validation_test
python3 setup_validation_test.py
```

### Step 3: Run cgDist Tests

> If you installed via `cargo install cgdist`, the binary is on your
> `PATH`: replace `../target/release/cgdist` with `cgdist` in the
> commands below.

```bash
# Test all distance modes
../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_hamming.tsv --mode hamming --hasher-type crc32

../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_snps.tsv --mode snps --hasher-type crc32 --hamming-fallback

../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_snps_indel_contiguous.tsv --mode snps-indel-contiguous --hasher-type crc32

../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/crc32_snps_indel_bases.tsv --mode snps-indel-bases --hasher-type crc32

# Test alignment saving with gaps
../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/test_alignments.tsv --mode snps-indel-bases --hasher-type crc32 --save-alignments results/alignments_with_gaps.tsv --force-recompute
```

### Step 4: Validate Results
```bash
python3 run_validation.py
```

### Step 5: Validate Cache Integrity
```bash
# Generate cache file
../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/cache_test.tsv --mode snps-indel-bases --hasher-type crc32 --cache-file results/validation_cache.lz4 --cache-note "Validation test cache" --force-recompute

# Validate cache consistency
python3 validate_cache.py
```

### Step 6: Test cache-only mode & recombination-candidate flagging
```bash
# Test cache-only mode and the enrich -> recombination_candidate_analyzer workflow
python3 test_cache_only_and_recombination.py
```

### Step 7: Comprehensive CLI argument smoke-test (optional)
```bash
# Exercise every cgdist argument/combination plus the analyzer
python3 test_cli_arguments.py
```

Expected output:
```
🎉 ALL VALIDATION TESTS PASSED!
✅ cgDist is correctly calculating SNPs, InDel events, and InDel bases
✅ Mathematical invariants are maintained
✅ Parasail alignment integration is working correctly

🎉 ALL CACHE VALIDATION TESTS PASSED!
✅ Cache is consistent across all distance modes
✅ Cache metadata is correct and complete
✅ Cache provides expected performance benefits

🎉 ALL NEW FEATURE TESTS PASSED!
✅ Cache-only mode working correctly
✅ Recombination detection functional
✅ CSV output format validated
```

## New Features (2025-09-02)

### Cache-Only Mode
- **Purpose**: Pre-compute alignments without generating distance matrix
- **Usage**: `--cache-only --cache-file cache.lz4` 
- **Benefit**: Separate alignment computation from distance matrix generation for large datasets

### Recombination Detection
- **Purpose**: Identify potential recombination events between alleles
- **Usage**: `--recombination-log events.csv --recombination-threshold 20`
- **Default threshold**: 20 SNPs+InDel bases (based on 3-4% divergence for typical MLST loci)
- **Output**: CSV log with locus, sample pairs, divergence percentages, and sequence lengths

## Scientific Significance

This validation demonstrates that:

1. **cgDist accurately counts mutations**: The algorithm correctly distinguishes between SNPs, InDel events, and InDel bases
2. **Parasail integration works correctly**: Global sequence alignment produces biologically meaningful results
3. **Mathematical properties are preserved**: The fundamental ordering relationship (cgDist ≥ Hamming) is maintained
4. **Cache architecture is sound**: Unified cache enables rapid switching between distance modes without recomputation
5. **Recombination detection**: Scientists can identify potentially recombined alleles that may skew phylogenetic analysis

## For the Scientific Community

This validation suite provides:
- **Reproducible test cases** with known ground truth
- **Transparent methodology** for verifying algorithmic correctness
- **Confidence in results** through controlled, testable scenarios
- **Foundation for further testing** with organism-specific data

The successful validation confirms that cgDist can be trusted for epidemiological analysis, outbreak investigation, and phylogenetic studies requiring nucleotide-level resolution.

---

**Note**: This validation uses synthetic sequences for controlled testing. For organism-specific validation, consider testing with known outbreak datasets where epidemiological relationships are well-established.