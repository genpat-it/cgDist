#!/usr/bin/env python3
"""
Corrected validation script for cgDist with accurate expected values.
"""

import pandas as pd
import sys
import os
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def find_cgdist():
    """Locate the cgdist binary: $CGDIST_BIN, local release build, or PATH."""
    candidates = []
    env_bin = os.environ.get('CGDIST_BIN')
    if env_bin:
        candidates.append(env_bin)
    candidates.append(os.path.join(SCRIPT_DIR, '..', 'target', 'release', 'cgdist'))
    on_path = shutil.which('cgdist')
    if on_path:
        candidates.append(on_path)
    for cand in candidates:
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)
    return None


def generate_matrices(cgdist):
    """(Re)generate the four distance matrices that the validation loads.

    Running the binary here makes the suite a single command: it validates
    whatever cgdist is installed, whether built locally (../target/release)
    or installed via `cargo install cgdist` (on PATH).
    """
    schema = 'schema_crc32'
    profiles = 'profiles/test_profiles_crc32.tsv'
    os.makedirs('results', exist_ok=True)
    runs = [
        ('results/crc32_hamming.tsv', ['--mode', 'hamming']),
        ('results/crc32_snps.tsv', ['--mode', 'snps', '--hamming-fallback']),
        ('results/crc32_snps_indel_contiguous.tsv', ['--mode', 'snps-indel-contiguous']),
        ('results/crc32_snps_indel_bases.tsv', ['--mode', 'snps-indel-bases']),
    ]
    for out, mode_args in runs:
        cmd = [cgdist, '--schema', schema, '--profiles', profiles,
               '--output', out, '--hasher-type', 'crc32'] + mode_args
        print('  $ ' + ' '.join(cmd))
        subprocess.run(cmd, check=True)

def load_distance_matrix(filepath):
    """Load a distance matrix from TSV file."""
    with open(filepath, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]
    
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t', index_col=0)
    return df

def manual_sequence_analysis():
    """Manually analyze sequences to determine correct expected values."""
    
    sequences = {
        'locus1': {
            'ref': 'ATCGATCGATCGATCG',        # Reference
            'snp1': 'ATCGATCGATCGATGG',       # 1 SNP: C→G at position 15 
            'del1': 'ATCGATCGATCGTCG',        # 1 deletion (A missing at pos 11)
            'del3': 'ATCGATCGTCG',            # 5 deletions (GATCG missing)  
            'ins1': 'ATCGATCGATCGATCGA',      # 1 insertion (A at end)
            'ins3': 'ATCGATCGATCGATCGAAA',    # 3 insertions (AAA at end)
        },
        'locus2': {
            'ref': 'GGGGGGGGGGGGGGGG',        # Reference (16 G's)
            'snp1': 'GGGGGGGGGTGGGGGG',       # 1 SNP: G→T at position 9
            'del2': 'GGGGGGGGGGGGGG',         # 2 deletions at end
            'del5': 'GGGGGGGGGGG',            # 5 deletions at end
            'ins2': 'GGGGGGGGGGGGGGGGCC',     # 2 insertions (CC at end)
        },
        'locus3': {
            'ref': 'CCCCCCCCCCCCCCCC',        # Reference (16 C's)
            'snp3': 'CCCCCTTCCCCCCCCC',       # 2 SNPs: CC→TT at positions 5-6
            'del1': 'CCCCCCCCCCCCCCC',        # 1 deletion at end
            'ins4': 'CCCCCCCCCCCCCCCCTTTT',   # 4 insertions (TTTT at end)
            'complex1': 'CCCAAACCCCCCCCC',     # 3 SNPs (CCC→AAA) + 1 deletion
            'large_del': 'CCCCCCCC',          # 8 deletions (half missing)
        }
    }
    
    # Map sample profiles to sequences (based on create_validation_test.py)
    sample_profiles = {
        'Sample_Ref': ('ref', 'ref', 'ref'),
        'Sample_Identical': ('ref', 'ref', 'ref'),  # Note: identical alleles have same CRC32
        'Sample_SNPs_Only': ('snp1', 'snp1', 'snp3'),
        'Sample_Dels_Only': ('del1', 'del2', 'del1'),
        'Sample_Ins_Only': ('ins1', 'ins2', 'ins4'),
        'Sample_Mixed1': ('snp2', 'del5', 'complex1'),  # Need to check snp2
    }
    
    # Calculate expected distances
    expected_distances = {}
    
    # Sample_Ref vs Sample_SNPs_Only
    # locus1: ref vs snp1 = 1 SNP
    # locus2: ref vs snp1 = 1 SNP  
    # locus3: ref vs snp3 = 2 SNPs
    # Total: 4 SNPs
    expected_distances[('Sample_Ref', 'Sample_SNPs_Only')] = {
        'hamming': 3,
        'snps': 4,
        'snps_indel_contiguous': 4,
        'snps_indel_bases': 4
    }
    
    # Sample_Ref vs Sample_Dels_Only  
    # locus1: ref vs del1 = 0 SNPs, 1 del event, 1 del base
    # locus2: ref vs del2 = 0 SNPs, 1 del event, 2 del bases
    # locus3: ref vs del1 = 0 SNPs, 1 del event, 1 del base
    # Total: 0 SNPs (→3 Hamming fallback), 3 del events, 4 del bases
    expected_distances[('Sample_Ref', 'Sample_Dels_Only')] = {
        'hamming': 3,
        'snps': 3,  # Hamming fallback
        'snps_indel_contiguous': 3,  # 0 SNPs + 3 InDel events
        'snps_indel_bases': 4   # 0 SNPs + 4 InDel bases
    }
    
    # Sample_Ref vs Sample_Ins_Only
    # locus1: ref vs ins1 = 0 SNPs, 1 ins event, 1 ins base
    # locus2: ref vs ins2 = 0 SNPs, 1 ins event, 2 ins bases  
    # locus3: ref vs ins4 = 0 SNPs, 1 ins event, 4 ins bases
    # Total: 0 SNPs (→3 Hamming fallback), 3 ins events, 7 ins bases
    expected_distances[('Sample_Ref', 'Sample_Ins_Only')] = {
        'hamming': 3,
        'snps': 3,  # Hamming fallback
        'snps_indel_contiguous': 3,  # 0 SNPs + 3 InDel events
        'snps_indel_bases': 7   # 0 SNPs + 7 InDel bases
    }
    
    return expected_distances

def validate_cgdist():
    """Main validation function."""
    
    # Load matrices
    matrices = {
        'hamming': load_distance_matrix('results/crc32_hamming.tsv'),
        'snps': load_distance_matrix('results/crc32_snps.tsv'),
        'snps_indel_contiguous': load_distance_matrix('results/crc32_snps_indel_contiguous.tsv'),
        'snps_indel_bases': load_distance_matrix('results/crc32_snps_indel_bases.tsv')
    }
    
    # Get corrected expected values
    expected_distances = manual_sequence_analysis()
    
    print("=" * 80)
    print("CORRECTED cgDist VALIDATION RESULTS")
    print("=" * 80)
    
    all_passed = True
    
    # Test identical samples (must be 0)
    test_pair = ('Sample_Ref', 'Sample_Identical')
    print(f"\n📊 Testing: {test_pair[0]} vs {test_pair[1]}")
    print("   Expected: All distances = 0 (identical sequences)")
    print("-" * 60)
    
    for mode, matrix in matrices.items():
        actual = int(matrix.loc[test_pair[0], test_pair[1]])
        if actual == 0:
            status = "✅ PASS"
        else:
            status = "❌ FAIL"
            all_passed = False
        print(f"   {mode:20s}: Expected=0, Actual={actual} {status}")
    
    # Test expected distance pairs
    for pair, expected in expected_distances.items():
        print(f"\n📊 Testing: {pair[0]} vs {pair[1]}")
        print("-" * 60)
        
        test_passed = True
        for mode, exp_val in expected.items():
            actual = int(matrices[mode].loc[pair[0], pair[1]])
            
            if actual == exp_val:
                status = "✅ PASS"
            else:
                status = "❌ FAIL"
                test_passed = False
                all_passed = False
            
            print(f"   {mode:20s}: Expected={exp_val:3d}, Actual={actual:3d} {status}")
        
        if test_passed:
            print("   Overall: ✅ TEST PASSED")
        else:
            print("   Overall: ❌ TEST FAILED")
    
    # Mathematical invariant check (cgDist ≥ Hamming)
    print("\n" + "=" * 80)
    print("MATHEMATICAL INVARIANT CHECK (cgDist ≥ Hamming):")
    print("=" * 80)
    
    invariant_violations = 0
    for sample1 in matrices['hamming'].index:
        for sample2 in matrices['hamming'].columns:
            if sample1 != sample2:  # Skip diagonal
                hamming_dist = int(matrices['hamming'].loc[sample1, sample2])
                snps_dist = int(matrices['snps'].loc[sample1, sample2])
                
                if snps_dist < hamming_dist:
                    print(f"❌ VIOLATED: {sample1} vs {sample2}: SNPs({snps_dist}) < Hamming({hamming_dist})")
                    invariant_violations += 1
    
    if invariant_violations == 0:
        print("✅ Mathematical invariant maintained for all pairs")
    else:
        print(f"❌ Found {invariant_violations} invariant violations")
        all_passed = False
    
    # Summary of all distance patterns
    print("\n" + "=" * 80)
    print("DISTANCE MATRIX SUMMARY:")
    print("=" * 80)
    
    print("\\nDistance from Sample_Ref to all samples:")
    for mode, matrix in matrices.items():
        print(f"\\n{mode.upper()}:")
        ref_distances = matrix.loc['Sample_Ref'].sort_values()
        for sample, dist in ref_distances.items():
            if sample != 'Sample_Ref':
                print(f"  {sample:18s}: {int(dist):3d}")
    
    # Final result
    print("\n" + "=" * 80)
    if all_passed:
        print("🎉 ALL VALIDATION TESTS PASSED!")
        print("✅ cgDist is correctly calculating SNPs, InDel events, and InDel bases")
        print("✅ Mathematical invariants are maintained")
        print("✅ Parasail alignment integration is working correctly")
    else:
        print("⚠️  SOME TESTS FAILED")
    print("=" * 80)
    
    return all_passed

if __name__ == '__main__':
    # Run from this script's directory so the relative paths below (schema_crc32,
    # profiles/, results/) resolve regardless of the caller's working directory.
    os.chdir(SCRIPT_DIR)

    cgdist = find_cgdist()
    if cgdist is None:
        print("ERROR: could not find the 'cgdist' binary.")
        print("  Fix it with ONE of:")
        print("    - build locally:   (cd .. && cargo build --release)")
        print("    - install it:      cargo install cgdist   (puts it on your PATH)")
        print("    - point to it:     CGDIST_BIN=/path/to/cgdist python3 run_validation.py")
        sys.exit(2)

    print(f"Using cgdist binary: {cgdist}")
    print("Generating distance matrices...")
    generate_matrices(cgdist)
    print()

    success = validate_cgdist()
    sys.exit(0 if success else 1)