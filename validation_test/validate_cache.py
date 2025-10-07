#!/usr/bin/env python3
"""
cgDist Cache Validator

Validates cache consistency by:
1. Loading cache with different distance modes
2. Computing distances with and without cache
3. Comparing results for consistency
4. Checking cache integrity and metadata
"""

import subprocess
import pandas as pd
import sys
import json
import os
from pathlib import Path

def run_cgdist(args):
    """Run cgDist and return stdout, stderr, returncode."""
    cmd = ['../target/release/cgdist'] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr, result.returncode

def load_distance_matrix(filepath):
    """Load a distance matrix from TSV file."""
    with open(filepath, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]
    
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t', index_col=0)
    return df

def validate_cache_consistency():
    """Test that cache produces consistent results across distance modes."""
    
    print("=" * 80)
    print("CACHE CONSISTENCY VALIDATION")
    print("=" * 80)
    
    cache_file = "results/validation_cache.lz4"
    if not os.path.exists(cache_file):
        print("❌ Cache file not found. Run with --cache-file first.")
        return False
    
    schema_dir = "schema_crc32"
    profiles = "profiles/test_profiles_crc32.tsv"
    
    distance_modes = [
        ("hamming", "Hamming distance"),
        ("snps", "SNPs only"), 
        ("snps-indel-events", "SNPs + InDel events"),
        ("snps-indel-bases", "SNPs + InDel bases")
    ]
    
    results = {}
    all_passed = True
    
    for mode, description in distance_modes:
        print(f"\n📊 Testing {description} mode...")
        
        # Test 1: With cache
        output_with_cache = f"results/cache_test_{mode}_with.tsv"
        stdout, stderr, rc = run_cgdist([
            '--schema', schema_dir,
            '--profiles', profiles, 
            '--output', output_with_cache,
            '--mode', mode,
            '--hasher-type', 'crc32',
            '--cache-file', cache_file
        ])
        
        if rc != 0:
            print(f"❌ Failed with cache: {stderr}")
            all_passed = False
            continue
            
        # Test 2: Without cache (force recompute)
        output_without_cache = f"results/cache_test_{mode}_without.tsv"
        stdout, stderr, rc = run_cgdist([
            '--schema', schema_dir,
            '--profiles', profiles,
            '--output', output_without_cache, 
            '--mode', mode,
            '--hasher-type', 'crc32',
            '--force-recompute'
        ])
        
        if rc != 0:
            print(f"❌ Failed without cache: {stderr}")
            all_passed = False
            continue
            
        # Compare results
        try:
            matrix_with = load_distance_matrix(output_with_cache)
            matrix_without = load_distance_matrix(output_without_cache)
            
            # Check if matrices are identical
            if matrix_with.equals(matrix_without):
                print(f"   ✅ PASS - Matrices identical with/without cache")
                results[mode] = True
            else:
                print(f"   ❌ FAIL - Matrices differ with/without cache")
                # Show differences
                diff_mask = matrix_with != matrix_without
                if diff_mask.any().any():
                    print("   Differences found:")
                    for i in matrix_with.index:
                        for j in matrix_with.columns:
                            if diff_mask.loc[i, j]:
                                print(f"     {i} vs {j}: with_cache={matrix_with.loc[i, j]}, without_cache={matrix_without.loc[i, j]}")
                results[mode] = False
                all_passed = False
                
        except Exception as e:
            print(f"   ❌ Error comparing matrices: {e}")
            results[mode] = False
            all_passed = False
            
        # Clean up test files
        try:
            os.remove(output_with_cache)
            os.remove(output_without_cache)
        except:
            pass
    
    return all_passed, results

def validate_cache_metadata():
    """Validate cache metadata and structure."""
    
    print("\n" + "=" * 80)
    print("CACHE METADATA VALIDATION")
    print("=" * 80)
    
    cache_file = "results/validation_cache.lz4"
    
    # Use cgdist inspector
    stdout, stderr, rc = run_cgdist(['--inspector', cache_file])
    
    if rc != 0:
        print(f"❌ Failed to inspect cache: {stderr}")
        return False
    
    print("📋 Cache inspection results:")
    print(stdout)
    
    # Check for expected content
    checks = [
        ("Version: 0.1.0", "✅ Version check"),
        ("Hasher type: crc32", "✅ Hasher type check"), 
        ("Distance mode: snps-indel-bases", "✅ Distance mode check"),
        ("Total entries: 70", "✅ Entry count check"),
        ("Unique loci: 3", "✅ Loci count check"),
        ("Validation test cache", "✅ User note check")
    ]
    
    all_passed = True
    for check_string, pass_msg in checks:
        if check_string in stdout:
            print(f"   {pass_msg}")
        else:
            print(f"   ❌ Missing: {check_string}")
            all_passed = False
    
    return all_passed

def validate_cache_performance():
    """Test cache performance benefits."""
    
    print("\n" + "=" * 80) 
    print("CACHE PERFORMANCE VALIDATION")
    print("=" * 80)
    
    cache_file = "results/validation_cache.lz4"
    schema_dir = "schema_crc32"
    profiles = "profiles/test_profiles_crc32.tsv"
    
    # Test with cache - should be fast
    print("🚀 Testing with cache (should be fast)...")
    stdout, stderr, rc = run_cgdist([
        '--schema', schema_dir,
        '--profiles', profiles,
        '--output', 'results/perf_test_cached.tsv',
        '--mode', 'snps-indel-bases',
        '--hasher-type', 'crc32', 
        '--cache-file', cache_file
    ])
    
    if rc != 0:
        print(f"❌ Failed with cache: {stderr}")
        return False
    
    # Look for cache hit information in output
    if "Already in cache: 70 (100.0%)" in stdout:
        print("   ✅ PASS - 100% cache hit rate achieved")
        cache_performance = True
    else:
        print("   ❌ FAIL - Expected 100% cache hit rate")
        cache_performance = False
    
    # Clean up
    try:
        os.remove('results/perf_test_cached.tsv')
    except:
        pass
        
    return cache_performance

def main():
    """Main validation function."""
    
    print("🔍 cgDist Cache Validation Suite")
    print("=" * 80)
    print("Validating cache consistency, metadata, and performance...")
    print()
    
    all_tests_passed = True
    
    # Test 1: Cache consistency across distance modes
    try:
        consistency_passed, mode_results = validate_cache_consistency()
        if not consistency_passed:
            all_tests_passed = False
    except Exception as e:
        print(f"❌ Cache consistency test failed: {e}")
        all_tests_passed = False
        
    # Test 2: Cache metadata validation
    try:
        metadata_passed = validate_cache_metadata()
        if not metadata_passed:
            all_tests_passed = False
    except Exception as e:
        print(f"❌ Cache metadata test failed: {e}")
        all_tests_passed = False
        
    # Test 3: Cache performance validation
    try:
        performance_passed = validate_cache_performance()
        if not performance_passed:
            all_tests_passed = False
    except Exception as e:
        print(f"❌ Cache performance test failed: {e}")
        all_tests_passed = False
    
    # Final summary
    print("\n" + "=" * 80)
    if all_tests_passed:
        print("🎉 ALL CACHE VALIDATION TESTS PASSED!")
        print("✅ Cache is consistent across all distance modes")
        print("✅ Cache metadata is correct and complete")
        print("✅ Cache provides expected performance benefits")
        print("✅ Scientific community can trust cache integrity")
    else:
        print("⚠️  SOME CACHE VALIDATION TESTS FAILED")
        print("Review the detailed output above for specific issues")
    print("=" * 80)
    
    return all_tests_passed

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)