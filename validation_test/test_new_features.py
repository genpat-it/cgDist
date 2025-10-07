#!/usr/bin/env python3
"""
Test script for new cgDist features: cache-only mode and recombination detection.
"""

import subprocess
import os
import sys
import pandas as pd

def run_command(cmd, description):
    """Run a command and check its exit status."""
    print(f"\n🧪 Testing: {description}")
    print(f"Command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            print(f"❌ FAILED: {description}")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return False
        else:
            print(f"✅ PASSED: {description}")
            return True
    except subprocess.TimeoutExpired:
        print(f"⏰ TIMEOUT: {description}")
        return False
    except Exception as e:
        print(f"❌ ERROR: {description} - {e}")
        return False

def test_cache_only_mode():
    """Test cache-only functionality."""
    print("\n" + "="*70)
    print("TESTING CACHE-ONLY MODE")
    print("="*70)
    
    # Clean up previous cache files
    cache_files = [
        "results/test_cache_only_new.lz4",
        "results/test_from_cache_new.tsv"
    ]
    for f in cache_files:
        if os.path.exists(f):
            os.remove(f)
    
    tests_passed = 0
    total_tests = 3
    
    # Test 1: Cache-only mode
    cmd = ("../target/release/cgdist --schema schema_crc32 "
           "--profiles profiles/test_profiles_crc32.tsv "
           "--mode snps-indel-bases --hasher-type crc32 "
           "--cache-file results/test_cache_only_new.lz4 --cache-only "
           "--cache-note 'Testing cache-only mode'")
    
    if run_command(cmd, "Cache-only mode execution"):
        if os.path.exists("results/test_cache_only_new.lz4"):
            print("✅ Cache file created successfully")
            tests_passed += 1
        else:
            print("❌ Cache file not found")
    
    # Test 2: Using cached data
    cmd = ("../target/release/cgdist --schema schema_crc32 "
           "--profiles profiles/test_profiles_crc32.tsv "
           "--output results/test_from_cache_new.tsv "
           "--mode snps-indel-bases --hasher-type crc32 "
           "--cache-file results/test_cache_only_new.lz4")
    
    if run_command(cmd, "Using cached data"):
        if os.path.exists("results/test_from_cache_new.tsv"):
            print("✅ Matrix created from cache successfully")
            tests_passed += 1
        else:
            print("❌ Output matrix not found")
    
    # Test 3: Validate cache error handling (cache-only without cache-file)
    cmd = ("../target/release/cgdist --schema schema_crc32 "
           "--profiles profiles/test_profiles_crc32.tsv "
           "--mode snps-indel-bases --hasher-type crc32 --cache-only")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0 and "requires --cache-file" in result.stderr:
        print("✅ Cache-only validation works correctly")
        tests_passed += 1
    else:
        print("❌ Cache-only validation failed")
    
    return tests_passed, total_tests

def test_recombination_detection():
    """Test recombination detection functionality."""
    print("\n" + "="*70)
    print("TESTING RECOMBINATION DETECTION")
    print("="*70)
    
    tests_passed = 0
    total_tests = 4
    
    # Clean up previous files
    log_files = [
        "results/recomb_low_threshold.csv",
        "results/recomb_high_threshold.csv",
        "results/test_recomb_matrix.tsv"
    ]
    for f in log_files:
        if os.path.exists(f):
            os.remove(f)
    
    # Test 1: Low threshold (should detect many events)
    cmd = ("../target/release/cgdist --schema schema_crc32 "
           "--profiles profiles/test_profiles_crc32.tsv "
           "--output results/test_recomb_matrix.tsv "
           "--mode snps-indel-bases --hasher-type crc32 "
           "--recombination-log results/recomb_low_threshold.csv "
           "--recombination-threshold 3")
    
    if run_command(cmd, "Recombination detection (low threshold)"):
        if os.path.exists("results/recomb_low_threshold.csv"):
            try:
                df = pd.read_csv("results/recomb_low_threshold.csv")
                if len(df) > 0:
                    print(f"✅ Detected {len(df)} recombination events (threshold=3)")
                    tests_passed += 1
                else:
                    print("❌ No events detected with low threshold")
            except Exception as e:
                print(f"❌ Error reading recombination log: {e}")
        else:
            print("❌ Recombination log not created")
    
    # Test 2: High threshold (should detect fewer/no events)  
    cmd = ("../target/release/cgdist --schema schema_crc32 "
           "--profiles profiles/test_profiles_crc32.tsv "
           "--output results/test_recomb_matrix.tsv "
           "--mode snps-indel-bases --hasher-type crc32 "
           "--recombination-log results/recomb_high_threshold.csv "
           "--recombination-threshold 50")
    
    if run_command(cmd, "Recombination detection (high threshold)"):
        if os.path.exists("results/recomb_high_threshold.csv"):
            try:
                df = pd.read_csv("results/recomb_high_threshold.csv")
                print(f"✅ Detected {len(df)} recombination events (threshold=50)")
                tests_passed += 1
            except Exception as e:
                print(f"❌ Error reading recombination log: {e}")
        else:
            print("❌ Recombination log not created")
    
    # Test 3: Validate CSV format
    if os.path.exists("results/recomb_low_threshold.csv"):
        try:
            df = pd.read_csv("results/recomb_low_threshold.csv")
            expected_columns = [
                'locus', 'sample1', 'sample2', 'allele1_hash', 'allele2_hash',
                'snps_indel_bases', 'threshold', 'seq_length1', 'seq_length2', 
                'divergence_percent'
            ]
            
            if all(col in df.columns for col in expected_columns):
                print("✅ CSV format is correct")
                tests_passed += 1
            else:
                print("❌ CSV format is incorrect")
                print(f"Expected: {expected_columns}")
                print(f"Found: {list(df.columns)}")
        except Exception as e:
            print(f"❌ Error validating CSV format: {e}")
    
    # Test 4: Validate divergence percentages are calculated
    if os.path.exists("results/recomb_low_threshold.csv"):
        try:
            df = pd.read_csv("results/recomb_low_threshold.csv")
            if len(df) > 0:
                divergence_values = df['divergence_percent']
                if all(0 <= val <= 100 for val in divergence_values):
                    print("✅ Divergence percentages are valid")
                    tests_passed += 1
                else:
                    print("❌ Invalid divergence percentages found")
            else:
                print("⚠️ No data to validate divergence percentages")
                tests_passed += 1  # Not a failure if no events detected
        except Exception as e:
            print(f"❌ Error validating divergence percentages: {e}")
    
    return tests_passed, total_tests

def main():
    """Main test runner."""
    print("🧪 cgDist New Features Validation Test")
    print("=" * 70)
    
    total_passed = 0
    total_tests = 0
    
    # Test cache-only mode
    cache_passed, cache_total = test_cache_only_mode()
    total_passed += cache_passed
    total_tests += cache_total
    
    # Test recombination detection
    recomb_passed, recomb_total = test_recombination_detection()
    total_passed += recomb_passed
    total_tests += recomb_total
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"Cache-only mode: {cache_passed}/{cache_total} tests passed")
    print(f"Recombination detection: {recomb_passed}/{recomb_total} tests passed")
    print(f"\nOverall: {total_passed}/{total_tests} tests passed")
    
    if total_passed == total_tests:
        print("🎉 ALL NEW FEATURE TESTS PASSED!")
        return True
    else:
        print(f"⚠️ {total_tests - total_passed} tests failed")
        return False

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)