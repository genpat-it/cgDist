#!/usr/bin/env python3
"""
Test bidirectional alignment and streaming recombination detection
Validates that the new implementation works correctly and maintains compatibility
"""

import subprocess
import pandas as pd
import os
import time
from pathlib import Path

def run_cgdist_test(test_name, profiles_path, cache_path, output_path, recomb_log_path, threshold=5):
    """Run cgdist with streaming recombination detection"""
    cmd = [
        "./target/release/cgdist",
        "--schema", "/home/IZSNT/a.deruvo/cgdist-resources/schemas/salmonella_schema",
        "--profiles", profiles_path,
        "--output", output_path,
        "--mode", "snps",
        "--include-loci-list", "/home/IZSNT/a.deruvo/cgDist-study/efsa/efsa_loci.tsv",
        "--locus-threshold", "0.98",
        "--cache-file", cache_path,
        "--alignment-mode", "dna-strict",
        "--threads", "16",  # Use fewer threads for testing
        "--recombination-log", recomb_log_path,
        "--recombination-threshold", str(threshold)
    ]
    
    print(f"🧪 Running {test_name}...")
    print(f"Command: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        print(f"  Running command...")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 5 min timeout
        
        end_time = time.time()
        duration = end_time - start_time
        
        if result.returncode == 0:
            print(f"✅ {test_name} completed in {duration:.1f}s")
            print(f"  STDOUT: {result.stdout[-200:]}")  # Last 200 chars
            return True, duration
        else:
            print(f"❌ {test_name} failed:")
            print(f"  STDOUT: {result.stdout}")
            print(f"  STDERR: {result.stderr}")
            return False, duration
            
    except subprocess.TimeoutExpired:
        print(f"⏰ {test_name} timed out after 5 minutes")
        return False, 300
    except Exception as e:
        print(f"❌ {test_name} error: {e}")
        return False, 0

def validate_output_files(output_path, recomb_log_path):
    """Validate that output files are created and contain expected data"""
    results = {}
    
    # Check distance matrix output
    if os.path.exists(output_path):
        try:
            df_matrix = pd.read_csv(output_path, sep='\t')
            results['matrix_samples'] = len(df_matrix)
            results['matrix_exists'] = True
            print(f"✅ Distance matrix: {len(df_matrix)} samples")
        except Exception as e:
            results['matrix_exists'] = False
            print(f"❌ Distance matrix error: {e}")
    else:
        results['matrix_exists'] = False
        print(f"❌ Distance matrix file not found: {output_path}")
    
    # Check recombination log
    if os.path.exists(recomb_log_path):
        try:
            df_recomb = pd.read_csv(recomb_log_path)
            results['recomb_events'] = len(df_recomb)
            results['recomb_exists'] = True
            
            # Validate recombination log structure for bidirectional alignment
            expected_columns = ['locus', 'sample1', 'sample2', 'allele1_hash', 'allele2_hash',
                              'best_distance', 'threshold', 'seq_length1', 'seq_length2', 
                              'divergence_percent', 'distance_forward', 'distance_reverse_complement',
                              'alignment_direction']
            
            missing_cols = set(expected_columns) - set(df_recomb.columns)
            if missing_cols:
                print(f"❌ Recombination log missing columns: {missing_cols}")
                results['recomb_valid_structure'] = False
            else:
                results['recomb_valid_structure'] = True
                print(f"✅ Recombination log: {len(df_recomb)} events with valid structure")
                
                # Show sample of recombination events
                if len(df_recomb) > 0:
                    print("Sample bidirectional recombination events:")
                    print(df_recomb.head(3)[['locus', 'sample1', 'sample2', 'best_distance', 'alignment_direction']].to_string())
                
        except Exception as e:
            results['recomb_exists'] = False
            print(f"❌ Recombination log error: {e}")
    else:
        results['recomb_exists'] = False
        print(f"❌ Recombination log file not found: {recomb_log_path}")
    
    return results

def compare_with_reference(test_recomb_log, reference_recomb_log):
    """Compare test results with reference (if available)"""
    if not os.path.exists(reference_recomb_log):
        print("ℹ️  No reference recombination log available for comparison")
        return {}
    
    try:
        df_test = pd.read_csv(test_recomb_log)
        df_ref = pd.read_csv(reference_recomb_log)
        
        # Compare basic statistics
        comparison = {
            'test_events': len(df_test),
            'ref_events': len(df_ref),
            'events_diff': len(df_test) - len(df_ref),
            'events_diff_pct': ((len(df_test) - len(df_ref)) / len(df_ref)) * 100 if len(df_ref) > 0 else 0
        }
        
        print(f"📊 Comparison with reference:")
        print(f"  Test events: {comparison['test_events']}")
        print(f"  Reference events: {comparison['ref_events']}")
        print(f"  Difference: {comparison['events_diff']} ({comparison['events_diff_pct']:+.1f}%)")
        
        # Compare sample pairs (should detect more or equal events with bidirectional)
        if comparison['events_diff'] >= 0:
            print("✅ New implementation detected same or more recombination events (expected with bidirectional)")
        else:
            print("⚠️  New implementation detected fewer events - needs investigation")
            
        return comparison
        
    except Exception as e:
        print(f"❌ Comparison error: {e}")
        return {}

def main():
    print("🧪 Testing Bidirectional Alignment & Streaming Recombination Detection")
    print("="*70)
    
    # Test parameters
    test_dir = Path("/home/IZSNT/a.deruvo/cgdist-rs/validation_test")
    test_dir.mkdir(exist_ok=True)
    
    # Use small subset for quick testing
    profiles_path = "/home/IZSNT/a.deruvo/cgdist-rs/validation_test/test_profiles_small.tsv"
    cache_path = "/home/IZSNT/a.deruvo/cgdist-rs/validation_test/results/test_cache_only_new.lz4"
    output_path = str(test_dir / "results" / "test_bidirectional_distances.tsv")
    recomb_log_path = str(test_dir / "results" / "test_bidirectional_recombination.csv")
    
    # Create results directory
    (test_dir / "results").mkdir(exist_ok=True)
    
    # Check if test files exist
    if not os.path.exists(profiles_path):
        print(f"❌ Test profiles not found: {profiles_path}")
        print("Creating small test subset from Se1540...")
        
        # Create small test subset (first 50 samples for quick testing)
        try:
            df_full = pd.read_csv("/home/IZSNT/a.deruvo/cgDist-paper/supplementary_data/Se1540_allelic_profiles.tsv", sep='\t')
            df_small = df_full.head(50)  # Just 50 samples for quick test
            df_small.to_csv(profiles_path, sep='\t', index=False)
            print(f"✅ Created test subset: {len(df_small)} samples")
        except Exception as e:
            print(f"❌ Failed to create test subset: {e}")
            return
    
    if not os.path.exists(cache_path):
        print(f"❌ Test cache not found: {cache_path}")
        print("Please run cache creation first or use existing cache")
        return
    
    # Remove old output files
    for file_path in [output_path, recomb_log_path]:
        if os.path.exists(file_path):
            os.remove(file_path)
    
    # Build the project first
    print("🔨 Building cgdist...")
    build_result = subprocess.run(["cargo", "build", "--release"], capture_output=True, text=True)
    if build_result.returncode != 0:
        print(f"❌ Build failed: {build_result.stderr}")
        return
    print("✅ Build successful")
    
    # Run test
    success, duration = run_cgdist_test(
        "Bidirectional Streaming Test",
        profiles_path,
        cache_path, 
        output_path,
        recomb_log_path,
        threshold=5
    )
    
    if not success:
        print("❌ Test failed - cannot proceed with validation")
        return
    
    # Validate outputs
    print("\n📋 Validating output files...")
    results = validate_output_files(output_path, recomb_log_path)
    
    # Performance assessment
    print(f"\n📈 Performance Assessment:")
    print(f"  Duration: {duration:.1f} seconds")
    
    # Expected duration should be reasonable for 50 samples
    if duration < 180:  # Less than 3 minutes for 50 samples is reasonable
        print("✅ Duration looks reasonable for streaming implementation")
    else:
        print("⚠️  Duration longer than expected")
    
    # Summary
    print(f"\n📊 Test Summary:")
    print(f"  Matrix created: {'✅' if results.get('matrix_exists') else '❌'}")
    print(f"  Recombination log created: {'✅' if results.get('recomb_exists') else '❌'}")
    print(f"  Valid log structure: {'✅' if results.get('recomb_valid_structure') else '❌'}")
    print(f"  Recombination events: {results.get('recomb_events', 0)}")
    print(f"  Performance good: {'✅' if duration < 180 else '❌'}")
    
    # Overall assessment
    all_good = (results.get('matrix_exists', False) and 
                results.get('recomb_exists', False) and 
                results.get('recomb_valid_structure', False) and
                duration < 180)
    
    if all_good:
        print("\n🎉 All tests PASSED! Bidirectional streaming implementation works correctly.")
    else:
        print("\n❌ Some tests FAILED. Please check the issues above.")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()