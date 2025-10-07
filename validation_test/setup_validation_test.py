#!/usr/bin/env python3
"""
Create a controlled validation test for cgDist with known SNPs and InDels.
This script generates FASTA schema files and allelic profiles with CRC32 hashes.
"""

import zlib
import os

def crc32_hash(sequence):
    """Calculate CRC32 hash of a sequence."""
    return zlib.crc32(sequence.encode('utf-8')) & 0xffffffff

def create_test_cases():
    """Create controlled test cases with known differences."""
    
    # Define test sequences with known differences
    test_cases = {
        'locus1': [
            ('ref', 'ATCGATCGATCGATCG'),  # Reference
            ('snp1', 'ATCGATCGATCGATGG'),  # 2 SNPs at end
            ('snp2', 'GTCGATCGATCGATCG'),  # 1 SNP at start
            ('del1', 'ATCGATCGATCGTCG'),   # 1 base deletion (A missing)
            ('del3', 'ATCGATCGTCG'),       # 3 base deletion (GAT missing)
            ('ins1', 'ATCGATCGATCGATCGA'), # 1 base insertion (A at end)
            ('ins3', 'ATCGATCGATCGATCGAAA'), # 3 base insertion (AAA at end)
            ('complex', 'ATCGATCCCCGATCG'),  # Complex: 3 SNPs in middle
        ],
        'locus2': [
            ('ref', 'GGGGGGGGGGGGGGGG'),  # Reference
            ('identical', 'GGGGGGGGGGGGGGGG'),  # Identical to ref
            ('snp1', 'GGGGGGGGGTGGGGGG'),  # 1 SNP in middle (G→T)
            ('snp2', 'AGGGGGGGGGGGGGGG'),  # 1 SNP at start (G→A)
            ('del2', 'GGGGGGGGGGGGGG'),    # 2 base deletion at end
            ('del5', 'GGGGGGGGGGG'),       # 5 base deletion at end
            ('ins2', 'GGGGGGGGGGGGGGGGCC'), # 2 base insertion (CC at end)
            ('mixed', 'GGGGTGGGGGGGGG'),   # 1 SNP + 2 deletions
        ],
        'locus3': [
            ('ref', 'CCCCCCCCCCCCCCCC'),  # Reference
            ('snp3', 'CCCCCTTCCCCCCCCC'),  # 2 SNPs in middle (CC→TT)
            ('identical', 'CCCCCCCCCCCCCCCC'),  # Identical to ref
            ('del1', 'CCCCCCCCCCCCCCC'),   # 1 base deletion at end
            ('ins4', 'CCCCCCCCCCCCCCCCTTTT'), # 4 base insertion (TTTT at end)
            ('complex1', 'CCCAAACCCCCCCCC'), # 3 SNPs (CCC→AAA) + 1 deletion
            ('complex2', 'CCCCCCCCCCCCCCCCAAAA'), # 4 base insertion with different bases
            ('large_del', 'CCCCCCCC'),      # 8 base deletion (half the sequence)
        ]
    }
    
    # Create schema directory
    os.makedirs('schema_crc32', exist_ok=True)
    
    # Create FASTA files and collect CRC32 hashes
    hashes = {}
    
    for locus, sequences in test_cases.items():
        fasta_path = f'schema_crc32/{locus}.fasta'
        hashes[locus] = {}
        
        with open(fasta_path, 'w') as f:
            for name, seq in sequences:
                hash_val = crc32_hash(seq)
                hashes[locus][name] = hash_val
                f.write(f'>{hash_val}\n{seq}\n')
        
        print(f'Created {fasta_path} with {len(sequences)} alleles')
    
    return hashes, test_cases

def create_test_profiles(hashes):
    """Create allelic profiles matrix with specific test scenarios."""
    
    # Define test samples with specific allele combinations
    samples = [
        ('Sample_Ref', 'ref', 'ref', 'ref'),  # All reference alleles
        ('Sample_Identical', 'ref', 'identical', 'identical'),  # Should be distance 0 from Ref
        ('Sample_SNPs_Only', 'snp1', 'snp1', 'snp3'),  # Only SNPs (2+1+2=5 SNPs total)
        ('Sample_Dels_Only', 'del1', 'del2', 'del1'),  # Only deletions (1+2+1=4 del bases)
        ('Sample_Ins_Only', 'ins1', 'ins2', 'ins4'),  # Only insertions (1+2+4=7 ins bases)
        ('Sample_Mixed1', 'snp2', 'del5', 'complex1'),  # Mixed: 1 SNP + 5 dels + (3 SNPs + 1 del)
        ('Sample_Mixed2', 'del3', 'snp2', 'ins4'),  # Mixed: 3 dels + 1 SNP + 4 ins
        ('Sample_Complex', 'complex', 'mixed', 'complex2'),  # Complex patterns
        ('Sample_Large_Dels', 'del3', 'del5', 'large_del'),  # Large deletions (3+5+8=16)
        ('Sample_Large_Ins', 'ins3', 'ins2', 'complex2'),  # Large insertions (3+2+4=9)
    ]
    
    # Create profiles file
    with open('profiles/test_profiles_crc32.tsv', 'w') as f:
        f.write('sample\tlocus1\tlocus2\tlocus3\n')
        for sample_data in samples:
            sample_name = sample_data[0]
            locus1_allele = hashes['locus1'][sample_data[1]]
            locus2_allele = hashes['locus2'][sample_data[2]]
            locus3_allele = hashes['locus3'][sample_data[3]]
            f.write(f'{sample_name}\t{locus1_allele}\t{locus2_allele}\t{locus3_allele}\n')
    
    print(f'\nCreated profiles/test_profiles_crc32.tsv with {len(samples)} samples')
    
    return samples

def generate_expected_results(test_cases, samples):
    """Generate documentation of expected results."""
    
    doc = """# Expected Results for CRC32 Validation Test

## Test Sequences

### Locus 1:
"""
    for name, seq in test_cases['locus1']:
        doc += f"- {name}: `{seq}`\n"
    
    doc += "\n### Locus 2:\n"
    for name, seq in test_cases['locus2']:
        doc += f"- {name}: `{seq}`\n"
    
    doc += "\n### Locus 3:\n"
    for name, seq in test_cases['locus3']:
        doc += f"- {name}: `{seq}`\n"
    
    doc += """

## Expected Key Distances

### Sample_Ref vs Sample_Identical:
- **Expected**: Hamming=0, SNPs=0, SNPs+InDel_events=0, SNPs+InDel_bases=0
- **Reason**: All alleles are identical

### Sample_Ref vs Sample_SNPs_Only:
- Locus1: ref vs snp1 = 2 SNPs (CG→GG at positions 14-15)
- Locus2: ref vs snp1 = 1 SNP (G→T at position 9)
- Locus3: ref vs snp3 = 2 SNPs (CC→TT at positions 6-7)
- **Expected**: Hamming=3, SNPs=5, SNPs+InDel_events=5, SNPs+InDel_bases=5

### Sample_Ref vs Sample_Dels_Only:
- Locus1: ref vs del1 = 0 SNPs, 1 deletion event, 1 deletion base
- Locus2: ref vs del2 = 0 SNPs, 1 deletion event, 2 deletion bases
- Locus3: ref vs del1 = 0 SNPs, 1 deletion event, 1 deletion base
- **Expected**: Hamming=3, SNPs=3 (Hamming fallback), SNPs+InDel_events=6, SNPs+InDel_bases=7

### Sample_Ref vs Sample_Ins_Only:
- Locus1: ref vs ins1 = 0 SNPs, 1 insertion event, 1 insertion base
- Locus2: ref vs ins2 = 0 SNPs, 1 insertion event, 2 insertion bases
- Locus3: ref vs ins4 = 0 SNPs, 1 insertion event, 4 insertion bases
- **Expected**: Hamming=3, SNPs=3 (Hamming fallback), SNPs+InDel_events=6, SNPs+InDel_bases=10

### Sample_Ref vs Sample_Mixed1:
- Locus1: ref vs snp2 = 1 SNP (A→G at position 0)
- Locus2: ref vs del5 = 0 SNPs, 1 deletion event, 5 deletion bases
- Locus3: ref vs complex1 = 3 SNPs + 1 deletion event + 1 deletion base
- **Expected**: Hamming=3, SNPs=5, SNPs+InDel_events=7, SNPs+InDel_bases=11

### Sample_SNPs_Only vs Sample_Dels_Only:
- All loci differ, each with different mutation types
- **Expected**: Hamming=3, complex alignment results

## Notes:
- Hamming fallback: When SNPs=0 but InDels exist, SNPs mode returns 1 per locus
- The exact alignment results depend on Parasail's global alignment algorithm
- Gap penalties affect how InDels are counted as events vs bases
"""
    
    with open('EXPECTED_RESULTS_CRC32.md', 'w') as f:
        f.write(doc)
    
    print('\nCreated EXPECTED_RESULTS_CRC32.md')

def main():
    """Main function to create the validation test."""
    print("Creating cgDist validation test with CRC32 hashes...")
    
    # Create test cases
    hashes, test_cases = create_test_cases()
    
    # Create profiles
    samples = create_test_profiles(hashes)
    
    # Generate expected results documentation
    generate_expected_results(test_cases, samples)
    
    print("\n✅ Validation test created successfully!")
    print("\nTo run the test:")
    print("  ../target/release/cgdist --schema schema_crc32 --profiles profiles/test_profiles_crc32.tsv --output results/test.tsv --mode snps --hasher-type crc32")

if __name__ == '__main__':
    main()