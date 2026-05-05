# Expected Results for CRC32 Validation Test

## Test Sequences

### Locus 1:
- ref: `ATCGATCGATCGATCG`
- snp1: `ATCGATCGATCGATGG`
- snp2: `GTCGATCGATCGATCG`
- del1: `ATCGATCGATCGTCG`
- del3: `ATCGATCGTCG`
- ins1: `ATCGATCGATCGATCGA`
- ins3: `ATCGATCGATCGATCGAAA`
- complex: `ATCGATCCCCGATCG`

### Locus 2:
- ref: `GGGGGGGGGGGGGGGG`
- identical: `GGGGGGGGGGGGGGGG`
- snp1: `GGGGGGGGGTGGGGGG`
- snp2: `AGGGGGGGGGGGGGGG`
- del2: `GGGGGGGGGGGGGG`
- del5: `GGGGGGGGGGG`
- ins2: `GGGGGGGGGGGGGGGGCC`
- mixed: `GGGGTGGGGGGGGG`

### Locus 3:
- ref: `CCCCCCCCCCCCCCCC`
- snp3: `CCCCCTTCCCCCCCCC`
- identical: `CCCCCCCCCCCCCCCC`
- del1: `CCCCCCCCCCCCCCC`
- ins4: `CCCCCCCCCCCCCCCCTTTT`
- complex1: `CCCAAACCCCCCCCC`
- complex2: `CCCCCCCCCCCCCCCCAAAA`
- large_del: `CCCCCCCC`


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
