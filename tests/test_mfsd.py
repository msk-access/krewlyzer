import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import pysam
from krewlyzer.mfsd import calc_mfsd, classify_read, parse_input_file

# Mock pysam.AlignedSegment
class MockRead:
    def __init__(self, query_sequence, reference_start, reference_end, template_length, mapping_quality=60, is_duplicate=False, is_unmapped=False, is_secondary=False):
        self.query_sequence = query_sequence
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.template_length = template_length
        self.mapping_quality = mapping_quality
        self.is_duplicate = is_duplicate
        self.is_unmapped = is_unmapped
        self.is_secondary = is_secondary
        self.query_length = len(query_sequence)
        
    def get_aligned_pairs(self, matches_only=True):
        # Simple 1-to-1 mapping for mock
        pairs = []
        for i in range(len(self.query_sequence)):
            pairs.append((i, self.reference_start + i))
        return pairs

def test_classify_read():
    # Ref: A, Alt: T at pos 10
    # Read covers 0-20.
    # Base at 10 is T -> Mutant
    read_mut = MockRead("A" * 10 + "T" + "A" * 9, 0, 20, 150)
    assert classify_read(read_mut, 10, "A", "T") == "Mutant"
    
    # Base at 10 is A -> WildType
    read_wt = MockRead("A" * 20, 0, 20, 160)
    assert classify_read(read_wt, 10, "A", "T") == "WildType"
    
    # Base at 10 is C -> Other
    read_other = MockRead("A" * 10 + "C" + "A" * 9, 0, 20, 150)
    assert classify_read(read_other, 10, "A", "T") == "Other"
    
    # Read does not cover pos
    read_out = MockRead("A" * 10, 20, 30, 150)
    assert classify_read(read_out, 10, "A", "T") == "Unknown"

def test_parse_input_file(tmp_path):
    # Test VCF
    vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\t.\t.
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    
    df_vcf = parse_input_file(vcf_file, "vcf")
    assert len(df_vcf) == 1
    assert df_vcf.iloc[0]['chrom'] == "chr1"
    assert df_vcf.iloc[0]['pos'] == 99 # 0-based
    assert df_vcf.iloc[0]['ref'] == "A"
    assert df_vcf.iloc[0]['alt'] == "T"
    
    # Test MAF
    maf_content = """#version 2.4
Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele2
GENE\tchr2\t200\tC\tG
"""
    maf_file = tmp_path / "test.maf"
    maf_file.write_text(maf_content)
    
    df_maf = parse_input_file(maf_file, "maf")
    assert len(df_maf) == 1
    assert df_maf.iloc[0]['chrom'] == "chr2"
    assert df_maf.iloc[0]['pos'] == 199 # 0-based
    assert df_maf.iloc[0]['ref'] == "C"
    assert df_maf.iloc[0]['alt'] == "G"

def test_calc_mfsd_integration(tmp_path, mocker):
    # Mock pysam.AlignmentFile and fetch
    mock_bam = mocker.Mock()
    mocker.patch("pysam.AlignmentFile", return_value=mock_bam)
    
    # Mock reads at chr1:99 (pos 100 in VCF)
    # Mutant reads: short (100bp)
    # WT reads: long (160bp)
    
    reads = []
    # 10 Mutant reads
    for _ in range(10):
        reads.append(MockRead("A" * 99 + "T" + "A" * 50, 0, 150, 100)) # covers 99 with T
    # 10 WT reads
    for _ in range(10):
        reads.append(MockRead("A" * 150, 0, 150, 160)) # covers 99 with A (ref)
        
    mock_bam.fetch.return_value = reads
    
    # Input VCF
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\t.\t.
""")
    
    output_file = tmp_path / "output.tsv"
    
    calc_mfsd(Path("dummy.bam"), vcf_file, output_file, input_format="vcf")
    
    # Verify output
    df = pd.read_csv(output_file, sep='\t')
    assert len(df) == 1
    row = df.iloc[0]
    assert row['Chrom'] == "chr1"
    assert row['Pos'] == 100
    assert row['Mut_Count'] == 10
    assert row['WT_Count'] == 10
    assert row['Mut_MeanSize'] == 100.0
    assert row['WT_MeanSize'] == 160.0
    assert row['Delta_Size'] == 60.0
    # KS test should be significant (p < 0.05) as distributions are disjoint
    assert row['KS_Pval'] < 0.001
