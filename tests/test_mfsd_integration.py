import pytest
from pathlib import Path
import pysam
import pandas as pd
from krewlyzer.cli import app
from typer.testing import CliRunner

def create_mock_bam(path):
    # Create BAM with header
    header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 5000, 'SN': 'chr1'}] }
    
    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        # Create mutant read at pos 1000 (0-based)
        # Variant: chr1:1001 A->T (1-based) => 0-based 1000.
        # Ref A, Alt T.
        
        # Read 1: Mutant (T)
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "TTTTT" # T at pos 0
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 1000
        a.mapping_quality = 60
        a.cigar = ((0, 5),) # 5M
        a.template_length = 150
        outf.write(a)
        
        # Read 2: Wildtype (A)
        b = pysam.AlignedSegment()
        b.query_name = "read2"
        b.query_sequence = "AAAAA" # A at pos 0
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 1000
        b.mapping_quality = 60
        b.cigar = ((0, 5),) 
        b.template_length = 200
        outf.write(b)

    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))

def create_mock_vcf(path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        # chr1, 1001, ., A, T
        f.write("chr1\t1001\t.\tA\tT\t.\t.\t.\n")

def test_mfsd_integration(tmp_path):
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    # New CLI: mfsd <bam> -i <variants> -o <output_dir> -s <sample_name>
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test"
    ])
    
    # Debug info if fails
    if result.exit_code != 0:
        print(result.stdout)
        print(result.exception)
        
    assert result.exit_code == 0
    
    output_file = output_dir / "test.mFSD.tsv"
    assert output_file.exists()
