import pytest
from pathlib import Path
import pysam
import pandas as pd
from krewlyzer.cli import app
from typer.testing import CliRunner

def create_mock_bam_uxm(path):
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 5000, 'SN': 'chr1'}]}
    
    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        # Read 1: Methylated (M) - XM: ZZZ (all methylated)
        a = pysam.AlignedSegment()
        a.query_name = "read_M"
        a.query_sequence = "CGATA"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 1000
        a.mapping_quality = 60
        a.cigar = ((0, 5),) 
        a.set_tag("XM", "ZZZ")  # 3 CpGs methylated
        outf.write(a)
        
        # Read 2: Unmethylated (U) - XM: zzz
        b = pysam.AlignedSegment()
        b.query_name = "read_U"
        b.query_sequence = "CGATA"
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 1000
        b.mapping_quality = 60
        b.cigar = ((0, 5),)
        b.set_tag("XM", "zzz")  # 3 CpGs unmethylated
        outf.write(b)

    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))

def test_uxm_integration(tmp_path):
    # Single BAM file (not directory)
    bam_file = tmp_path / "sample.bam"
    create_mock_bam_uxm(bam_file)
    
    mark_file = tmp_path / "markers.bed"
    with open(mark_file, "w") as f:
        f.write("chr1\t1000\t1005\n")
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    # New CLI: uxm <bam> -o <output_dir> -s <sample_name> -m <markers>
    result = runner.invoke(app, [
        "uxm", str(bam_file),
        "-o", str(output_dir),
        "-s", "sample",
        "-m", str(mark_file)
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    # Check output
    output_file = output_dir / "sample.UXM.tsv"
    assert output_file.exists()
