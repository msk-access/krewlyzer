import pytest
from pathlib import Path
import pysam
from krewlyzer.ocf import ocf

def create_mock_bedgz(path):
    # Create a small BED file
    data = []
    for i in range(100):
        data.append(f"chr1\t{1000+i*10}\t{1000+i*10+160}\t0.5")
    
    with open(path.with_suffix(""), "w") as f:
        f.write("\n".join(data))
    
    pysam.tabix_compress(str(path.with_suffix("")), str(path), force=True)
    pysam.tabix_index(str(path), preset="bed", force=True)

def test_ocf_integration(tmp_path):
    # Setup - create single BED.gz file
    bed_file = tmp_path / "test_sample.bed.gz"
    create_mock_bedgz(bed_file)
    
    # OCR file format: chrom, start, end, label
    ocr_file = tmp_path / "ocr.bed"
    with open(ocr_file, "w") as f:
        f.write("chr1\t1000\t2000\ttissue1\n")
    
    output_dir = tmp_path / "output"
    
    # Run CLI
    from typer.testing import CliRunner
    from krewlyzer.cli import app
    runner = CliRunner()
    
    # New CLI: ocf <input.bed.gz> -o <output_dir> --sample-name <name> -r <ocr>
    result = runner.invoke(app, [
        "ocf", str(bed_file), 
        "-o", str(output_dir),
        "-s", "test_sample",
        "-r", str(ocr_file)
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    # Verify output files exist
    ocf_file = output_dir / "test_sample.OCF.csv"
    sync_file = output_dir / "test_sample.OCF.sync.tsv"
    # At least one should exist
    assert ocf_file.exists() or sync_file.exists(), "No OCF output files generated"
