"""
Integration tests for FSC (Fragment Size Coverage) CLI.

Tests the standalone fsc command.
"""
import pytest
import pysam
import gzip
from pathlib import Path
from typer.testing import CliRunner
from krewlyzer.cli import app


runner = CliRunner()


@pytest.fixture
def sample_bedgz(tmp_path):
    """Create sample BED.gz for FSC testing."""
    bed = tmp_path / "sample.bed"
    with open(bed, "w") as f:
        # Fragments in different size categories
        f.write("chr1\t100\t180\t0.45\n")   # 80bp - ultra-short
        f.write("chr1\t200\t340\t0.50\n")   # 140bp - short
        f.write("chr1\t400\t600\t0.48\n")   # 200bp - intermediate
        f.write("chr1\t700\t1000\t0.52\n")  # 300bp - long
    
    bedgz = tmp_path / "sample.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.fixture
def sample_bins(tmp_path):
    """Create bins file for FSC."""
    bins = tmp_path / "bins.bed"
    bins.write_text("chr1\t0\t500\tBin1\nchr1\t500\t1500\tBin2\n")
    return bins


@pytest.mark.integration
def test_fsc_cli_help():
    """Test FSC CLI help output."""
    result = runner.invoke(app, ["fsc", "--help"], color=False)
    assert result.exit_code == 0
    assert "--bin-input" in result.output
    assert "--output" in result.output


@pytest.mark.integration
def test_fsc_basic_run(tmp_path, sample_bedgz, sample_bins):
    """Test basic FSC execution."""
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "fsc", str(sample_bedgz),
        "-b", str(sample_bins),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    # Check output exists
    fsc_out = output_dir / "test_sample.FSC.tsv"
    assert fsc_out.exists(), "FSC output missing"


@pytest.mark.integration
def test_fsc_output_format(tmp_path, sample_bedgz, sample_bins):
    """Test FSC output format has expected columns."""
    import pandas as pd
    
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "fsc", str(sample_bedgz),
        "-b", str(sample_bins),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0
    
    fsc_out = output_dir / "test_sample.FSC.tsv"
    df = pd.read_csv(fsc_out, sep="\t")
    
    # Check expected columns (FSC z-score outputs)
    expected_cols = ["region", "short-fragment-zscore"]
    for col in expected_cols:
        assert col in df.columns, f"Missing column: {col}"
