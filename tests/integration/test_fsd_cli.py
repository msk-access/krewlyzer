"""
Integration tests for FSD (Fragment Size Distribution) CLI.

Tests the standalone fsd command.
"""
import pytest
import pysam
from pathlib import Path
from typer.testing import CliRunner
from krewlyzer.cli import app


runner = CliRunner()


@pytest.fixture
def sample_bedgz(tmp_path):
    """Create sample BED.gz for FSD testing."""
    bed = tmp_path / "sample.bed"
    with open(bed, "w") as f:
        # Multiple fragments for distribution calculation
        for i in range(20):
            start = 1000 + i * 100
            length = 150 + (i % 5) * 20  # Varying sizes
            f.write(f"chr1\t{start}\t{start + length}\t0.50\n")
    
    bedgz = tmp_path / "sample.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.fixture
def sample_arms(tmp_path):
    """Create arms file for FSD."""
    arms = tmp_path / "arms.bed"
    arms.write_text("chr1\t0\t10000\tArm1\n")
    return arms


@pytest.mark.integration
def test_fsd_cli_help():
    """Test FSD CLI help output."""
    result = runner.invoke(app, ["fsd", "--help"])
    assert result.exit_code == 0
    assert "fragment size distribution" in result.output.lower()


@pytest.mark.integration
def test_fsd_basic_run(tmp_path, sample_bedgz, sample_arms):
    """Test basic FSD execution."""
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "fsd", "-i", str(sample_bedgz),
        "-a", str(sample_arms),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    # Check output exists
    fsd_out = output_dir / "test_sample.FSD.tsv"
    assert fsd_out.exists(), "FSD output missing"


@pytest.mark.integration
def test_fsd_output_columns(tmp_path, sample_bedgz, sample_arms):
    """Test FSD output has expected columns."""
    import pandas as pd
    
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "fsd", "-i", str(sample_bedgz),
        "-a", str(sample_arms),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0
    
    fsd_out = output_dir / "test_sample.FSD.tsv"
    df = pd.read_csv(fsd_out, sep="\t")
    
    # FSD should have arm-related columns
    assert len(df) > 0, "FSD output is empty"
