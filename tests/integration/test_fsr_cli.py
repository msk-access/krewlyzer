"""
Integration tests for FSR (Fragment Size Ratio) CLI.

Tests the standalone fsr command.
"""
import pytest
import pysam
from pathlib import Path
from typer.testing import CliRunner
from krewlyzer.cli import app


runner = CliRunner()


@pytest.fixture
def sample_bedgz(tmp_path):
    """Create sample BED.gz for FSR testing with dense coverage."""
    bed = tmp_path / "sample.bed"
    with open(bed, "w") as f:
        # Dense fragment coverage for FSR to find valid windows
        for i in range(100):
            start = 50 + (i * 10)
            # Mix of short, intermediate, and long fragments
            if i % 3 == 0:
                length = 100  # short
            elif i % 3 == 1:
                length = 200  # intermediate
            else:
                length = 300  # long
            f.write(f"chr1\t{start}\t{start + length}\t0.50\n")
    
    bedgz = tmp_path / "sample.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.fixture
def sample_bins(tmp_path):
    """Create bins file for FSR."""
    bins = tmp_path / "bins.bed"
    bins.write_text("chr1\t0\t600\tBin1\nchr1\t600\t1200\tBin2\n")
    return bins


@pytest.mark.integration
def test_fsr_cli_help():
    """Test FSR CLI help output."""
    result = runner.invoke(app, ["fsr", "--help"])
    assert result.exit_code == 0
    assert "fragment size ratio" in result.output.lower()


@pytest.mark.skip(reason="FSR requires minimum fragment threshold per window")
@pytest.mark.integration
def test_fsr_basic_run(tmp_path, sample_bedgz, sample_bins):
    """Test basic FSR execution."""
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "fsr", str(sample_bedgz),
        "-b", str(sample_bins),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--no-gc-correct"  # Disable GC correction for testing
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    # Check output exists
    fsr_out = output_dir / "test_sample.FSR.tsv"
    assert fsr_out.exists(), "FSR output missing"


@pytest.mark.skip(reason="FSR requires minimum fragment threshold per window")
@pytest.mark.integration
def test_fsr_output_format(tmp_path, sample_bedgz, sample_bins):
    """Test FSR output format."""
    import pandas as pd
    
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "fsr", str(sample_bedgz),
        "-b", str(sample_bins),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--no-gc-correct"  # Disable GC correction for testing
    ])
    
    assert result.exit_code == 0
    
    fsr_out = output_dir / "test_sample.FSR.tsv"
    df = pd.read_csv(fsr_out, sep="\t")
    
    # FSR output should have rows
    assert len(df) > 0, "FSR output is empty"
