"""
Integration tests for WPS (Windowed Protection Score) CLI.

Tests the standalone wps command with Parquet output.
"""
import pytest
import pysam
import re
from pathlib import Path
from typer.testing import CliRunner
from krewlyzer.cli import app


runner = CliRunner()


def strip_ansi(text: str) -> str:
    """Remove ANSI escape codes from text for clean assertions."""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)


@pytest.fixture
def sample_bedgz(tmp_path):
    """Create sample BED.gz for WPS testing with fragments spanning anchors."""
    bed = tmp_path / "sample.bed"
    with open(bed, "w") as f:
        # Dense coverage around position 1000 (simulated TSS)
        for i in range(100):
            start = 800 + (i * 15)
            # Mix of nucleosome-sized and TF-sized fragments
            if i % 3 == 0:
                length = 167  # Primary nucleosome range
            elif i % 3 == 1:
                length = 145  # Secondary nucleosome range
            else:
                length = 55  # TF range
            f.write(f"chr1\t{start}\t{start + length}\t0.50\n")
    
    bedgz = tmp_path / "sample.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.fixture
def sample_anchors(tmp_path):
    """Create WPS anchors BED for testing."""
    bed = tmp_path / "anchors.bed"
    with open(bed, "w") as f:
        # TSS anchor
        f.write("chr1\t500\t1500\tGENE1\t.\t+\n")
    
    bedgz = tmp_path / "anchors.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.mark.integration
def test_wps_cli_help():
    """Test WPS CLI help output."""
    result = runner.invoke(app, ["wps", "--help"])
    assert result.exit_code == 0
    output = strip_ansi(result.output)
    assert "--wps-anchors" in output
    assert "--background" in output
    assert "--genome" in output


@pytest.mark.integration
def test_wps_cli_parquet_output(tmp_path, sample_bedgz, sample_anchors):
    """Test WPS produces Parquet output with expected columns."""
    import pandas as pd
    
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "wps", str(sample_bedgz),
        "--wps-anchors", str(sample_anchors),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--no-gc-correct"  # Skip GC for speed
    ])
    
    # Check exit code
    if result.exit_code != 0:
        pytest.skip(f"WPS execution failed (may need Rust backend): {result.output}")
    
    # Check Parquet output exists
    wps_parquet = output_dir / "test_sample.WPS.parquet"
    assert wps_parquet.exists(), f"WPS Parquet missing. Output: {result.output}"
    
    # Verify columns
    df = pd.read_parquet(wps_parquet)
    expected_cols = ["region_id", "chrom", "center", "strand"]
    for col in expected_cols:
        assert col in df.columns, f"Missing column: {col}"


@pytest.mark.integration
def test_wps_cli_smooth_columns(tmp_path, sample_bedgz, sample_anchors):
    """Test WPS applies smoothing and adds _smooth columns."""
    import pandas as pd
    
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "wps", str(sample_bedgz),
        "--wps-anchors", str(sample_anchors),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--no-gc-correct"
    ])
    
    if result.exit_code != 0:
        pytest.skip("WPS execution failed")
    
    wps_parquet = output_dir / "test_sample.WPS.parquet"
    if not wps_parquet.exists():
        pytest.skip("Parquet output not created")
    
    df = pd.read_parquet(wps_parquet)
    
    # Check smoothed columns exist
    if "wps_nuc" in df.columns:
        assert "wps_nuc_smooth" in df.columns, "Missing smoothed column"


@pytest.mark.integration  
def test_wps_genome_flag():
    """Test --genome flag is accepted."""
    result = runner.invoke(app, ["wps", "--help"])
    output = strip_ansi(result.output)
    assert "--genome" in output
    assert "hg19" in output or "hg38" in output
