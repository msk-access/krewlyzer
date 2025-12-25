"""
Integration tests for WPS (Windowed Protection Score) CLI.

Tests the standalone wps command.
"""
import pytest
import pysam
from pathlib import Path
from typer.testing import CliRunner
from krewlyzer.cli import app


runner = CliRunner()


@pytest.fixture
def sample_bedgz(tmp_path):
    """Create sample BED.gz for WPS testing with fragments spanning gene."""
    bed = tmp_path / "sample.bed"
    with open(bed, "w") as f:
        # Many fragments spanning TSS (at 1000) and TES (at 2000)
        # WPS needs fragments that span windows for protection scores
        for i in range(50):
            start = 800 + (i * 30)  # Dense coverage from 800 to 2300
            length = 160  # Long fragments for WPS
            f.write(f"chr1\t{start}\t{start + length}\t0.50\n")
    
    bedgz = tmp_path / "sample.bed.gz"
    pysam.tabix_compress(str(bed), str(bedgz), force=True)
    pysam.tabix_index(str(bedgz), preset="bed", force=True)
    return bedgz


@pytest.fixture
def sample_transcripts(tmp_path):
    """Create transcripts TSV for WPS."""
    tsv = tmp_path / "transcripts.tsv"
    tsv.write_text("Gene1\tchr1\t1000\t2000\t+\n")
    return tsv


@pytest.mark.integration
def test_wps_cli_help():
    """Test WPS CLI help output."""
    result = runner.invoke(app, ["wps", "--help"])
    assert result.exit_code == 0
    assert "windowed protection score" in result.output.lower() or "wps" in result.output.lower()


@pytest.mark.integration
def test_wps_basic_run(tmp_path, sample_bedgz, sample_transcripts):
    """Test basic WPS execution."""
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "wps", str(sample_bedgz),
        "--tsv-input", str(sample_transcripts),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    # Check output exists (may be .gz)
    wps_out = output_dir / "test_sample.WPS.tsv.gz"
    wps_out_tsv = output_dir / "test_sample.WPS.tsv"
    assert wps_out.exists() or wps_out_tsv.exists(), "WPS output missing"


@pytest.mark.integration
def test_wps_output_columns(tmp_path, sample_bedgz, sample_transcripts):
    """Test WPS output has expected columns."""
    import pandas as pd
    
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "wps", str(sample_bedgz),
        "--tsv-input", str(sample_transcripts),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0
    
    # Try both formats
    wps_out = output_dir / "test_sample.WPS.tsv.gz"
    if wps_out.exists():
        df = pd.read_csv(wps_out, sep="\t")
    else:
        wps_out = output_dir / "test_sample.WPS.tsv"
        df = pd.read_csv(wps_out, sep="\t")
    
    # WPS output should have gene_id column
    assert "gene_id" in df.columns
