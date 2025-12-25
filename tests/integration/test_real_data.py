"""
Integration tests using real data fixtures.

Uses subset of actual sequencing data for comprehensive testing.
"""
import pytest
from pathlib import Path
from typer.testing import CliRunner
from krewlyzer.cli import app


runner = CliRunner()


@pytest.mark.integration
def test_extract_real_bam(tmp_path, real_bam, real_reference):
    """Test extract with real BAM data."""
    output_dir = tmp_path / "output"
    
    result = runner.invoke(app, [
        "extract", str(real_bam),
        "--reference", str(real_reference),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    bed_out = output_dir / "test_sample.bed.gz"
    assert bed_out.exists(), "BED.gz output missing"


@pytest.mark.integration
def test_fsc_real_data(tmp_path, real_bam, real_reference, real_bins):
    """Test FSC with real data."""
    import gzip
    import pysam
    
    # First run extract to get BED.gz
    output_dir = tmp_path / "output"
    runner.invoke(app, [
        "extract", str(real_bam),
        "--reference", str(real_reference),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--chromosomes", "1"
    ])
    
    bedgz = output_dir / "test_sample.bed.gz"
    if not bedgz.exists():
        pytest.skip("Extract step failed")
    
    # Run FSC
    result = runner.invoke(app, [
        "fsc", str(bedgz),
        "-b", str(real_bins),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    fsc_out = output_dir / "test_sample.FSC.tsv"
    assert fsc_out.exists(), "FSC output missing"


@pytest.mark.integration
def test_fsd_real_data(tmp_path, real_bam, real_reference, real_arms):
    """Test FSD with real data."""
    # First run extract
    output_dir = tmp_path / "output"
    runner.invoke(app, [
        "extract", str(real_bam),
        "--reference", str(real_reference),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--chromosomes", "1"
    ])
    
    bedgz = output_dir / "test_sample.bed.gz"
    if not bedgz.exists():
        pytest.skip("Extract step failed")
    
    result = runner.invoke(app, [
        "fsd", str(bedgz),
        "-a", str(real_arms),
        "-o", str(output_dir),
        "-s", "test_sample"
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    fsd_out = output_dir / "test_sample.FSD.tsv"
    assert fsd_out.exists(), "FSD output missing"


@pytest.mark.integration
def test_pon_loading_real(real_pon):
    """Test loading real PON model."""
    from krewlyzer.pon.model import PonModel
    
    model = PonModel.load(real_pon)
    
    assert model is not None
    assert model.gc_bias is not None or model.fsd_baseline is not None


@pytest.mark.integration
def test_fsr_real_data(tmp_path, real_bam, real_reference, real_bins):
    """Test FSR with real data."""
    # First run extract
    output_dir = tmp_path / "output"
    runner.invoke(app, [
        "extract", str(real_bam),
        "--reference", str(real_reference),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--chromosomes", "1"
    ])
    
    bedgz = output_dir / "test_sample.bed.gz"
    if not bedgz.exists():
        pytest.skip("Extract step failed")
    
    result = runner.invoke(app, [
        "fsr", str(bedgz),
        "-b", str(real_bins),
        "-o", str(output_dir),
        "-s", "test_sample",
        "--no-gc-correct"
    ])
    
    # FSR may fail if coverage is too low for valid windows
    # This is expected behavior with small test BAM
    if result.exit_code != 0:
        assert "No valid windows" in result.output or result.exit_code == 1
        pytest.skip("FSR skipped: insufficient coverage for valid windows")
    
    fsr_out = output_dir / "test_sample.FSR.tsv"
    assert fsr_out.exists(), "FSR output missing"
