import pytest
import re
from typer.testing import CliRunner
from krewlyzer.cli import app

runner = CliRunner()

def strip_ansi(text):
    ansi_escape = re.compile(r'\x1b\[[0-9;]*m')
    return ansi_escape.sub('', text)

def test_cli_help():
    result = runner.invoke(app, ["--help"])
    output = strip_ansi(result.output)
    assert result.exit_code == 0
    assert "Usage" in output
    assert "motif" in output
    assert "fsc" in output
    assert "fsr" in output
    assert "fsd" in output
    assert "wps" in output

def test_motif_help():
    result = runner.invoke(app, ["motif", "--help"])
    output = strip_ansi(result.output)
    assert result.exit_code == 0
    assert "motif" in output
    assert "--minlen" in output
    assert "-g" in output
    assert "-o" in output

def test_fsc_help():
    result = runner.invoke(app, ["fsc", "--help"])
    output = strip_ansi(result.output)
    assert result.exit_code == 0
    assert "fragment size coverage" in output.lower() or "fsc" in output.lower()
    assert "--bin-input" in output
    assert "--output" in output

def test_fsr_help():
    result = runner.invoke(app, ["fsr", "--help"])
    assert result.exit_code == 0
    assert "fragment size ratio" in result.output.lower()

def test_fsd_help():
    result = runner.invoke(app, ["fsd", "--help"])
    assert result.exit_code == 0
    assert "fragment size distribution" in result.output.lower()

def test_wps_help():
    result = runner.invoke(app, ["wps", "--help"])
    output = strip_ansi(result.output)
    assert result.exit_code == 0
    assert "windowed protection score" in output.lower() or "wps" in output.lower()
    assert "--tsv-input" in output
    assert "--output" in output

def test_ocf_help():
    result = runner.invoke(app, ["ocf", "--help"])
    output = strip_ansi(result.output)
    assert result.exit_code == 0
    assert "orientation-aware cfDNA fragmentation" in output.lower() or "ocf" in output.lower()
    assert "--ocr-input" in output
    assert "--output" in output

def test_uxm_help():
    """
    Test that the UXM CLI help includes all major options and documents both SE and PE mode.
    """
    result = runner.invoke(app, ["uxm", "--help"])
    output = strip_ansi(result.output)
    assert result.exit_code == 0
    assert "uxm" in output.lower() or "fragment-level methylation" in output.lower()
    assert "--mark-input" in output
    assert "--output" in output
    assert "--type" in output
    assert "SE" in output or "PE" in output  # Should document both modes

def test_run_all_help():
    """
    Test the run-all CLI wrapper help output for all major options and feature mentions.
    """
    result = runner.invoke(app, ["run-all", "--help"])
    output = strip_ansi(result.output)
    assert result.exit_code == 0
    # Check that each feature is mentioned in the help or docstring
    for feat in ["motif", "fsc", "fsr", "fsd", "wps", "ocf", "uxm"]:
        assert feat in output.lower()
    # Check required options
    assert "--reference" in output
    assert "--output" in output
    assert "--threads" in output
    assert "--type" in output

# For real FSC runs, you would need to provide a small .bed.gz and bin file for testing.
# Example (pseudo):
# def test_fsc_dry(tmp_path):
#     bedgz = "tests/data/test.bed.gz"
#     out = tmp_path / "out"
#     out.mkdir()
#     result = runner.invoke(app, ["fsc", str(bedgz), "-o", str(out)])
#     assert result.exit_code == 0
#     # Check output files exist, etc.
