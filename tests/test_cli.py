import pytest
from typer.testing import CliRunner
from krewlyzer.cli import app

runner = CliRunner()

def test_cli_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Usage" in result.output
    assert "motif" in result.output
    assert "fsc" in result.output
    assert "fsr" in result.output
    assert "fsd" in result.output
    assert "wps" in result.output

def test_motif_help():
    result = runner.invoke(app, ["motif", "--help"])
    assert result.exit_code == 0
    assert "motif" in result.output
    assert "--minlen" in result.output
    assert "-g" in result.output
    assert "-o" in result.output

def test_fsc_help():
    result = runner.invoke(app, ["fsc", "--help"])
    assert result.exit_code == 0
    assert "fragment size coverage" in result.output.lower() or "fsc" in result.output.lower()
    assert "--bin-input" in result.output
    assert "--output" in result.output

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
    assert result.exit_code == 0
    assert "windowed protection score" in result.output.lower() or "wps" in result.output.lower()
    assert "--tsv-input" in result.output
    assert "--output" in result.output

# For real FSC runs, you would need to provide a small .bed.gz and bin file for testing.
# Example (pseudo):
# def test_fsc_dry(tmp_path):
#     bedgz = "tests/data/test.bed.gz"
#     out = tmp_path / "out"
#     out.mkdir()
#     result = runner.invoke(app, ["fsc", str(bedgz), "-o", str(out)])
#     assert result.exit_code == 0
#     # Check output files exist, etc.
