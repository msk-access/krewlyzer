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

# For real FSC runs, you would need to provide a small .bed.gz and bin file for testing.
# Example (pseudo):
# def test_fsc_dry(tmp_path):
#     bedgz = "tests/data/test.bed.gz"
#     out = tmp_path / "out"
#     out.mkdir()
#     result = runner.invoke(app, ["fsc", str(bedgz), "-o", str(out)])
#     assert result.exit_code == 0
#     # Check output files exist, etc.
