import pytest
from typer.testing import CliRunner
from krewlyzer.cli import app

runner = CliRunner()

def test_cli_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Usage" in result.output
    assert "motif" in result.output
    assert "extract-features" in result.output or "extract_features" in result.output
    assert "quality-control" in result.output or "quality_control" in result.output

def test_motif_help():
    result = runner.invoke(app, ["motif", "--help"])
    assert result.exit_code == 0
    assert "motif" in result.output
    assert "--minlen" in result.output
    assert "-g" in result.output
    assert "-o" in result.output

# For real motif runs, you would need to provide a small BAM and genome file for testing.
# Example (pseudo):
# def test_motif_dry(tmp_path):
#     bam = "tests/data/test.bam"
#     genome = "tests/data/test.fa"
#     out = tmp_path / "out"
#     result = runner.invoke(app, ["motif", str(bam), "-g", str(genome), "-o", str(out)])
#     assert result.exit_code == 0
#     assert (out / "EDM").exists()
