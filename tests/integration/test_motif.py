"""
Integration tests for Motif extraction.

Tests end-motif and breakpoint-motif extraction via CLI.
"""

import pytest
import pysam
from typer.testing import CliRunner
from krewlyzer.cli import app

runner = CliRunner()


@pytest.fixture
def mock_bam_for_motif(tmp_path):
    """Create BAM with proper pairs for motif extraction."""
    bam = tmp_path / "test.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 2000, "SN": "chr1"}]}

    with pysam.AlignmentFile(str(bam), "wb", header=header) as outf:
        # Proper pair
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "ACGT" * 25  # 100bp
        a.flag = 99
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigar = ((0, 100),)
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        outf.write(a)

        b = pysam.AlignedSegment()
        b.query_name = "read1"
        b.query_sequence = "TGCA" * 25
        b.flag = 147
        b.reference_id = 0
        b.reference_start = 200
        b.mapping_quality = 60
        b.cigar = ((0, 100),)
        b.next_reference_id = 0
        b.next_reference_start = 100
        b.template_length = -167
        outf.write(b)

    pysam.index(str(bam))
    return bam


@pytest.fixture
def mock_reference(tmp_path):
    """Create reference FASTA."""
    ref = tmp_path / "genome.fa"
    seq = "ACGTACGT" * 250  # 2kb
    ref.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(ref))
    return ref


import re


def strip_ansi(text: str) -> str:
    """Remove ANSI escape codes from text."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


@pytest.mark.integration
def test_motif_cli_help():
    """Test motif CLI help output."""
    result = runner.invoke(app, ["motif", "--help"])
    assert result.exit_code == 0
    output = strip_ansi(result.output)
    assert "--kmer" in output
    assert "-r" in output


@pytest.mark.integration
def test_motif_extraction(tmp_path, mock_bam_for_motif, mock_reference):
    """Test basic motif extraction."""
    output_dir = tmp_path / "output"

    result = runner.invoke(
        app,
        [
            "motif",
            "-i",
            str(mock_bam_for_motif),
            "-r",
            str(mock_reference),
            "-o",
            str(output_dir),
            "-s",
            "test",
            "--kmer",
            "4",
        ],
    )

    assert result.exit_code == 0, f"CLI failed: {result.output}"

    # Check outputs exist
    assert (output_dir / "test.EndMotif.tsv").exists()
    assert (output_dir / "test.BreakPointMotif.tsv").exists()
    assert (output_dir / "test.MDS.tsv").exists()


@pytest.mark.integration
def test_motif_mds_score(tmp_path, mock_bam_for_motif, mock_reference):
    """Test MDS score calculation."""
    output_dir = tmp_path / "output"

    result = runner.invoke(
        app,
        [
            "motif",
            "-i",
            str(mock_bam_for_motif),
            "-r",
            str(mock_reference),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    assert result.exit_code == 0

    mds_file = output_dir / "test.MDS.tsv"
    assert mds_file.exists()

    # MDS should be between 0 and 1
    import pandas as pd

    df = pd.read_csv(mds_file, sep="\t")
    if "mds" in df.columns:
        mds = df["mds"].iloc[0]
        assert 0 <= mds <= 1
