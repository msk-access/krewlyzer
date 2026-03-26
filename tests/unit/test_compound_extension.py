"""
Tests for compound-extension path safety.

Verifies that write_table/read_table correctly handle compound base names
like 'sample.MDS.exon', 'sample.EndMotif', 'sample.MDS.ontarget' —
which are vulnerable to Python's Path.with_suffix() bug.

Bug reference: with_suffix() replaces the LAST dot-segment, so:
    Path("P-XXX.MDS.exon").with_suffix(".tsv") → "P-XXX.MDS.tsv"  (WRONG)
    Path("P-XXX.EndMotif").with_suffix(".tsv") → "P-XXX.tsv"      (WRONG)

The safe pattern is: base.parent / (base.name + ".tsv")
"""

import pandas as pd
import pytest
from pathlib import Path

from krewlyzer.core.output_utils import write_table, read_table

# ── Compound base names used throughout krewlyzer ────────────────────────────

COMPOUND_BASES = [
    "sample.MDS.exon",
    "sample.MDS.gene",
    "sample.MDS.ontarget",
    "sample.EndMotif",
    "sample.EndMotif.ontarget",
    "sample.BreakPointMotif",
    "sample.BreakPointMotif.ontarget",
    "sample.EndMotif1mer",
    "sample.FSC.gene",
    "sample.FSC.regions",
    "sample.FSC.regions.e1only",
    "sample.OCF.sync",
    "sample.TFBS.sync",
]


@pytest.fixture
def df():
    """Minimal DataFrame for testing output format conversions."""
    return pd.DataFrame({"gene": ["TP53", "BRCA1"], "mds": [0.85, 0.92]})


# ── write_table compound-extension tests ─────────────────────────────────────


class TestCompoundExtensionSafety:
    """Verify write_table produces correct filenames for compound bases."""

    @pytest.mark.unit
    @pytest.mark.parametrize("base_name", COMPOUND_BASES)
    def test_tsv_output(self, tmp_path, df, base_name):
        """write_table with TSV format produces {base}.tsv, not {truncated}.tsv."""
        base = tmp_path / base_name
        write_table(df, base, output_format="tsv")
        expected = tmp_path / (base_name + ".tsv")
        assert (
            expected.exists()
        ), f"Expected {expected.name}, got: {[f.name for f in tmp_path.iterdir()]}"

    @pytest.mark.unit
    @pytest.mark.parametrize("base_name", COMPOUND_BASES)
    def test_parquet_output(self, tmp_path, df, base_name):
        """write_table with Parquet format produces {base}.parquet."""
        base = tmp_path / base_name
        write_table(df, base, output_format="parquet")
        expected = tmp_path / (base_name + ".parquet")
        assert (
            expected.exists()
        ), f"Expected {expected.name}, got: {[f.name for f in tmp_path.iterdir()]}"

    @pytest.mark.unit
    @pytest.mark.parametrize("base_name", COMPOUND_BASES)
    def test_both_output(self, tmp_path, df, base_name):
        """write_table with 'both' format produces {base}.tsv AND {base}.parquet."""
        base = tmp_path / base_name
        write_table(df, base, output_format="both")
        tsv_expected = tmp_path / (base_name + ".tsv")
        pq_expected = tmp_path / (base_name + ".parquet")
        assert tsv_expected.exists(), f"Missing {tsv_expected.name}"
        assert pq_expected.exists(), f"Missing {pq_expected.name}"

    @pytest.mark.unit
    @pytest.mark.parametrize("base_name", COMPOUND_BASES)
    def test_roundtrip(self, tmp_path, df, base_name):
        """write_table → read_table roundtrip preserves data for compound bases."""
        base = tmp_path / base_name
        write_table(df, base, output_format="both")
        tsv_path = tmp_path / (base_name + ".tsv")
        result = read_table(tsv_path)
        assert result is not None, f"read_table returned None for {tsv_path.name}"
        assert len(result) == len(df)
        assert list(result.columns) == list(df.columns)


# ── Region MDS path distinctness tests ───────────────────────────────────────


class TestRegionMdsOutputPaths:
    """Verify region_mds exon and gene paths are distinct."""

    @pytest.mark.unit
    def test_exon_gene_paths_distinct(self, tmp_path):
        """Safe pattern produces distinct exon and gene output paths."""
        sample = "P-0000113-T02-XS1"
        exon_base = tmp_path / f"{sample}.MDS.exon"
        gene_base = tmp_path / f"{sample}.MDS.gene"

        exon_path = exon_base.parent / (exon_base.name + ".tsv")
        gene_path = gene_base.parent / (gene_base.name + ".tsv")

        assert exon_path != gene_path, "Exon and gene paths must be distinct"
        assert exon_path.name == f"{sample}.MDS.exon.tsv"
        assert gene_path.name == f"{sample}.MDS.gene.tsv"

    @pytest.mark.unit
    def test_with_suffix_would_collide(self):
        """Document the bug: with_suffix() makes exon and gene paths identical."""
        exon_base = Path("P-XXX.MDS.exon")
        gene_base = Path("P-XXX.MDS.gene")

        # with_suffix replaces last dot-segment: both → P-XXX.MDS.tsv
        assert exon_base.with_suffix(".tsv") == gene_base.with_suffix(".tsv"), (
            "This test documents the Python with_suffix() bug: "
            "both .MDS.exon and .MDS.gene collapse to .MDS.tsv"
        )


# ── Motif tracking path tests ───────────────────────────────────────────────


class TestMotifTrackingPaths:
    """Verify motif output tracking paths resolve correctly."""

    @pytest.mark.unit
    def test_motif_paths_all_distinct(self, tmp_path):
        """EndMotif, BreakPointMotif, and MDS tracking paths are distinct."""
        sample = "P-0000113-T02-XS1"
        ext = ".tsv"

        edm_base = tmp_path / f"{sample}.EndMotif"
        bpm_base = tmp_path / f"{sample}.BreakPointMotif"
        mds_base = tmp_path / f"{sample}.MDS"

        edm_out = edm_base.parent / (edm_base.name + ext)
        bpm_out = bpm_base.parent / (bpm_base.name + ext)
        mds_out = mds_base.parent / (mds_base.name + ext)

        assert edm_out.name == f"{sample}.EndMotif.tsv"
        assert bpm_out.name == f"{sample}.BreakPointMotif.tsv"
        assert mds_out.name == f"{sample}.MDS.tsv"
        assert len({edm_out, bpm_out, mds_out}) == 3, "All three paths must be distinct"

    @pytest.mark.unit
    def test_with_suffix_would_collapse_motifs(self):
        """Document the bug: with_suffix(ext) makes all motif paths identical."""
        edm = Path("P-XXX.EndMotif")
        bpm = Path("P-XXX.BreakPointMotif")
        mds = Path("P-XXX.MDS")

        # with_suffix replaces .EndMotif / .BreakPointMotif / .MDS → all become P-XXX.tsv
        assert (
            edm.with_suffix(".tsv")
            == bpm.with_suffix(".tsv")
            == mds.with_suffix(".tsv")
        )

    @pytest.mark.unit
    def test_ontarget_mds_path(self):
        """On-target MDS path must NOT collapse to genome-wide MDS path."""
        mds_on_base = Path("P-XXX.MDS.ontarget")
        mds_base = Path("P-XXX.MDS")

        # Safe pattern: distinct paths
        mds_on_out = mds_on_base.parent / (mds_on_base.name + ".tsv")
        mds_out = mds_base.parent / (mds_base.name + ".tsv")

        assert mds_on_out != mds_out
        assert mds_on_out.name == "P-XXX.MDS.ontarget.tsv"
        assert mds_out.name == "P-XXX.MDS.tsv"

        # with_suffix strips ".ontarget" → "P-XXX.MDS.tsv" (wrong: reads global MDS)
        assert mds_on_base.with_suffix(".tsv").name == "P-XXX.MDS.tsv"
        assert mds_on_base.with_suffix(".tsv").name != "P-XXX.MDS.ontarget.tsv"
