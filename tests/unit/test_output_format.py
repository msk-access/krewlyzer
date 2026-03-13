"""
Unit tests for output format and compression support across processors.

Tests that FSD, OCF, and region entropy processors correctly honour
--output-format and --compress flags by producing the expected output
files (.tsv, .tsv.gz, .parquet) when write_table() is called.

These tests use mock Rust functions to avoid needing real BAM/BED data.
"""

import pandas as pd
import pytest

from krewlyzer.core.output_utils import write_table, read_table


# ── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture
def sample_df():
    """Minimal DataFrame for testing output format conversions."""
    return pd.DataFrame({"tissue": ["A", "B", "C"], "score": [1.0, 2.0, 3.0]})


@pytest.fixture
def fsd_df():
    """Minimal FSD DataFrame (arm-level fragment distribution)."""
    return pd.DataFrame(
        {
            "arm": ["1p", "1q", "2p"],
            "short": [100, 200, 150],
            "long": [300, 400, 250],
            "ratio": [0.333, 0.500, 0.600],
        }
    )


@pytest.fixture
def ocf_df():
    """Minimal OCF DataFrame (tissue-level OCF scores)."""
    return pd.DataFrame(
        {"tissue": ["Breast", "Liver", "Lung"], "OCF": [0.12, 0.05, -0.03]}
    )


@pytest.fixture
def entropy_df():
    """Minimal region entropy DataFrame."""
    return pd.DataFrame(
        {
            "label": ["CTCF", "FOXA1", "MYC"],
            "count": [500, 300, 200],
            "mean_size": [167.5, 170.2, 165.0],
            "entropy": [3.45, 2.89, 3.12],
        }
    )


# ── write_table core tests ──────────────────────────────────────────────────


class TestWriteTableFormats:
    """Test that write_table() produces the correct files for each format."""

    @pytest.mark.unit
    def test_tsv_only(self, tmp_path, sample_df):
        """output_format='tsv' → only .tsv file, no .parquet."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="tsv", compress=False)
        assert (tmp_path / "test.tsv").exists()
        assert not (tmp_path / "test.parquet").exists()

    @pytest.mark.unit
    def test_parquet_only(self, tmp_path, sample_df):
        """output_format='parquet' → only .parquet file."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="parquet", compress=False)
        assert (tmp_path / "test.parquet").exists()

    @pytest.mark.unit
    def test_both_formats(self, tmp_path, sample_df):
        """output_format='both' → .tsv AND .parquet files."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="both", compress=False)
        assert (tmp_path / "test.tsv").exists()
        assert (tmp_path / "test.parquet").exists()

    @pytest.mark.unit
    def test_tsv_compressed(self, tmp_path, sample_df):
        """compress=True with tsv → .tsv.gz file."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="tsv", compress=True)
        assert (tmp_path / "test.tsv.gz").exists()

    @pytest.mark.unit
    def test_both_compressed(self, tmp_path, sample_df):
        """compress=True with both → .tsv.gz AND .parquet files."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="both", compress=True)
        assert (tmp_path / "test.tsv.gz").exists()
        assert (tmp_path / "test.parquet").exists()


# ── Parquet roundtrip test ───────────────────────────────────────────────────


class TestParquetRoundtrip:
    """Verify data integrity through write_table → read_table roundtrip."""

    @pytest.mark.unit
    def test_roundtrip_tsv(self, tmp_path, sample_df):
        """write_table(tsv) → read_table() preserves data."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="tsv")
        result = read_table(out)
        assert result is not None
        pd.testing.assert_frame_equal(result, sample_df)

    @pytest.mark.unit
    def test_roundtrip_parquet(self, tmp_path, sample_df):
        """write_table(parquet) → read_table() preserves data."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="parquet")
        result = read_table(tmp_path / "test.parquet")
        assert result is not None
        pd.testing.assert_frame_equal(result, sample_df)

    @pytest.mark.unit
    def test_roundtrip_gzip(self, tmp_path, sample_df):
        """write_table(tsv, compress) → read_table() preserves data."""
        out = tmp_path / "test.tsv"
        write_table(sample_df, out, output_format="tsv", compress=True)
        result = read_table(tmp_path / "test.tsv.gz")
        assert result is not None
        pd.testing.assert_frame_equal(result, sample_df)


# ── FSD processor output format tests ────────────────────────────────────────


class TestFsdOutputFormat:
    """Test that process_fsd() correctly produces requested output formats."""

    @pytest.mark.unit
    def test_fsd_both_format(self, tmp_path, fsd_df):
        """process_fsd with output_format='both' produces both .tsv and .parquet."""
        from krewlyzer.core.fsd_processor import _write_fsd_output

        tsv = tmp_path / "sample.FSD.tsv"
        fsd_df.to_csv(tsv, sep="\t", index=False)

        _write_fsd_output(tsv, output_format="both", compress=False)

        assert (tmp_path / "sample.FSD.tsv").exists()
        assert (tmp_path / "sample.FSD.parquet").exists()

    @pytest.mark.unit
    def test_fsd_compress(self, tmp_path, fsd_df):
        """process_fsd with compress=True produces .tsv.gz."""
        from krewlyzer.core.fsd_processor import _write_fsd_output

        tsv = tmp_path / "sample.FSD.tsv"
        fsd_df.to_csv(tsv, sep="\t", index=False)

        _write_fsd_output(tsv, output_format="tsv", compress=True)

        assert (tmp_path / "sample.FSD.tsv.gz").exists()

    @pytest.mark.unit
    def test_fsd_both_compress(self, tmp_path, fsd_df):
        """process_fsd with both+compress produces .tsv.gz AND .parquet."""
        from krewlyzer.core.fsd_processor import _write_fsd_output

        tsv = tmp_path / "sample.FSD.tsv"
        fsd_df.to_csv(tsv, sep="\t", index=False)

        _write_fsd_output(tsv, output_format="both", compress=True)

        assert (tmp_path / "sample.FSD.tsv.gz").exists()
        assert (tmp_path / "sample.FSD.parquet").exists()

    @pytest.mark.unit
    def test_fsd_parquet_only_cleans_tsv(self, tmp_path, fsd_df):
        """process_fsd with parquet only removes intermediate TSV."""
        from krewlyzer.core.fsd_processor import _write_fsd_output

        tsv = tmp_path / "sample.FSD.tsv"
        fsd_df.to_csv(tsv, sep="\t", index=False)

        _write_fsd_output(tsv, output_format="parquet", compress=False)

        assert (tmp_path / "sample.FSD.parquet").exists()
        assert not tsv.exists(), "Intermediate TSV should be removed for parquet-only"


# ── OCF processor output format tests ────────────────────────────────────────


class TestOcfOutputFormat:
    """Test that OCF processor functions correctly produce requested formats."""

    @pytest.mark.unit
    def test_ocf_convert_both(self, tmp_path, ocf_df):
        """convert_ocf_output with 'both' produces .tsv + .parquet."""
        from krewlyzer.core.ocf_processor import convert_ocf_output

        tsv = tmp_path / "sample.OCF.tsv"
        ocf_df.to_csv(tsv, sep="\t", index=False)

        convert_ocf_output(tsv, output_format="both", compress=False)

        assert (tmp_path / "sample.OCF.tsv").exists()
        assert (tmp_path / "sample.OCF.parquet").exists()

    @pytest.mark.unit
    def test_ocf_convert_compress(self, tmp_path, ocf_df):
        """convert_ocf_output with compress produces .tsv.gz."""
        from krewlyzer.core.ocf_processor import convert_ocf_output

        tsv = tmp_path / "sample.OCF.tsv"
        ocf_df.to_csv(tsv, sep="\t", index=False)

        convert_ocf_output(tsv, output_format="tsv", compress=True)

        assert (tmp_path / "sample.OCF.tsv.gz").exists()

    @pytest.mark.unit
    def test_ocf_sync_both(self, tmp_path, ocf_df):
        """convert_ocf_output works on sync files too."""
        from krewlyzer.core.ocf_processor import convert_ocf_output

        tsv = tmp_path / "sample.OCF.sync.tsv"
        ocf_df.to_csv(tsv, sep="\t", index=False)

        convert_ocf_output(tsv, output_format="both", compress=True)

        assert (tmp_path / "sample.OCF.sync.tsv.gz").exists()
        assert (tmp_path / "sample.OCF.sync.parquet").exists()

    @pytest.mark.unit
    def test_ocf_parquet_only_cleans_tsv(self, tmp_path, ocf_df):
        """convert_ocf_output with parquet only removes intermediate TSV."""
        from krewlyzer.core.ocf_processor import _write_ocf_output

        tsv = tmp_path / "sample.OCF.tsv"
        ocf_df.to_csv(tsv, sep="\t", index=False)

        _write_ocf_output(tsv, output_format="parquet", compress=False)

        assert (tmp_path / "sample.OCF.parquet").exists()
        assert not tsv.exists(), "Intermediate TSV should be removed"


# ── Region entropy output format tests ───────────────────────────────────────


class TestRegionEntropyOutputFormat:
    """Test region entropy processor format handling in the with-PON branch."""

    @pytest.mark.unit
    def test_no_pon_both_format(self, tmp_path, entropy_df):
        """No-PON branch with output_format='both' produces tsv + parquet."""
        from krewlyzer.core.region_entropy_processor import process_region_entropy

        raw = tmp_path / "sample.TFBS.raw.tsv"
        out = tmp_path / "sample.TFBS.tsv"
        entropy_df.to_csv(raw, sep="\t", index=False)

        n = process_region_entropy(
            raw,
            out,
            pon_parquet_path=None,
            output_format="both",
            compress=False,
        )

        assert n == 0  # no PON, no z-scores
        assert (tmp_path / "sample.TFBS.tsv").exists()
        assert (tmp_path / "sample.TFBS.parquet").exists()

    @pytest.mark.unit
    def test_no_pon_compress(self, tmp_path, entropy_df):
        """No-PON branch with compress=True produces .tsv.gz."""
        from krewlyzer.core.region_entropy_processor import process_region_entropy

        raw = tmp_path / "sample.ATAC.raw.tsv"
        out = tmp_path / "sample.ATAC.tsv"
        entropy_df.to_csv(raw, sep="\t", index=False)

        process_region_entropy(
            raw,
            out,
            pon_parquet_path=None,
            output_format="tsv",
            compress=True,
        )

        assert (tmp_path / "sample.ATAC.tsv.gz").exists()


# ── fsc_counts format conversion tests ───────────────────────────────────────


class TestFscCountsFormat:
    """Test that fsc_counts respects output_format and compress."""

    @pytest.mark.unit
    def test_fsc_counts_both(self, tmp_path):
        """fsc_counts written as both TSV and Parquet."""
        df = pd.DataFrame({"bin": [100, 200], "count": [50, 30]})
        out = tmp_path / "sample.fsc_counts.tsv"
        df.to_csv(out, sep="\t", index=False)

        # Simulate format conversion (same as unified_processor.py line 523-535)
        loaded = read_table(out)
        assert loaded is not None
        write_table(loaded, out, output_format="both", compress=False)

        assert (tmp_path / "sample.fsc_counts.tsv").exists()
        assert (tmp_path / "sample.fsc_counts.parquet").exists()

    @pytest.mark.unit
    def test_fsc_counts_compress(self, tmp_path):
        """fsc_counts compressed to .tsv.gz."""
        df = pd.DataFrame({"bin": [100, 200], "count": [50, 30]})
        out = tmp_path / "sample.fsc_counts.tsv"
        df.to_csv(out, sep="\t", index=False)

        loaded = read_table(out)
        assert loaded is not None
        write_table(loaded, out, output_format="tsv", compress=True)

        assert (tmp_path / "sample.fsc_counts.tsv.gz").exists()
