"""
Unit tests for FSR (Fragment Size Ratio) processor.

Tests the Python-side FSR processing: bin aggregation, region label
computation, ratio calculation, and PoN normalization order.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from krewlyzer.core.fsr_processor import process_fsr


def _make_counts_df(
    n_bins: int = 50,
    chrom: str = "chr1",
    bin_size: int = 100000,
) -> pd.DataFrame:
    """Create a synthetic fsc_counts DataFrame with realistic structure.

    Args:
        n_bins: Number of bins to generate
        chrom: Chromosome name
        bin_size: Size of each bin in bp (default 100kb)

    Returns:
        DataFrame matching the fsc_counts schema
    """
    starts = np.arange(n_bins) * bin_size
    ends = starts + bin_size
    rng = np.random.default_rng(42)

    return pd.DataFrame(
        {
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "ultra_short": rng.integers(10, 100, n_bins),
            "core_short": rng.integers(10, 100, n_bins),
            "mono_nucl": rng.integers(50, 200, n_bins),
            "di_nucl": rng.integers(5, 50, n_bins),
            "long": rng.integers(5, 50, n_bins),
            "total": rng.integers(100, 400, n_bins),
            "mean_gc": rng.uniform(0.3, 0.6, n_bins),
        }
    )


# =========================================================================
# Region Label Tests — the core fix
# =========================================================================


@pytest.mark.unit
def test_wgs_mode_labels(tmp_path):
    """WGS mode: 50 bins × 100kb → 1 window with correct 5Mb label."""
    df = _make_counts_df(n_bins=50, bin_size=100000)
    output = tmp_path / "test.FSR.tsv"

    process_fsr(df, output, windows=100000, continue_n=50)

    result = pd.read_csv(output, sep="\t")
    assert len(result) == 1, "50 bins / continue_n=50 should give 1 window"

    # Label should reflect actual genomic coordinates: chr1:0-5000000
    region = result["region"].iloc[0]
    assert region == "chr1:0-5000000", f"Expected chr1:0-5000000, got {region}"


@pytest.mark.unit
def test_panel_mode_labels(tmp_path):
    """Panel mode: continue_n=1 → per-bin rows with correct 100kb labels.

    This is the primary regression test for the bug where panel mode
    produced labels like 'chr1:55-55' instead of 'chr1:5500000-5600000'.
    """
    df = _make_counts_df(n_bins=5, bin_size=100000)
    output = tmp_path / "test.FSR.tsv"

    process_fsr(df, output, windows=100000, continue_n=1)

    result = pd.read_csv(output, sep="\t")
    assert len(result) == 5, "5 bins / continue_n=1 should give 5 rows"

    # Each row should have correct genomic coordinates
    expected_regions = [
        "chr1:0-100000",
        "chr1:100000-200000",
        "chr1:200000-300000",
        "chr1:300000-400000",
        "chr1:400000-500000",
    ]
    assert result["region"].tolist() == expected_regions


@pytest.mark.unit
def test_custom_bin_sizes(tmp_path):
    """Non-100kb bins produce correct labels from actual DataFrame coords."""
    # Use 200kb bins — label should come from data, not formula
    df = _make_counts_df(n_bins=4, bin_size=200000)
    output = tmp_path / "test.FSR.tsv"

    # Even with windows=100000 (wrong!), labels should be correct
    # because FSR now reads real coordinates from the DataFrame
    process_fsr(df, output, windows=100000, continue_n=2)

    result = pd.read_csv(output, sep="\t")
    assert len(result) == 2, "4 bins / continue_n=2 should give 2 windows"

    # Labels reflect actual 200kb bins aggregated (0-400000, 400000-800000)
    assert result["region"].iloc[0] == "chr1:0-400000"
    assert result["region"].iloc[1] == "chr1:400000-800000"


# =========================================================================
# Output Column Tests
# =========================================================================


@pytest.mark.unit
def test_output_columns(tmp_path):
    """FSR output should have all 10 expected columns."""
    df = _make_counts_df(n_bins=50)
    output = tmp_path / "test.FSR.tsv"

    process_fsr(df, output, continue_n=50)

    result = pd.read_csv(output, sep="\t")
    expected_cols = [
        "region",
        "short_count",
        "long_count",
        "total_count",
        "short_norm",
        "long_norm",
        "short_long_ratio",
        "short_long_log2",
        "short_frac",
        "long_frac",
    ]
    assert list(result.columns) == expected_cols


# =========================================================================
# Edge Cases
# =========================================================================


@pytest.mark.unit
def test_empty_input(tmp_path):
    """Empty DataFrame → empty output with headers."""
    df = pd.DataFrame(
        columns=[
            "chrom",
            "start",
            "end",
            "ultra_short",
            "core_short",
            "mono_nucl",
            "di_nucl",
            "long",
            "total",
        ]
    )
    output = tmp_path / "test.FSR.tsv"

    result_path = process_fsr(df, output)

    assert result_path.exists()
    result = pd.read_csv(result_path, sep="\t")
    assert len(result) == 0
    assert "region" in result.columns


@pytest.mark.unit
def test_multi_chrom(tmp_path):
    """Multiple chromosomes are each aggregated independently."""
    df1 = _make_counts_df(n_bins=50, chrom="chr1")
    df2 = _make_counts_df(n_bins=50, chrom="chr2")
    df = pd.concat([df1, df2], ignore_index=True)
    output = tmp_path / "test.FSR.tsv"

    process_fsr(df, output, continue_n=50)

    result = pd.read_csv(output, sep="\t")
    assert len(result) == 2, "50 bins per chrom / continue_n=50 = 1 window each"

    regions = result["region"].tolist()
    assert regions[0].startswith("chr1:")
    assert regions[1].startswith("chr2:")


@pytest.mark.unit
def test_missing_required_columns():
    """Missing required columns raise ValueError."""
    df = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [100000]})

    with pytest.raises(ValueError, match="Missing required columns"):
        process_fsr(df, Path("/tmp/dummy.tsv"))
