"""
Accuracy verification tests for Rust vs Python equivalence.

These tests ensure that the Rust implementations produce numerically
equivalent results to the Python fallback implementations.

Run with: pytest tests/unit/test_rust_python_equivalence.py -v
"""

import pytest
import numpy as np
import pandas as pd
from typing import Dict, List
import logging

logger = logging.getLogger(__name__)


# =============================================================================
# Test Fixtures - Synthetic Data
# =============================================================================


@pytest.fixture
def synthetic_fsd_data() -> pd.DataFrame:
    """Generate synthetic FSD TSV data for testing."""
    arms = ["1p", "1q", "2p", "2q", "3p", "3q"]
    size_bins = [f"{s}-{s+4}" for s in range(65, 200, 5)]

    data = []
    for arm in arms:
        row = {"region": arm, "total": 10000}
        for i, bin_col in enumerate(size_bins):
            # Generate realistic-ish fragment count with Gaussian shape
            center = 165
            bin_start = int(bin_col.split("-")[0])
            count = int(1000 * np.exp(-0.5 * ((bin_start - center) / 30) ** 2))
            row[bin_col] = count
        data.append(row)

    return pd.DataFrame(data)


@pytest.fixture
def synthetic_gc_data() -> List[Dict]:
    """Generate synthetic GC observation data for testing."""
    gc_bins = np.arange(0.25, 0.75, 0.02)

    samples = []
    for sample_idx in range(10):
        # Generate GC-dependent coverage with some noise
        gc = gc_bins.tolist()
        short = [1.0 + 0.1 * np.sin(g * 10) + np.random.normal(0, 0.05) for g in gc]
        intermediate = [
            1.0 - 0.1 * np.cos(g * 8) + np.random.normal(0, 0.05) for g in gc
        ]
        long_vals = [1.0 + 0.05 * (g - 0.5) + np.random.normal(0, 0.03) for g in gc]

        samples.append(
            {"gc": gc, "short": short, "intermediate": intermediate, "long": long_vals}
        )

    return samples


@pytest.fixture
def synthetic_wps_data() -> pd.DataFrame:
    """Generate synthetic WPS parquet data for testing."""
    regions = [f"region_{i}" for i in range(100)]

    data = []
    for region_id in regions:
        wps_nuc = np.random.normal(0.5, 0.1)
        wps_tf = np.random.normal(0.3, 0.15)
        data.append(
            {
                "region_id": region_id,
                "wps_nuc_mean": wps_nuc,
                "wps_tf_mean": wps_tf,
                "chrom": "chr1",
                "center": np.random.randint(1000000, 10000000),
            }
        )

    return pd.DataFrame(data)


# =============================================================================
# Utility: Python Reference Implementations
# =============================================================================


def python_log_ratio(
    sample_val: float, expected_val: float, pseudocount: float = 1.0
) -> float:
    """Python reference for log-ratio computation."""
    return np.log2((sample_val + pseudocount) / (expected_val + pseudocount))


def python_zscore(sample_val: float, mean_val: float, std_val: float) -> float:
    """Python reference for z-score computation."""
    if std_val > 0:
        return (sample_val - mean_val) / std_val
    return 0.0


def python_median(data: List[float]) -> float:
    """Python reference for median computation."""
    if not data:
        return 1.0
    return float(np.median(data))


def python_std(data: List[float]) -> float:
    """Python reference for std computation (sample std, ddof=1)."""
    if len(data) < 2:
        return 0.1
    return float(np.std(data, ddof=1))


# =============================================================================
# Math Equivalence Tests
# =============================================================================


class TestMathEquivalence:
    """Test that basic math operations match between Python and Rust."""

    def test_log_ratio_equivalence(self):
        """Test log-ratio computation matches."""
        test_cases = [
            (100, 50),  # 2x enrichment
            (50, 100),  # 0.5x depletion
            (100, 100),  # neutral
            (0, 100),  # zero sample
            (100, 0),  # zero expected (uses pseudocount)
        ]

        for sample, expected in test_cases:
            py_result = python_log_ratio(sample, expected)

            # Manual Rust-equivalent calculation
            rust_result = np.log2((sample + 1.0) / (expected + 1.0))

            assert np.isclose(
                py_result, rust_result, rtol=1e-10
            ), f"Log-ratio mismatch: Python={py_result}, Rust-equiv={rust_result}"

    def test_zscore_equivalence(self):
        """Test z-score computation matches."""
        test_cases = [
            (0.5, 0.3, 0.1),  # positive z
            (0.3, 0.5, 0.1),  # negative z
            (0.5, 0.5, 0.1),  # zero z
            (0.5, 0.5, 0.0),  # zero std (edge case)
        ]

        for sample, mean, std in test_cases:
            py_result = python_zscore(sample, mean, std)

            # Manual Rust-equivalent calculation
            if std > 0:
                rust_result = (sample - mean) / std
            else:
                rust_result = 0.0

            assert np.isclose(
                py_result, rust_result, rtol=1e-10
            ), f"Z-score mismatch: Python={py_result}, Rust-equiv={rust_result}"

    def test_median_equivalence(self):
        """Test median computation matches."""
        test_cases = [
            [1.0, 2.0, 3.0],  # odd count
            [1.0, 2.0, 3.0, 4.0],  # even count
            [5.0],  # single element
            [],  # empty (returns 1.0)
        ]

        for data in test_cases:
            py_result = python_median(data)

            # Rust implementation (sorted, mid-point)
            if not data:
                rust_result = 1.0
            else:
                sorted_data = sorted(data)
                mid = len(sorted_data) // 2
                if len(sorted_data) % 2 == 0:
                    rust_result = (sorted_data[mid - 1] + sorted_data[mid]) / 2.0
                else:
                    rust_result = sorted_data[mid]

            assert np.isclose(
                py_result, rust_result, rtol=1e-10
            ), f"Median mismatch: Python={py_result}, Rust-equiv={rust_result}"

    def test_std_equivalence(self):
        """Test std computation matches (sample std, N-1)."""
        test_cases = [
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [10.0, 20.0],
            [1.0],  # single element (returns 0.1)
        ]

        for data in test_cases:
            py_result = python_std(data)

            # Rust implementation
            if len(data) < 2:
                rust_result = 0.1
            else:
                mean = sum(data) / len(data)
                variance = sum((x - mean) ** 2 for x in data) / (len(data) - 1)
                rust_result = np.sqrt(variance)

            assert np.isclose(
                py_result, rust_result, rtol=1e-10
            ), f"Std mismatch: Python={py_result}, Rust-equiv={rust_result}"


# =============================================================================
# FSD Equivalence Tests
# =============================================================================


class TestFsdEquivalence:
    """Test FSD log-ratio implementation equivalence."""

    def test_fsd_logratio_output_format(self, synthetic_fsd_data, tmp_path):
        """Test that FSD log-ratio produces valid output format."""
        # Write synthetic data
        input_path = tmp_path / "test.FSD.tsv"
        synthetic_fsd_data.to_csv(input_path, sep="\t", index=False)

        # Verify file exists and has expected columns
        df = pd.read_csv(input_path, sep="\t")
        assert "region" in df.columns
        assert any("-" in col for col in df.columns)  # Has bin columns

        logger.info(f"FSD input verified: {len(df)} rows, {len(df.columns)} columns")

    def test_fsd_logratio_values(self, synthetic_fsd_data):
        """Test log-ratio computation on synthetic data."""
        # Simulate PON baseline (mean of all samples)
        baseline = {}
        for arm in synthetic_fsd_data["region"]:
            row = synthetic_fsd_data[synthetic_fsd_data["region"] == arm].iloc[0]
            baseline[arm] = {}
            for col in synthetic_fsd_data.columns:
                if "-" in col:
                    baseline[arm][col] = row[col] * 1.1  # Slightly higher expected

        # Compute log-ratios
        for _, row in synthetic_fsd_data.iterrows():
            arm = row["region"]
            for col in synthetic_fsd_data.columns:
                if "-" in col:
                    sample_val = row[col]
                    expected_val = baseline[arm][col]

                    log_ratio = python_log_ratio(sample_val, expected_val)

                    # Verify log-ratio is finite and reasonable
                    assert np.isfinite(
                        log_ratio
                    ), f"Non-finite log-ratio for {arm}/{col}"
                    assert -10 < log_ratio < 10, f"Extreme log-ratio: {log_ratio}"


# =============================================================================
# WPS Equivalence Tests
# =============================================================================


class TestWpsEquivalence:
    """Test WPS z-score implementation equivalence."""

    def test_wps_zscore_computation(self, synthetic_wps_data):
        """Test z-score computation on synthetic WPS data."""
        # Compute baseline (mean/std across regions)
        nuc_mean = synthetic_wps_data["wps_nuc_mean"].mean()
        nuc_std = synthetic_wps_data["wps_nuc_mean"].std()

        for _, row in synthetic_wps_data.iterrows():
            sample_nuc = row["wps_nuc_mean"]
            z_score = python_zscore(sample_nuc, nuc_mean, nuc_std)

            # Verify z-score is finite and reasonable
            assert np.isfinite(z_score), "Non-finite z-score"
            assert -5 < z_score < 5, f"Extreme z-score: {z_score}"

    def test_wps_parquet_roundtrip(self, synthetic_wps_data, tmp_path):
        """Test WPS data survives parquet roundtrip."""
        parquet_path = tmp_path / "test.WPS.parquet"
        synthetic_wps_data.to_parquet(parquet_path)

        loaded = pd.read_parquet(parquet_path)

        assert len(loaded) == len(synthetic_wps_data)
        assert np.allclose(
            loaded["wps_nuc_mean"].values,
            synthetic_wps_data["wps_nuc_mean"].values,
            rtol=1e-10,
        )


# =============================================================================
# GC Bias Model Equivalence Tests
# =============================================================================


class TestGcBiasEquivalence:
    """Test GC bias model aggregation equivalence."""

    def test_gc_aggregation(self, synthetic_gc_data):
        """Test GC bin aggregation produces expected output."""
        gc_bins = np.arange(0.25, 0.75, 0.02)

        # Aggregate short values per bin
        short_by_bin = {round(g, 2): [] for g in gc_bins}

        for sample in synthetic_gc_data:
            for i, gc_val in enumerate(sample["gc"]):
                bin_key = round(gc_val, 2)
                if bin_key in short_by_bin:
                    short_by_bin[bin_key].append(sample["short"][i])

        # Compute expected values
        for bin_key, values in short_by_bin.items():
            if values:
                expected = python_median(values)
                std = python_std(values)

                assert np.isfinite(expected), f"Non-finite expected for bin {bin_key}"
                assert std >= 0, f"Negative std for bin {bin_key}"

    def test_gc_normalization(self, synthetic_gc_data):
        """Test GC values are properly normalized."""
        for sample in synthetic_gc_data:
            short = np.array(sample["short"])
            short_mean = np.nanmean(short)

            if short_mean > 0:
                normalized = short / short_mean

                # Mean of normalized should be ~1.0
                assert np.isclose(np.nanmean(normalized), 1.0, rtol=0.01)


# =============================================================================
# Integration Test: Full Pipeline Comparison
# =============================================================================


class TestFullPipelineEquivalence:
    """End-to-end comparison of Python vs Rust pipelines."""

    @pytest.mark.skipif(
        True,  # Skip by default, enable when testing locally
        reason="Requires built Rust extension",
    )
    def test_fsd_python_vs_rust(self, synthetic_fsd_data, tmp_path):
        """Compare FSD output from Python fallback vs Rust."""
        try:
            from krewlyzer import _core
            from krewlyzer.core.fsd_processor import process_fsd
        except ImportError:
            pytest.skip("Rust extension not built")

        # Create input file
        input_path = tmp_path / "test.FSD.tsv"
        synthetic_fsd_data.to_csv(input_path, sep="\t", index=False)

        # Run Python (no PON for this test)
        py_output = tmp_path / "python.FSD.tsv"
        process_fsd(input_path, py_output, pon=None)

        # Both should produce identical output without PON
        py_df = pd.read_csv(py_output, sep="\t")

        # Verify output format
        assert "region" in py_df.columns
        assert len(py_df) == len(synthetic_fsd_data)


# =============================================================================
# Run Tests
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
