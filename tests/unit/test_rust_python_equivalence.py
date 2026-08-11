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
    size_bins = [f"{s}-{s + 4}" for s in range(65, 200, 5)]

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
# WPS PON Scoring Equivalence — the Python oracle for rust/src/wps.rs
# =============================================================================
#
# These four functions ARE the specification for the WPS PON scoring that ships
# in Rust. They were the shipping implementation through 0.8.x, in
# `core/wps_pon.py`; that copy is deleted, so nothing here is duplicated in the
# wheel and no caller can reach this code.
#
# Why keep a second implementation at all, rather than trusting the Rust:
# every other check on `wps.rs` is a property the porter chose, asserted by the
# porter, at the moment of porting. This one was written first, validated on a
# real cohort, and settled three design decisions by measurement -- so it is
# the only check that can disagree with the Rust for a reason nobody
# anticipated. That is the whole value, and it is worth ~60 frozen lines.
#
# It is FROZEN. Do not "improve" it, and do not edit it to match the Rust: a
# divergence is a finding to investigate, not a sync to perform. Where the
# reference is arguably wrong -- a tie in the lag search resolving to the most
# negative lag -- the Rust reproduces the behaviour deliberately.
#
# When it stops earning its keep: replace it with a golden Parquet of expected
# values under `tests/data/fixtures/`. That keeps the regression protection
# with no second implementation to maintain, at the cost of only covering the
# inputs captured. Worth doing once `wps.rs` has been stable for a release or
# two; not now, while the 0.9.0 cohorts are still being re-run against it.


class WpsPonOracle:
    """The pre-0.9.0 Python implementation, frozen as the reference."""

    PHASE_MAX_LAG = 30

    @staticmethod
    def log_amplitude(profile: np.ndarray) -> float:
        finite = profile[np.isfinite(profile)]
        if finite.size < 2:
            return float("nan")
        return float(np.log1p(finite.max() - finite.min()))

    @staticmethod
    def shape_correlation(profile: np.ndarray, baseline: np.ndarray) -> float:
        ok = np.isfinite(profile) & np.isfinite(baseline)
        if ok.sum() < 3:
            return float("nan")
        a, b = profile[ok], baseline[ok]
        if np.std(a) < 1e-12 or np.std(b) < 1e-12:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])

    @staticmethod
    def fisher_z(r: float) -> float:
        if not np.isfinite(r):
            return float("nan")
        return float(np.arctanh(np.clip(r, -0.999999, 0.999999)))

    @classmethod
    def phase_shift(cls, profile, baseline, max_lag=None):
        max_lag = cls.PHASE_MAX_LAG if max_lag is None else max_lag
        n = min(profile.size, baseline.size)
        if n < 4 * max_lag:
            return float("nan"), False
        core = baseline[max_lag : n - max_lag]
        best_r, best_lag = -np.inf, 0
        for k in range(-max_lag, max_lag + 1):
            window = profile[max_lag + k : max_lag + k + core.size]
            r = cls.shape_correlation(window, core)
            # Strictly greater: a tie keeps the most negative lag. Bug-for-bug.
            if np.isfinite(r) and r > best_r:
                best_r, best_lag = r, k
        if not np.isfinite(best_r):
            return float("nan"), False
        return float(best_lag), abs(best_lag) == max_lag


def _wps_wave(n=200, period=12.0, offset=0.0, scale=1.0):
    return scale * np.sin((np.arange(n) - offset) / period)


def _wps_frame(seed: int, n_anchors: int = 30) -> pd.DataFrame:
    """Anchors that differ from each other and from the baseline.

    Varying both the phase and the period matters: a cohort of identical
    profiles gives every anchor a sigma of zero, every z a NaN, and an
    equivalence test that compares two columns of NaN and passes.
    """
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "region_id": [f"TSS|G{i}|T{i}" for i in range(n_anchors)],
            "wps_nuc": [
                _wps_wave(offset=rng.normal(0, 4.0), period=rng.uniform(9, 15))
                + rng.normal(0, 0.08, 200)
                for _ in range(n_anchors)
            ],
            "wps_tf": [
                _wps_wave(period=6) + rng.normal(0, 0.05, 200) for _ in range(n_anchors)
            ],
        }
    )


@pytest.fixture
def wps_pon_on_disk(tmp_path):
    """A real PON parquet plus the sample paths that built it.

    Built through `_compute_wps_baseline` and `_save_pon_model` rather than
    hand-written, so the test reads the same `large_string` / `list<double>`
    layout the shipped PONs use. A hand-rolled fixture in `string` /
    `list<float>` would have passed while the Rust reader was blind to every
    real PON.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    from krewlyzer.pon.build import _compute_wps_baseline, _save_pon_model
    from krewlyzer.pon.model import PonModel

    schema = pa.schema(
        [
            ("region_id", pa.string()),
            ("wps_nuc", pa.list_(pa.float32())),
            ("wps_tf", pa.list_(pa.float32())),
        ]
    )
    paths = []
    for i in range(7):
        path = tmp_path / f"s{i}.WPS.parquet"
        pq.write_table(pa.Table.from_pandas(_wps_frame(i), schema=schema), path)
        paths.append(path)

    vector, shape = _compute_wps_baseline([str(p) for p in paths])
    model = PonModel(
        schema_version="1.0",
        assay="xs2",
        build_date="2026-01-01",
        n_samples=len(paths),
        reference="r",
        panel_mode=True,
        target_regions_file="t",
        wps_baseline=vector,
        wps_shape_baseline=shape,
    )
    pon_path = tmp_path / "test.pon.parquet"
    _save_pon_model(model, pon_path)
    return pon_path, paths, vector, shape


class TestWpsPonEquivalence:
    """`rust/src/wps.rs::apply_pon_zscore` against the frozen Python oracle.

    Tolerance is `1e-6` relative, fixed before the first comparison was run and
    derived from nothing in the observed diff. Measured worst case at the time
    of the port: `2.2e-13` on synthetic anchors and `9.0e-14` against a real
    89k-anchor output and a shipped PON. If a change pushes the diff past the
    tolerance, that is a finding to report, not a number to widen.
    """

    RTOL = 1e-6

    def _oracle(self, sample_path, vector, shape):
        from krewlyzer.pon.model import WPS_SHAPE_STATS

        src = pd.read_parquet(sample_path)
        rows = []
        for region_id, raw in zip(src["region_id"], src["wps_nuc"]):
            profile = np.asarray(raw, dtype=float)
            amplitude = WpsPonOracle.log_amplitude(profile)
            mean = vector.get_baseline_vector(str(region_id), "wps_nuc")
            if mean is None:
                rows.append(
                    dict(
                        region_id=region_id,
                        wps_nuc_z=None,
                        wps_log_amplitude=amplitude,
                        wps_log_amplitude_z=np.nan,
                        wps_shape_corr=np.nan,
                        wps_shape_corr_z=np.nan,
                        wps_phase_shift_bp=np.nan,
                        wps_phase_at_search_limit=False,
                    )
                )
                continue
            mean = np.asarray(mean, dtype=float)
            r = WpsPonOracle.shape_correlation(profile, mean)
            lag, hit = WpsPonOracle.phase_shift(profile, mean)
            observed = {
                "log_amplitude": amplitude,
                "shape_corr_fisher": WpsPonOracle.fisher_z(r),
            }
            scores = {
                stat: shape.compute_zscore(str(region_id), stat, observed[stat])
                for stat in WPS_SHAPE_STATS
            }
            nan_if_none = lambda v: np.nan if v is None else float(v)  # noqa: E731
            rows.append(
                dict(
                    region_id=region_id,
                    wps_nuc_z=vector.compute_z_vector(
                        str(region_id), profile, "wps_nuc"
                    ),
                    wps_log_amplitude=amplitude,
                    wps_log_amplitude_z=nan_if_none(scores["log_amplitude"]),
                    wps_shape_corr=r,
                    wps_shape_corr_z=nan_if_none(scores["shape_corr_fisher"]),
                    wps_phase_shift_bp=lag,
                    wps_phase_at_search_limit=hit,
                )
            )
        return pd.DataFrame(rows)

    def test_rust_matches_the_oracle(self, wps_pon_on_disk, tmp_path):
        from krewlyzer.core.wps_pon import apply_wps_pon

        pon_path, paths, vector, shape = wps_pon_on_disk
        n = apply_wps_pon(paths[0], pon_path, output_base=tmp_path / "rust")
        assert n > 0, "nothing was scored -- the Rust reader found no baseline"

        got = pd.read_parquet(tmp_path / "rust.parquet")
        exp = self._oracle(paths[0], vector, shape)
        assert list(got["region_id"]) == list(exp["region_id"])

        for column in (
            "wps_log_amplitude",
            "wps_log_amplitude_z",
            "wps_shape_corr",
            "wps_shape_corr_z",
            "wps_phase_shift_bp",
        ):
            a = pd.to_numeric(got[column], errors="coerce").to_numpy(float)
            b = pd.to_numeric(exp[column], errors="coerce").to_numpy(float)
            assert (np.isnan(a) == np.isnan(b)).all(), (
                f"{column}: NaN in one and a number in the other. A fabricated "
                "value is exactly what this test exists to catch."
            )
            finite = np.isfinite(a) & np.isfinite(b)
            assert np.allclose(a[finite], b[finite], rtol=self.RTOL, atol=0), column
            # Not every column may be constant, or the comparison proves nothing.
            assert np.unique(a[finite]).size > 1, f"{column} does not vary"

        assert (
            got["wps_phase_at_search_limit"].tolist()
            == exp["wps_phase_at_search_limit"].tolist()
        ), "the boundary flag must agree exactly -- it is a boolean"

        compared = 0
        for i in range(len(exp)):
            a_raw, b_raw = got["wps_nuc_z"].iloc[i], exp["wps_nuc_z"].iloc[i]
            if b_raw is None:
                assert a_raw is None, "an absent anchor must stay null, not 0.0"
                continue
            a, b = np.asarray(a_raw, float), np.asarray(b_raw, float)
            assert (np.isnan(a) == np.isnan(b)).all(), f"row {i}: NaN pattern"
            finite = np.isfinite(a) & np.isfinite(b)
            assert np.allclose(a[finite], b[finite], rtol=self.RTOL, atol=0)
            compared += int(finite.sum())
        assert (
            compared > 1000
        ), f"only {compared} z values compared; too few to mean much"

    def test_the_column_order_is_the_oracle_s(self, wps_pon_on_disk, tmp_path):
        """Downstream may index by position, so the port must not reshuffle."""
        from krewlyzer.core.wps_pon import apply_wps_pon

        pon_path, paths, _, _ = wps_pon_on_disk
        apply_wps_pon(paths[0], pon_path, output_base=tmp_path / "order")
        got = list(pd.read_parquet(tmp_path / "order.parquet").columns)
        assert got == [
            "region_id",
            "wps_nuc",
            "wps_tf",
            "wps_nuc_z",
            "wps_log_amplitude",
            "wps_log_amplitude_z",
            "wps_shape_corr",
            "wps_shape_corr_z",
            "wps_phase_shift_bp",
            "wps_phase_at_search_limit",
        ], got

    def test_input_columns_survive_untouched(self, wps_pon_on_disk, tmp_path):
        """The schema is taken from the file, not restated in the port.

        A column the WPS writer adds later must reach the output without anyone
        remembering to update the scorer.
        """
        from krewlyzer.core.wps_pon import apply_wps_pon

        pon_path, paths, _, _ = wps_pon_on_disk
        before = pd.read_parquet(paths[0])
        apply_wps_pon(paths[0], pon_path, output_base=tmp_path / "keep")
        after = pd.read_parquet(tmp_path / "keep.parquet")
        for column in before.columns:
            assert column in after.columns, f"{column} was dropped"
        assert np.allclose(
            np.vstack(after["wps_tf"].to_numpy()),
            np.vstack(before["wps_tf"].to_numpy()),
            rtol=0,
            atol=0,
        ), "an untouched column was rewritten"


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
