"""When a PON baseline is missing, the scored column must be absent, not zero.

`zscore_or_nan` removed this fabrication from nine baseline classes. The apply
path kept five of its own, found by tracing every block from build to scoring:

    fsc_processor    {ch}_log2 = 0.0        when the PON has no gc_bias
    fsc_processor    {ch}_reliability = 1.0 when it has no variance
    fsr_processor    short_long_ratio = s   when there are no long fragments
    fsr_processor    short_long_log2 = 0.0  when the ratio is not positive
    region_entropy   baseline std = 0.0     from a single donor

Zero is never the cautious choice. A log2 ratio of zero says "this sample sits
exactly at the healthy baseline" and a z-score of zero says the same -- the most
confident statement either column can make, asserted precisely when there was
nothing to compare against. A reliability of 1.0 is not neutral either: with
``RELIABILITY_K`` at 0.01 it is what a variance of 0.99 would produce.

The entropy case shows why it matters most. A single-donor baseline stored
sigma 0.0, the Rust consumer treats sigma at or below 1e-9 as "cannot divide"
and emits 0.0, so every sample read as exactly average at that label -- while a
label *missing* from the baseline correctly read NaN. Two kinds of absence, and
the fabricating one fired more often.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from krewlyzer.pon.model import PonModel

RELIABILITY_K = 0.01


def _pon_without_gc_bias() -> PonModel:
    return PonModel(
        schema_version="1.0",
        assay="xs1",
        build_date="2026-08-04",
        n_samples=47,
        reference="hg19",
        panel_mode=True,
        target_regions_file="targets.bed",
    )


# ---------------------------------------------------------------------------
# the model accessors the FSC/FSR paths depend on
# ---------------------------------------------------------------------------


def test_a_pon_without_gc_bias_offers_no_mean_or_variance():
    """The premise for the two FSC branches below."""
    pon = _pon_without_gc_bias()
    assert pon.gc_bias is None
    assert pon.get_mean("short") is None
    assert pon.get_variance("short") is None


# ---------------------------------------------------------------------------
# FSC
# ---------------------------------------------------------------------------


def _fsc_frame() -> pd.DataFrame:
    """Windows whose raw values differ fourfold, so a constant output shows."""
    return pd.DataFrame(
        {
            "chrom": ["1", "1", "1"],
            "start": [0, 100_000, 200_000],
            "end": [100_000, 200_000, 300_000],
            "ultra_short": [10.0, 20.0, 40.0],
            "core_short": [100.0, 200.0, 400.0],
            "mono_nucl": [500.0, 1000.0, 2000.0],
            "di_nucl": [50.0, 100.0, 200.0],
            "long": [25.0, 50.0, 100.0],
            # `CHANNELS` has six entries -- omitting this one passed the
            # no-baseline test (that branch never indexes the column) and
            # failed only once a baseline was present.
            "ultra_long": [5.0, 10.0, 20.0],
        }
    )


def test_fsc_log2_is_absent_when_the_pon_has_no_gc_bias():
    """Not 0.0, which claims the sample sits exactly at the baseline.

    Drives the real `add_pon_channel_columns` rather than restating its
    branches -- a test that reimplements the code it checks agrees with itself
    forever, whatever the code goes on to do.
    """
    from krewlyzer.core.fsc_processor import CHANNELS, add_pon_channel_columns

    frame = add_pon_channel_columns(_fsc_frame(), _pon_without_gc_bias())
    for channel in CHANNELS:
        assert frame[f"{channel}_log2"].isna().all(), f"{channel}: fabricated log2"
        assert frame[f"{channel}_reliability"].isna().all(), f"{channel}: fabricated"


def test_fsc_log2_is_still_computed_when_the_baseline_is_there():
    """The fix must not turn a working column into NaN."""
    from krewlyzer.core.fsc_processor import add_pon_channel_columns

    class _Pon:
        def get_mean(self, channel):
            return 200.0

        def get_variance(self, channel):
            return 0.04

    frame = add_pon_channel_columns(_fsc_frame(), _Pon())
    log2 = frame["core_short_log2"]
    assert np.isfinite(log2).all()
    # ...and it varies with the data, which is the whole point of the column.
    assert log2.nunique() == 3
    assert np.isfinite(frame["core_short_reliability"]).all()


# ---------------------------------------------------------------------------
# FSR
# ---------------------------------------------------------------------------


def _short_long(s_norm: float, l_norm: float):
    """The two branches from `fsr_processor`, kept in step with the source."""
    ratio = s_norm / l_norm if l_norm > 0 else float("nan")
    log2 = float(np.log2(ratio)) if (np.isfinite(ratio) and ratio > 0) else float("nan")
    return ratio, log2


def test_no_long_fragments_means_no_ratio_not_a_ratio_of_the_numerator():
    """`short_long_ratio = s_norm` substituted the numerator for a divide by zero.

    A window with 40 short and 0 long read out as a ratio of 40 -- a specific,
    plausible, entirely invented number, which then flowed through log2 as
    though it had been measured.
    """
    ratio, log2 = _short_long(40.0, 0.0)
    assert math.isnan(ratio), f"invented a ratio of {ratio}"
    assert math.isnan(log2)


def test_an_empty_window_has_no_ratio():
    ratio, log2 = _short_long(0.0, 0.0)
    assert math.isnan(ratio) and math.isnan(log2)


def test_a_real_ratio_survives():
    ratio, log2 = _short_long(40.0, 10.0)
    assert ratio == pytest.approx(4.0)
    assert log2 == pytest.approx(2.0)


def test_the_fsr_source_no_longer_carries_the_zero_fallbacks():
    from krewlyzer.core import fsr_processor

    text = open((fsr_processor.__file__ or "").replace(".pyc", ".py")).read()
    assert "short_long_ratio = s_norm if s_norm > 0 else 0.0" not in text
    assert "short_long_log2 = 0.0" not in text


# ---------------------------------------------------------------------------
# region entropy -- the fifth sample_std_or_nan site
# ---------------------------------------------------------------------------


def test_a_single_donor_entropy_baseline_has_no_spread():
    """Not 0.0, which the Rust consumer turns into a z-score of exactly 0.0."""
    from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline

    baseline = RegionEntropyBaseline.from_samples([{"CTCF": 0.71}])
    mean, std = baseline.data["CTCF"]
    assert mean == pytest.approx(0.71), "the centre is still measurable"
    assert math.isnan(std), f"expected NaN, got {std}"


def test_a_real_entropy_cohort_still_gets_a_spread():
    """So the fix cannot be 'return NaN always'."""
    from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline

    values = [0.71, 0.74, 0.69]
    baseline = RegionEntropyBaseline.from_samples([{"CTCF": v} for v in values])
    mean, std = baseline.data["CTCF"]
    assert mean == pytest.approx(np.mean(values))
    assert std == pytest.approx(np.std(values, ddof=1))


def test_why_zero_was_the_worst_answer_here():
    """A NaN sigma cannot be divided by. A zero one gets treated as 'skip'.

    The Rust branch reads `if entropy_std > 1e-9 { z } else { 0.0 }`, so a
    stored zero produced a confident z of exactly 0.0 for every sample at that
    label.
    """
    observed, mean = 0.85, 0.71
    assert math.isnan((observed - mean) / float("nan"))
    # The old stored sigma took the other branch and reported this sample --
    # 0.14 above the baseline -- as sitting exactly on it.
    old_sigma = 0.0
    assert not old_sigma > 1e-9, "the Rust guard would have emitted 0.0 here"
