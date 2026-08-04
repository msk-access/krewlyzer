"""No spread is not a small spread — it is no measurement.

`pandas.std()` and `np.std(ddof=1)` both return **0.0** for identical values,
and NaN only below two observations. 0.0 is the single worst answer available:
a z-score divided by it is infinite rather than absent.

Six places in `pon/build.py` computed a spread and **one** converted zero to
NaN. Three of the other five carried comments asserting they did. The models
built from them failed their own `validate-pon` gate:

    xs1/duplex      fsc_region_baseline.depth_std   2 non-positive
                    region_mds_exon.mds_std          2 non-positive
                    wps_background.nrl_std           1 non-positive
    xs2/all_unique  wps_background.nrl_std           3 non-positive
    xs2/duplex      wps_background.nrl_std           4 non-positive

So the rule lives in one function now, and every builder is pinned to it
individually — a comment is not a test, which is exactly what those three
proved.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from krewlyzer.pon.build import (
    MIN_SAMPLES_PER_KEY,
    _compute_fsc_gene_baseline,
    _compute_fsc_region_baseline,
    _compute_mds_baseline,
    _compute_region_mds_exon_baseline,
    _compute_wps_background_baseline,
    sample_std_or_nan,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# the helper
# ---------------------------------------------------------------------------


def test_a_real_spread_uses_the_sample_correction():
    """ddof=1: a cohort is a sample of donors, not the population.

    The population form understates the spread and inflates every z built from
    it — by 2.5% at n=21.
    """
    assert sample_std_or_nan([2.0, 4.0, 6.0]) == pytest.approx(2.0)
    assert sample_std_or_nan([2.0, 4.0, 6.0]) != pytest.approx(np.std([2.0, 4.0, 6.0]))


@pytest.mark.parametrize(
    "values,why",
    [
        ([3.0, 3.0, 3.0], "identical values measure no spread"),
        ([5.0], "one observation has none to measure"),
        ([], "nothing at all"),
        ([float("nan"), float("nan")], "nothing finite"),
    ],
)
def test_no_measurable_spread_is_nan(values, why):
    assert math.isnan(sample_std_or_nan(values)), why


def test_a_tiny_but_real_spread_survives():
    """The counterpart: no floor sneaks in.

    `region_mds` sigmas run ~0.005 and whole-sample MDS ~0.0008. A guard that
    rejected "small" would silence the most precise baselines in the model.
    """
    sd = sample_std_or_nan([1.0, 1.000001])
    assert 0 < sd < 1e-5


def test_non_finite_values_are_dropped_not_propagated():
    """One NaN sample must not erase a group that has real observations."""
    assert sample_std_or_nan([1.0, float("nan"), 3.0]) == pytest.approx(
        np.std([1.0, 3.0], ddof=1)
    )


# ---------------------------------------------------------------------------
# every builder, individually
# ---------------------------------------------------------------------------


def _exon_paths(tmp_path, value):
    paths = []
    for i in range(MIN_SAMPLES_PER_KEY):
        p = tmp_path / f"s{i}.MDS.exon.tsv"
        pd.DataFrame({"gene": ["TP53"], "name": ["ex1"], "mds": [value]}).to_csv(
            p, sep="\t", index=False
        )
        paths.append(str(p))
    return paths


def _wps_bg_paths(tmp_path, value):
    paths = []
    for i in range(MIN_SAMPLES_PER_KEY):
        p = tmp_path / f"s{i}.WPS_background.parquet"
        pd.DataFrame(
            {
                "group_id": ["Global_All"],
                "nrl_bp": [value],
                "periodicity_score": [0.4],
            }
        ).to_parquet(p)
        paths.append(str(p))
    return paths


def test_region_mds_exon_gives_nan_for_identical_values(tmp_path):
    """Its comment claimed NaN; `float(row["std"])` returned 0.0."""
    baseline = _compute_region_mds_exon_baseline(_exon_paths(tmp_path, 0.9))
    assert math.isnan(baseline.exon_baseline[("TP53", "ex1")]["mds_std"])


def test_wps_background_gives_nan_for_identical_values(tmp_path):
    """The one that actually reached three shipped models."""
    row = _compute_wps_background_baseline(_wps_bg_paths(tmp_path, 190.0)).groups.iloc[
        0
    ]
    assert math.isnan(row.nrl_std)
    assert math.isnan(row.periodicity_std), "periodicity has the same problem"
    assert row.nrl_mean == pytest.approx(190.0), "the mean is still measurable"


def test_fsc_gene_gives_nan_for_identical_values():
    baseline = _compute_fsc_gene_baseline([{"TP53": 100.0} for _ in range(3)])
    assert math.isnan(baseline.data["TP53"][1])


def test_fsc_region_gives_nan_for_identical_values():
    """Its comment claimed NaN; the value went through unchanged."""
    baseline = _compute_fsc_region_baseline([{"1:100-200": 50.0} for _ in range(3)])
    assert math.isnan(baseline.data["1:100-200"][1])


def test_the_mds_baselines_give_nan_for_identical_values():
    """Both the whole-sample score and the 256 k-mer frequencies."""
    samples = [{"mds": 0.95, "kmers": {"ACGT": 0.25}} for _ in range(3)]
    baseline = _compute_mds_baseline(samples)
    assert math.isnan(baseline.mds_std)
    assert math.isnan(baseline.kmer_std["ACGT"])


def test_a_varying_cohort_still_gets_a_real_sigma(tmp_path):
    """The other direction, so the fix cannot be 'return NaN always'."""
    paths = []
    for i, v in enumerate((185.0, 190.0, 196.0)):
        p = tmp_path / f"v{i}.WPS_background.parquet"
        pd.DataFrame(
            {"group_id": ["Global_All"], "nrl_bp": [v], "periodicity_score": [0.4 + i]}
        ).to_parquet(p)
        paths.append(str(p))
    row = _compute_wps_background_baseline(paths).groups.iloc[0]
    assert row.nrl_std == pytest.approx(np.std([185.0, 190.0, 196.0], ddof=1))


def test_identity_is_tested_on_the_values_not_the_result():
    """`std([0.95, 0.95, 0.95])` is 1.36e-16, not 0.0.

    Cancellation in the variance sum leaves floating-point residue, so a
    `sd <= 0` guard misses the exact case it exists for — and a z divided by
    1.36e-16 is 1e16, which is worse than the 0.0 it was meant to prevent.
    Found by this file's own builder tests before it could reach a model.
    """
    naive = float(pd.Series([0.95, 0.95, 0.95]).std(ddof=1))
    assert naive > 0.0, "the premise: pandas does not return exactly zero here"
    assert naive < 1e-15, f"and it is residue, not signal: {naive}"
    assert math.isnan(sample_std_or_nan([0.95, 0.95, 0.95]))


def test_values_that_genuinely_differ_keep_their_spread():
    """The identity test must not become a tolerance in disguise.

    Two values differing by one ULP have a real, if tiny, spread. Rejecting it
    would be the floor this whole effort removes, reintroduced as a numerical
    convenience.
    """
    a = 0.95
    b = math.nextafter(a, 1.0)
    assert b != a
    assert not math.isnan(sample_std_or_nan([a, b]))
