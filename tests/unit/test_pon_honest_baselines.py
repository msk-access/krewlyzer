"""A PON baseline must be measured, or absent. Never fabricated.

Every defect found in the shipped models is the same shape: a value that is
present, plausible, and never derived from data.

``wps_background`` was the worst of them. The builder looked for columns named
``nrl`` and ``period_score``; ``WPS_background`` writes ``nrl_bp`` and
``periodicity_score``. Neither ever matched, so the fallback fired on every
build and all four shipped PONs carry an identical hardcoded
``167.0 / 5.0 / 0.0 / 1.0`` across all 28 groups — from cohorts of 21 and 47
samples, logging ``28 groups`` as though it had fitted them. Combined with the
0.8.x ``nrl_bp`` degeneracy, every sample ever produced scored
``nrl_z = (150 - 167) / 5 = -3.4``: a constant, moderately extreme z-score
presented as a measurement.

The other half is the σ floors. ``max(std, 0.001)`` does not make a z-score
conservative — it makes it arbitrarily large, because the divisor is a number
nothing measured. That is ``nrl_at_band_limit`` one level down: a boundary
value reported as a result. The honest output is NaN, which propagates to an
absent z.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from krewlyzer.pon.build import (
    _compute_mds_baseline,
    _compute_wps_background_baseline,
)


def _write(directory: Path, name: str, **columns) -> str:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}.parquet"
    pd.DataFrame(columns).to_parquet(path)
    return str(path)


# ---------------------------------------------------------------------------
# wps_background
# ---------------------------------------------------------------------------


def test_the_columns_the_builder_used_to_guess_are_refused(tmp_path):
    """The exact shape that produced four PONs of hardcoded constants.

    A missing source column must stop the build. Substituting a literal is
    worse than failing: nothing downstream can distinguish the result from a
    fitted one, and no schema check will ever flag it.
    """
    path = _write(
        tmp_path, "s", group_id=["Global_All"], nrl=[190.0], periodicity=[0.4]
    )
    with pytest.raises(ValueError, match="nrl_bp"):
        _compute_wps_background_baseline([path])


def test_the_refusal_names_what_it_found(tmp_path):
    """So the fix is obvious from the error, not from reading the builder."""
    path = _write(tmp_path, "s", group_id=["G"], nrl_bp=[190.0])
    with pytest.raises(ValueError) as excinfo:
        _compute_wps_background_baseline([path])
    message = str(excinfo.value)
    assert "periodicity_score" in message, "does not say which column is missing"
    assert "nrl_bp" in message, "does not list what it did find"


def test_a_real_baseline_is_fitted_from_the_data(tmp_path):
    """The values must track the input, which the hardcoded version never did."""
    paths = [
        _write(
            tmp_path,
            f"s{i}",
            group_id=["Global_All"],
            nrl_bp=[nrl],
            periodicity_score=[per],
        )
        for i, (nrl, per) in enumerate([(185.0, 0.41), (190.0, 0.44), (196.0, 0.39)])
    ]
    row = _compute_wps_background_baseline(paths).groups.iloc[0]
    assert row.nrl_mean == pytest.approx(190.333, abs=1e-3)
    assert row.nrl_std == pytest.approx(np.std([185.0, 190.0, 196.0], ddof=1))
    assert row.n_samples == 3
    assert row.nrl_mean != 167.0 and row.nrl_std != 5.0


def test_two_different_cohorts_give_two_different_baselines(tmp_path):
    """The anti-degeneracy invariant, applied to the PON itself.

    This is the single assertion that would have caught the placeholder: the
    shipped baseline was byte-identical across four models built from two
    different cohorts.
    """

    def fit(values):
        d = tmp_path / f"c{hash(values) & 0xFFFF}"
        paths = [
            _write(
                d,
                f"s{i}",
                group_id=["G"],
                nrl_bp=[v],
                periodicity_score=[0.4 + i / 100],
            )
            for i, v in enumerate(values)
        ]
        return _compute_wps_background_baseline(paths).groups.iloc[0]

    a, b = fit((180.0, 185.0, 190.0)), fit((195.0, 200.0, 205.0))
    assert a.nrl_mean != b.nrl_mean
    assert not np.isclose(a.nrl_mean, b.nrl_mean)


def test_an_unmeasurable_spread_is_absent_not_floored(tmp_path):
    """One sample has no spread. 0.1 is not a smaller answer, it is a wrong one."""
    path = _write(
        tmp_path, "only", group_id=["G"], nrl_bp=[190.0], periodicity_score=[0.4]
    )
    row = _compute_wps_background_baseline([path]).groups.iloc[0]
    assert row.nrl_mean == pytest.approx(190.0), "the mean is still measurable"
    assert np.isnan(row.nrl_std), f"expected NaN, got {row.nrl_std}"


def test_a_nan_sigma_produces_no_z_rather_than_an_enormous_one():
    """The reason NaN is the right answer: it propagates to absence.

    A floor of 0.01 would have turned this same input into z = 500.
    """
    observed, mean, sigma = 195.0, 190.0, float("nan")
    assert np.isnan((observed - mean) / sigma)
    assert (observed - mean) / 0.01 == pytest.approx(500.0)


# ---------------------------------------------------------------------------
# scalar baselines
# ---------------------------------------------------------------------------


def test_an_unfittable_mds_baseline_is_nan_not_a_standard_normal():
    """`mean=0.0, std=1.0` makes z equal the raw MDS value.

    MDS runs ~0.95, so the fabricated z is ~0.95 — an entirely unremarkable
    number, which is precisely why it would never have been questioned.
    """
    baseline = _compute_mds_baseline([{"mds": None, "kmers": {}}])
    assert np.isnan(baseline.mds_mean)
    assert np.isnan(baseline.mds_std)


def test_the_mds_baseline_uses_the_sample_standard_deviation():
    """ddof=1: these are a sample of donors, not the population.

    np.std defaults to ddof=0, which understates the spread and inflates every
    z built from it — by 2.5% at n=21.
    """
    values = [0.95, 0.96, 0.97]
    baseline = _compute_mds_baseline([{"mds": v, "kmers": {}} for v in values])
    assert baseline.mds_std == pytest.approx(np.std(values, ddof=1))
    assert baseline.mds_std != pytest.approx(np.std(values))
