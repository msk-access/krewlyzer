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
            nrl_at_band_limit=[False],
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
                nrl_at_band_limit=[False],
            )
            for i, v in enumerate(values)
        ]
        return _compute_wps_background_baseline(paths).groups.iloc[0]

    a, b = fit((180.0, 185.0, 190.0)), fit((195.0, 200.0, 205.0))
    assert a.nrl_mean != b.nrl_mean
    assert not np.isclose(a.nrl_mean, b.nrl_mean)


def test_an_unmeasurable_spread_is_absent_not_floored(tmp_path):
    """Identical donors have no spread. 0.1 is not a smaller answer, it is a wrong one.

    Uses `MIN_SAMPLES_PER_KEY` donors deliberately: this test is about the
    *spread* being unmeasurable while the centre is not, so the cohort has to
    clear the sample floor. One donor is a different failure, pinned by
    `test_too_few_measured_rows_drops_the_mean_too_not_just_the_spread`.
    """
    from krewlyzer.pon.build import MIN_SAMPLES_PER_KEY

    paths = [
        _bg(tmp_path, f"same{i}", [("G", 190.0, 0.4, False)])
        for i in range(MIN_SAMPLES_PER_KEY)
    ]
    row = _compute_wps_background_baseline(paths).groups.iloc[0]
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


def test_an_empty_wps_baseline_says_why(tmp_path):
    """The message has to name the sample floor, not blame the backend.

    A cohort below `MIN_SAMPLES_PER_KEY` has every anchor dropped, so the Rust
    call correctly returns nothing — and the old wording, "no data returned
    from Rust", sent the reader hunting a backend bug. Met in practice while
    reproducing a build with two samples.
    """
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq

    from krewlyzer.pon.build import MIN_SAMPLES_PER_KEY, _compute_wps_baseline

    schema = pa.schema(
        [
            ("region_id", pa.string()),
            ("wps_nuc", pa.list_(pa.float32())),
            ("wps_tf", pa.list_(pa.float32())),
        ]
    )
    paths = []
    for i in range(MIN_SAMPLES_PER_KEY - 1):
        frame = pd.DataFrame(
            {
                "region_id": ["TSS|A|1"],
                "wps_nuc": [list(np.sin(np.arange(200) / 12.0 + i))],
                "wps_tf": [list(np.sin(np.arange(200) / 6.0))],
            }
        )
        path = tmp_path / f"s{i}.WPS.parquet"
        pq.write_table(pa.Table.from_pandas(frame, schema=schema), path)
        paths.append(str(path))

    with pytest.raises(RuntimeError) as excinfo:
        _compute_wps_baseline(paths)
    message = str(excinfo.value)
    assert str(MIN_SAMPLES_PER_KEY) in message, "does not name the floor"
    assert "Rust" not in message, "still blames the backend"


# ---------------------------------------------------------------------------
# the NRL search-band edge
# ---------------------------------------------------------------------------


def _bg(directory: Path, name: str, rows):
    """One sample's WPS_background: (group_id, nrl_bp, periodicity, at_limit)."""
    return _write(
        directory,
        name,
        group_id=[r[0] for r in rows],
        nrl_bp=[r[1] for r in rows],
        periodicity_score=[r[2] for r in rows],
        nrl_at_band_limit=[r[3] for r in rows],
    )


def test_the_band_limit_flag_is_required(tmp_path):
    """As `nrl_bp` and `periodicity_score` became in 4cd634b.

    Without it the builder cannot tell a repeat length from the edge of the
    window it searched, and silently averages the two.
    """
    path = _write(
        tmp_path, "s", group_id=["G"], nrl_bp=[250.0], periodicity_score=[0.5]
    )
    with pytest.raises(ValueError, match="nrl_at_band_limit"):
        _compute_wps_background_baseline([path])


def test_a_group_measured_by_nobody_keeps_its_row_with_no_nrl(tmp_path):
    """The decision: absent, not dropped.

    A group xs1 can measure and xs2 cannot is information when the two models
    are compared. Dropping the row destroys it and leaves the reader unable to
    distinguish "unmeasurable here" from "never attempted".
    """
    paths = [
        _bg(tmp_path, f"s{i}", [("Chr8_All", 250.0, 0.50 + i / 100, True)])
        for i in range(5)
    ]
    groups = _compute_wps_background_baseline(paths).groups
    row = groups[groups.group_id == "Chr8_All"]
    assert len(row) == 1, "the row was dropped instead of emptied"
    row = row.iloc[0]
    assert np.isnan(row.nrl_mean), f"reported the band edge as a mean: {row.nrl_mean}"
    assert np.isnan(row.nrl_std)
    assert row.n_at_band_limit == 5, "does not record why the NRL is absent"
    assert row.n_nrl_fitted == 0


def test_a_partly_limited_group_is_fitted_from_the_rows_that_measured_one(tmp_path):
    """Twelve of 28 xs2 duplex groups are partial, so this is the common case."""
    real = [181.0, 185.0, 190.0]
    rows = [("Global_All", 250.0, 0.5, True)] * 2 + [
        ("Global_All", v, 0.5 + i / 100, False) for i, v in enumerate(real)
    ]
    paths = [_bg(tmp_path, f"s{i}", [r]) for i, r in enumerate(rows)]
    row = _compute_wps_background_baseline(paths).groups.iloc[0]

    assert row.nrl_mean == pytest.approx(np.mean(real))
    assert row.nrl_std == pytest.approx(np.std(real, ddof=1))
    assert row.n_at_band_limit == 2 and row.n_nrl_fitted == 3
    assert row.nrl_mean < 250.0, "the band edge leaked into the mean"


def test_too_few_measured_rows_drops_the_mean_too_not_just_the_spread(tmp_path):
    """A cohort baseline averaged over one donor is a fabrication.

    `sample_std_or_nan` alone would leave `nrl_mean` present and `nrl_std` NaN,
    which reads as "we know the centre but not the spread" — from n=1. Five xs2
    duplex groups land here.
    """
    from krewlyzer.pon.build import MIN_SAMPLES_PER_KEY

    rows = [("Chr13_All", 250.0, 0.5, True)] * 6 + [("Chr13_All", 178.0, 0.5, False)]
    paths = [_bg(tmp_path, f"s{i}", [r]) for i, r in enumerate(rows)]
    row = _compute_wps_background_baseline(paths).groups.iloc[0]

    assert row.n_nrl_fitted == 1 < MIN_SAMPLES_PER_KEY
    assert np.isnan(row.nrl_mean), "one donor reported as a cohort baseline"
    assert np.isnan(row.nrl_std)


def test_periodicity_is_fitted_from_every_row_including_limited_ones(tmp_path):
    """Only the peak's *position* hit the edge. Its strength was measured.

    Across the xs2 duplex cohort the 174 band-limited rows carry 174 distinct
    periodicity values spanning 0.37-0.86 — real data. Excluding them with the
    NRL would discard 30% of the cohort for no reason, which is the mirror of
    the error this change fixes.
    """
    per = [0.41, 0.52, 0.47, 0.60]
    rows = [("G", 250.0, per[0], True), ("G", 250.0, per[1], True)] + [
        ("G", 180.0, per[2], False),
        ("G", 186.0, per[3], False),
    ]
    paths = [_bg(tmp_path, f"s{i}", [r]) for i, r in enumerate(rows)]
    row = _compute_wps_background_baseline(paths).groups.iloc[0]

    assert row.periodicity_mean == pytest.approx(np.mean(per))
    assert row.periodicity_std == pytest.approx(np.std(per, ddof=1))
    assert row.n_nrl_fitted == 2, "the NRL fit must still exclude them"


def test_the_periodicity_reader_no_longer_invents_a_standard_normal(tmp_path):
    """`row.get("periodicity_mean", 0), row.get("periodicity_std", 1)`.

    Those defaults made z equal the raw periodicity — ~0.47, an unremarkable
    number nobody would question. The same fabricated baseline as the block
    itself, one level down in the reader.
    """
    from krewlyzer.pon.model import WpsBackgroundBaseline

    baseline = WpsBackgroundBaseline(groups=pd.DataFrame({"group_id": ["G"]}))
    with pytest.raises(KeyError):
        baseline.get_periodicity_stats("G")
