"""The #77 donor-support measurement (Step 0b; reused unchanged in 6c).

Each test pins one rule from `.agents/plans/2026-11-19-step0b-77-pilot.md`
("Definitions"), on synthetic arrays small enough to reason about by hand.
"""

from __future__ import annotations

import numpy as np
import pytest

from krewlyzer.pon.support import (
    assign_bins,
    bin_edges,
    bootstrap_inflation,
    cell_spread,
    choose_threshold,
    donor_invariant,
    affordability,
    granularity,
    inflation_table,
    support_depth,
    support_noise,
    support_present,
)

N_PAT, N_BINS = 20, 20  # bins per anchor kept small; the module does not assume 200


def _cohort(support, sd_of_anchor, seed=0, mask=None):
    """z[patients, anchors, bins] ~ N(0, sd(anchor)), all anchors off-target unless masked."""
    rng = np.random.default_rng(seed)
    sd = np.asarray(sd_of_anchor, dtype=float)[None, :, None]
    z = rng.normal(0.0, 1.0, (N_PAT, len(support), N_BINS)) * sd
    off = np.ones(len(support), bool) if mask is None else mask
    return z, off


def _log_uniform(n, lo=0, hi=4, seed=1):
    return 10 ** np.random.default_rng(seed).uniform(lo, hi, n)


def _table(z, support, off, min_anchors=50):
    edges = bin_edges(support, min_anchors=min_anchors)
    spread, n_valid = cell_spread(z, k=10, seed=0)
    return edges, inflation_table(z, spread, n_valid, support, edges, off)


# --- support ---------------------------------------------------------------


def test_a_donor_with_no_row_counts_as_zero_support():
    ld = np.array([[np.nan, 5.0], [np.nan, 7.0], [3.0, np.nan]])
    assert support_depth(ld).tolist() == [0.0, 5.0]
    assert support_present(ld).tolist() == [1, 2]


def test_zero_support_is_its_own_bin_and_small_bins_merge():
    support = np.concatenate([np.zeros(300), _log_uniform(1000)])
    edges = bin_edges(support, min_anchors=200)
    bins = assign_bins(support, edges)
    assert set(bins[:300]) == {0}, "every zero-support anchor in bin 0"
    assert 0 not in set(bins[300:]), "no positive support shares bin 0"
    counts = np.bincount(bins)
    assert counts.min() >= 200, f"a bin under the minimum survived: {counts}"


# --- the rule --------------------------------------------------------------


def test_a_homogeneous_cohort_inflates_nowhere_and_withholds_nothing():
    support = _log_uniform(1200)
    z, off = _cohort(support, np.ones(len(support)))
    edges, table = _table(z, support, off)
    assert table.attrs["plateau_ok"]
    assert (table["inflation"].between(0.6, 1.6)).all(), table[["lo", "inflation"]]
    assert choose_threshold(table) == edges[0]


def test_a_noisy_low_bin_is_found_and_the_threshold_sits_above_it():
    support = _log_uniform(1200)
    sd = np.where(support < 10, 10.0, 1.0)
    z, off = _cohort(support, sd)
    edges, table = _table(z, support, off)
    low = table[table["hi"] <= 10]
    assert (low["inflation"] > 6).all(), low[["lo", "hi", "inflation"]]
    s_star = choose_threshold(table)
    assert s_star is not None and 9 < s_star <= 10**1.5


def test_a_low_bin_passing_by_chance_below_a_failing_one_does_not_set_the_threshold():
    # #77's own case: 0.30 next to 1.69. Bin 10-100 is noisy, bin 1-10 is clean.
    support = _log_uniform(1600)
    sd = np.where((support >= 10) & (support < 100), 10.0, 1.0)
    z, off = _cohort(support, sd)
    _, table = _table(z, support, off)
    s_star = choose_threshold(table)
    assert s_star is not None and s_star >= 100 * 0.99, s_star


def test_the_defect_itself_cannot_hide_in_dropped_cells():
    # Low-support anchors: most cells are ordinary, so the bin's spread looks
    # fine; but in a few, one patient carries a huge z and the rest have no row
    # -- #77's "a trace of coverage". The spread needs 10 patients, so those
    # cells are dropped and only the pooled tail can see them.
    support = _log_uniform(1200)
    z, off = _cohort(support, np.ones(len(support)))
    defect = np.flatnonzero(support < 10)[::5]  # a fifth of the low anchors
    z[:, defect, :] = np.nan
    z[0, defect, :] = 5_000.0
    edges, table = _table(z, support, off)
    low = table[table["hi"] <= 10]
    assert (
        low["inflation"] < 2
    ).all(), "the spread alone must not see it, or this test proves nothing"
    assert (low["dropped_cells"] > 0).all() and (low["dropped_abs_z_gt_100"] > 0).all()
    assert (low["pooled_tail"] > 100).all()
    s_star = choose_threshold(table)
    assert (
        s_star is not None and s_star >= 10 * 0.99
    ), "the bin hiding the defect must not pass"


def test_the_plateau_ignores_on_target_anchors():
    support = _log_uniform(1200)
    on = (support > 1000) & (np.arange(len(support)) % 2 == 0)
    sd = np.where(on, 5.0, 1.0)  # captured cancer genes: wide by biology
    z, _ = _cohort(support, sd)
    _, table = _table(z, support, ~on)
    assert table.attrs["plateau"] < 2.0, table.attrs["plateau"]


def test_a_plateau_whose_bins_disagree_is_not_a_plateau():
    support = _log_uniform(1200)
    sd = np.where(support > 10**3.5, 4.0, 1.0)
    z, off = _cohort(support, sd)
    _, table = _table(z, support, off)
    assert not table.attrs["plateau_ok"]
    assert choose_threshold(table) is None


def test_fewer_than_four_bins_is_refused():
    support = _log_uniform(300, lo=0, hi=1)
    z, off = _cohort(support, np.ones(len(support)))
    with pytest.raises(ValueError, match="four"):
        _table(z, support, off, min_anchors=100)


def test_the_bootstrap_resamples_patients():
    # One patient carries the whole tail of the low bin: resampling patients
    # makes that bin's interval wide; resampling cells would not.
    support = _log_uniform(1200)
    z, off = _cohort(support, np.ones(len(support)))
    low = support < 10
    z[0, low, :] *= 40.0
    edges = bin_edges(support, min_anchors=50)
    boot = bootstrap_inflation(z, support, edges, off, b=60, seed=3)
    low_rows = boot[boot["hi"] <= 10]
    high_rows = boot[boot["lo"] >= 1000]

    def width(rows):
        return (rows["inflation_hi"] / rows["inflation_lo"]).median()

    assert width(low_rows) > 2 * width(high_rows), boot


def test_two_different_cohorts_give_two_different_tables():
    support = _log_uniform(1200)
    za, off = _cohort(support, np.ones(len(support)), seed=1)
    zb, _ = _cohort(support, np.where(support < 100, 3.0, 1.0), seed=2)
    _, ta = _table(za, support, off)
    _, tb = _table(zb, support, off)
    assert not np.allclose(ta["p99"], tb["p99"])


# --- donor invariant -------------------------------------------------------


def test_the_donor_invariant_is_an_exact_identity_with_per_cell_counts():
    rng = np.random.default_rng(5)
    x = rng.normal(10, 2, (12, 40, N_BINS))
    x[rng.random(x.shape) < 0.3] = np.nan  # donors without a row: n varies per cell
    n = np.sum(~np.isnan(x), axis=0)
    mu = np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0, ddof=1)
    sd[n < 2] = np.nan
    z = (x - mu) / sd
    observed, required = donor_invariant(z)
    assert abs(observed - required) < 1e-12
    assert abs(required - np.sqrt(np.sum(n[n >= 2] - 1) / np.sum(n[n >= 2]))) < 1e-12


# --- granularity, affordability, the noise key -----------------------------


def test_granularity_is_the_share_of_extreme_z_below_the_threshold():
    support = np.array([1.0, 2.0, 50.0, 80.0])
    z = np.zeros((3, 4, 2))
    z[0, 0, 0] = z[1, 1, 1] = z[2, 1, 0] = 500.0  # three extremes below S* = 10
    z[0, 3, 1] = -300.0  # one above
    g = granularity(z, support, s_star=10.0)
    assert g["n_extreme"] == 4
    assert g["share_extreme_below"] == pytest.approx(0.75)
    assert g["share_cells_below"] == pytest.approx(0.5)


def test_affordability_compares_cells_with_cells_over_scorable_cells():
    support = np.array([1.0, 100.0])
    pon_mean = np.array([[0.05, 5.0], [0.01, 3.0]])
    pon_sigma = np.array([[1.0, 1.0], [1.0, np.nan]])  # one cell not scorable
    a = affordability(support, pon_mean, pon_sigma, s_star=10.0)
    assert a["scorable_cells"] == 3
    assert a["share_withheld_support"] == pytest.approx(2 / 3)  # anchor 0, both cells
    assert a["share_withheld_mu_cut"] == pytest.approx(
        2 / 3
    )  # |mu| < 0.1 at two scorable cells


def test_the_noise_key_is_donor_sigma_over_single_sample_noise():
    rng = np.random.default_rng(9)
    pon_mean = rng.normal(0, 5, (30, 50))
    wps = pon_mean[None] + rng.normal(0, 2.0, (8, 30, 50))
    donor_sigma = np.full((30, 50), 4.0)
    ratio = support_noise(donor_sigma, wps, pon_mean)
    assert ratio.shape == (30,)
    assert np.allclose(ratio, 2.0, rtol=0.15)
