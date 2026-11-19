"""Donor support for WPS z: where the PON baseline can serve as a yardstick (#77).

#77: a WPS z-score is computed against a donor baseline at every anchor, but at
anchors no donor covered the baseline is not a measurement, and a patient with a
trace of coverage there gets |z| in the thousands. This module measures where
that happens, so a support threshold ``S*`` can be set from data and the
withheld z can be labelled rather than shipped.

Pure functions over dense arrays -- no paths, no I/O -- so the Step 0b pilot,
the 6c decision on the final code, and Step 6's real-data test compute the same
quantities. Definitions follow ``.agents/plans/2026-11-19-step0b-77-pilot.md``:

* arrays are dense over an assay's anchor universe; NaN means "no row" -- the
  WPS writer omits an anchor with no coverage rather than writing zeros;
* ``z`` is ``[patients, anchors, bins]``; a NaN z is "no row", or a cell the
  PON cannot score.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

#: |z| above which a value cannot be a deviation from a healthy cohort; #77's
#: own marker of the defect.
EXTREME_Z = 100.0


# --------------------------------------------------------------------------
# Support keys
# --------------------------------------------------------------------------


def support_depth(local_depth: np.ndarray) -> np.ndarray:
    """Median donor ``local_depth`` per anchor, over the **full** donor list.

    ``local_depth`` is ``[donors, anchors]`` with NaN where a donor has no row.
    No row means no coverage, so it counts as 0: a median over only the donors
    present could never be 0, and the anchors #77 is about -- covered by almost
    no donor -- would look supported by the one that did.
    """
    return np.median(np.nan_to_num(local_depth, nan=0.0), axis=0)


def support_present(local_depth: np.ndarray) -> np.ndarray:
    """Number of donors with a row at each anchor."""
    return np.sum(~np.isnan(local_depth), axis=0)


def support_noise(
    donor_sigma: np.ndarray, patient_wps: np.ndarray, pon_mean: np.ndarray
) -> np.ndarray:
    """Donor σ against single-sample noise, per anchor (#77 requirement 2).

    Single-sample noise is a patient's SD across an anchor's bins after the PON
    mean profile is removed; the median over patients is the anchor's noise.
    Returns the median over bins of ``donor_sigma / noise``.
    """
    resid = patient_wps - pon_mean[None]
    with np.errstate(invalid="ignore", divide="ignore"):
        noise = np.nanmedian(np.nanstd(resid, axis=2, ddof=1), axis=0)
        return np.nanmedian(donor_sigma / noise[:, None], axis=1)


# --------------------------------------------------------------------------
# Bins over support
# --------------------------------------------------------------------------


def bin_edges(
    support: np.ndarray, min_anchors: int = 200, step: float = 0.5
) -> List[float]:
    """Lower edges of half-decade bins of ``support``, plus ``inf``.

    Zero support is its own bin, ``[0, first positive edge)``, never
    ``log10(0)``. A bin with fewer than ``min_anchors`` anchors is merged into
    its neighbour towards the median -- counted in anchors, because the bins of
    one anchor are nearly identical and a cell count would overstate the data.
    """
    support = np.asarray(support, dtype=float)
    positive = support[support > 0]
    bins: List[List[float]] = []  # [lower edge, anchor count]
    n_zero = int(np.sum(support <= 0))
    if n_zero:
        bins.append([0.0, n_zero])
    if positive.size:
        keys = np.floor(np.log10(positive) / step).astype(int)
        for k in range(int(keys.min()), int(keys.max()) + 1):
            bins.append([10 ** (k * step), int(np.sum(keys == k))])
    median = float(np.median(support)) if support.size else 0.0
    while len(bins) > 1:
        small = [i for i, (_, n) in enumerate(bins) if n < min_anchors]
        if not small:
            break
        i = small[0]
        his = [b[0] for b in bins[1:]] + [np.inf]
        m = next(
            j
            for j, (lo, hi) in enumerate(zip([b[0] for b in bins], his))
            if lo <= median < hi
        )
        if i < m or (i == m and i + 1 < len(bins)):
            bins[i + 1] = [
                bins[i][0],
                bins[i][1] + bins[i + 1][1],
            ]  # upper neighbour grows down
        else:
            bins[i - 1][1] += bins[i][1]  # lower neighbour grows up
        del bins[i]
    return [b[0] for b in bins] + [float("inf")]


def assign_bins(support: np.ndarray, edges: List[float]) -> np.ndarray:
    """Bin index of each anchor."""
    return (
        np.searchsorted(
            np.asarray(edges), np.asarray(support, dtype=float), side="right"
        )
        - 1
    )


# --------------------------------------------------------------------------
# Spread, inflation, threshold
# --------------------------------------------------------------------------


def cell_spread(
    z: np.ndarray, k: int = 10, seed: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """SD of z over exactly ``k`` patients per cell, drawn at random from those
    with a value; NaN where fewer than ``k`` have one.

    A fixed ``k`` because a p99 of SDs grows with the number of patients, and
    cells differ in how many patients have a row. Returns ``(spread, n_valid)``.
    """
    valid = ~np.isnan(z)
    n_valid = valid.sum(axis=0)
    if z.shape[0] < k:
        return np.full(z.shape[1:], np.nan), n_valid
    priority = np.random.default_rng(seed).random(z.shape, dtype=np.float32)
    priority[~valid] = np.inf
    pick = np.argpartition(priority, k - 1, axis=0)[:k]
    sample = np.take_along_axis(z, pick, axis=0)
    with np.errstate(invalid="ignore"):
        spread = np.std(sample, axis=0, ddof=1)
    spread[n_valid < k] = np.nan
    return spread, n_valid


def inflation_table(
    z: np.ndarray,
    spread: np.ndarray,
    n_valid: np.ndarray,
    support: np.ndarray,
    edges: List[float],
    off_target: np.ndarray,
    k: int = 10,
    plateau_bins: int = 3,
    plateau_tolerance: float = 1.5,
    min_bins: int = 4,
) -> pd.DataFrame:
    """Per support bin: the p99 cell spread, its inflation over the plateau,
    the pooled tail, and what the k-patient rule dropped.

    The plateau is the p99 spread of the top ``plateau_bins`` bins, over
    **off-target** anchors only: the best-supported anchors are the captured
    genes, where patient biology widens z. ``attrs["plateau_ok"]`` is False
    when those bins disagree by more than ``plateau_tolerance`` -- then there
    is nothing stable to normalise against.

    The pooled tail (p99.9 of |z| over every non-null value) drops nothing, so
    a bin whose worst cells lack ``k`` patients -- #77's single patient with a
    trace of coverage -- still shows its tail.
    """
    n_bins = len(edges) - 1
    if n_bins < min_bins:
        raise ValueError(
            f"{n_bins} support bins after merging; at least four are needed"
        )
    bins = assign_bins(support, edges)
    absz = np.abs(z)
    rows = []
    for i in range(n_bins):
        a = bins == i
        cells = spread[a]
        finite = cells[np.isfinite(cells)]
        vals = absz[:, a, :]
        vals = vals[~np.isnan(vals)]
        dropped = (n_valid[a] > 0) & (n_valid[a] < k)
        rows.append(
            {
                "bin": i,
                "lo": edges[i],
                "hi": edges[i + 1],
                "anchors": int(a.sum()),
                "share_anchors": float(a.mean()),
                "cells": int(finite.size),
                "p99": float(np.percentile(finite, 99)) if finite.size else np.nan,
                "pooled_tail": (
                    float(np.percentile(vals, 99.9)) if vals.size else np.nan
                ),
                "share_abs_z_gt_100": (
                    float(np.mean(vals > EXTREME_Z)) if vals.size else np.nan
                ),
                "dropped_cells": int(dropped.sum()),
                "dropped_abs_z_gt_100": int(
                    np.sum(absz[:, a, :][:, dropped] > EXTREME_Z)
                ),
            }
        )
    table = pd.DataFrame(rows)

    top = bins >= n_bins - plateau_bins
    plateau_cells = top & off_target
    pooled = spread[plateau_cells]
    pooled = pooled[np.isfinite(pooled)]
    plateau = float(np.percentile(pooled, 99)) if pooled.size else np.nan
    tail_vals = absz[:, plateau_cells, :]
    tail_vals = tail_vals[~np.isnan(tail_vals)]
    plateau_tail = float(np.percentile(tail_vals, 99.9)) if tail_vals.size else np.nan
    per_bin = []
    for i in range(n_bins - plateau_bins, n_bins):
        c = spread[(bins == i) & off_target]
        c = c[np.isfinite(c)]
        per_bin.append(float(np.percentile(c, 99)) if c.size else np.nan)
    ratios = np.array(per_bin) / plateau
    plateau_ok = bool(
        np.isfinite(plateau)
        and np.all(np.isfinite(ratios))
        and np.all((ratios <= plateau_tolerance) & (ratios >= 1 / plateau_tolerance))
    )
    table["inflation"] = table["p99"] / plateau
    table["tail_ratio"] = table["pooled_tail"] / plateau_tail
    table.attrs.update(
        plateau=plateau, plateau_tail=plateau_tail, plateau_ok=plateau_ok
    )
    return table


def choose_threshold(
    table: pd.DataFrame, max_inflation: float = 2.0, max_tail_ratio: float = 2.0
) -> Optional[float]:
    """``S*``: the lower edge of the lowest bin such that it **and every higher
    bin** has inflation and tail ratio within bounds; None when no bin
    qualifies or there is no plateau.

    "Every higher bin" because a low bin can pass by chance beneath a failing
    one (#77 saw 0.30 next to 1.69). A NaN inflation -- no cell had enough
    patients -- fails: the absence of a measurement is not a pass.
    """
    if not table.attrs.get("plateau_ok", False):
        return None
    t = table.sort_values("lo")
    passes = (t["inflation"] <= max_inflation) & (t["tail_ratio"] <= max_tail_ratio)
    ok = passes.fillna(False).to_numpy(dtype=bool)
    for i in range(len(ok)):
        if ok[i:].all():
            return float(t["lo"].iloc[i])
    return None


def bootstrap_inflation(
    z: np.ndarray,
    support: np.ndarray,
    edges: List[float],
    off_target: np.ndarray,
    b: int = 200,
    seed: int = 0,
    k: int = 10,
) -> pd.DataFrame:
    """95% percentile intervals on each bin's inflation, resampling **patients**.

    Patients, not cells: the cells of one patient are correlated, so
    resampling cells would report an interval far narrower than the data
    supports, and miss a tail carried by one patient.
    """
    rng = np.random.default_rng(seed)
    draws = []
    for rep in range(b):
        idx = rng.integers(0, z.shape[0], z.shape[0])
        zb = z[idx]
        spread, n_valid = cell_spread(zb, k=k, seed=rep)
        draws.append(
            inflation_table(zb, spread, n_valid, support, edges, off_target, k=k)[
                "inflation"
            ].to_numpy()
        )
    arr = np.vstack(draws)
    out = pd.DataFrame({"lo": edges[:-1], "hi": edges[1:]})
    out["inflation_lo"] = np.nanpercentile(arr, 2.5, axis=0)
    out["inflation_hi"] = np.nanpercentile(arr, 97.5, axis=0)
    return out


# --------------------------------------------------------------------------
# Donor invariant, granularity, affordability
# --------------------------------------------------------------------------


def donor_invariant(z_donors: np.ndarray) -> Tuple[float, float]:
    """Observed and required RMS(z) of a PON's own donors scored against it.

    With ``z = (x − mean) / sd`` and ``sd`` the ddof-1 SD over a cell's ``n_c``
    donors, ``Σ z² = n_c − 1`` in every cell: an identity, so the required
    value uses each cell's own donor count -- not the cohort size, since the
    builder sees only the donors with a row there.
    """
    finite = np.isfinite(z_donors)
    n = finite.sum(axis=0)
    scored = n >= 2
    observed = float(np.sqrt(np.mean(z_donors[finite & scored[None]] ** 2)))
    required = float(np.sqrt(np.sum(n[scored] - 1) / np.sum(n[scored])))
    return observed, required


def granularity(z: np.ndarray, support: np.ndarray, s_star: float) -> Dict[str, float]:
    """Of patient × cell values with |z| > 100, the share in anchors below ``S*``.

    Read beside ``share_cells_below``: when nearly every cell is below ``S*``
    the first share is high by construction and says nothing.
    """
    below = np.asarray(support) < s_star
    extreme = np.abs(np.nan_to_num(z, nan=0.0)) > EXTREME_Z
    n_extreme = int(extreme.sum())
    scorable = np.any(~np.isnan(z), axis=0)
    return {
        "n_extreme": n_extreme,
        "share_extreme_below": (
            float(extreme[:, below, :].sum() / n_extreme) if n_extreme else float("nan")
        ),
        "share_cells_below": (
            float(scorable[below].sum() / scorable.sum())
            if scorable.any()
            else float("nan")
        ),
    }


def affordability(
    support: np.ndarray,
    pon_mean: np.ndarray,
    pon_sigma: np.ndarray,
    s_star: float,
    z: Optional[np.ndarray] = None,
    mu_cut: float = 0.1,
) -> Dict[str, float]:
    """Cells withheld at ``S*`` beside cells the ``|μ| < mu_cut`` cut withholds,
    both over the PON's scorable cells (finite σ), so the two are in one unit."""
    scorable = np.isfinite(pon_sigma)
    by_support = scorable & (np.asarray(support) < s_star)[:, None]
    by_mu = scorable & (np.abs(pon_mean) < mu_cut)
    n = int(scorable.sum())
    out: Dict[str, float] = {
        "scorable_cells": n,
        "share_withheld_support": float(by_support.sum() / n) if n else float("nan"),
        "share_withheld_mu_cut": float(by_mu.sum() / n) if n else float("nan"),
    }
    if z is not None:
        has = ~np.isnan(z)
        total = int(has.sum())
        out["share_values_withheld_support"] = (
            float((has & by_support[None]).sum() / total) if total else float("nan")
        )
        out["share_values_withheld_mu_cut"] = (
            float((has & by_mu[None]).sum() / total) if total else float("nan")
        )
    return out
