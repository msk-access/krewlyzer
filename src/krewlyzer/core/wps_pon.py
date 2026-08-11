"""Apply the PON to WPS output — the largest baseline, previously unread.

`wps_baseline` is ~128k anchors of 200-element mean and sigma vectors, roughly
90% of every PON file. Until this module its only consumer was a log line
appending ``"WPS"`` to a list of available components: the baseline was built,
stored and shipped, and nothing downstream carried a single value derived from
it.

Two kinds of output, and the distinction is the whole design.

**Per-position z vectors** (``wps_nuc_z``, ``wps_tf_z``). Elementwise
``(x - mean) / std`` against the baseline profile. This is the raw, correct
object: the *shape* is the measurement, and a shape survives no summary.

**Derived shape quantities**, each z-scored against a baseline of *itself*.
Never a reduction of the z vector — adjacent WPS positions have lag-1
autocorrelation **0.986**, because a fragment spans ~167 bp and contributes to
many consecutive positions. So:

- a mean of z over 200 positions has nothing like ``sigma/sqrt(200)``
  precision, and reporting it as a z would be badly overconfident;
- a max of \\|z\\| over 200 positions has an expected value of 2.97 under pure
  noise, so a ``|z| > 2`` rule would flag nearly every anchor.

Both were the obvious first design and both are wrong. Derive the biological
quantity first, then z-score that.

**Displacement is measured but not z-scored.** ``wps_phase_shift_bp`` is
emitted because it is cheap and genuinely non-redundant (correlation -0.24 and
-0.28 with the two scored statistics), but it gets no baseline: measured on a
real cohort, per-sample mean lag varies by 0.26 bp against a within-sample
spread of 8.43, so there is no whole-sample phasing signal, and per anchor the
intraclass correlation is 0.479 -- about half of any lag is noise, optimistically,
since that estimate used a baseline containing the samples being scored. A
z-score of an integer-valued statistic that is half noise is a plausible number
and nothing more.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from krewlyzer.pon.model import WPS_SHAPE_STATS

from .output_utils import read_exact_table, write_table

logger = logging.getLogger("core.wps_pon")

#: How far the phase search looks, in positions.
#:
#: Python-only: the builder no longer computes a phase baseline, so there is no
#: Rust twin this has to agree with. Measured on the real cohort, the search
#: terminates on its own edge for 3.3% of anchors at +/-20 and 1.8% at +/-30.
PHASE_MAX_LAG = 30


def log_amplitude(profile: np.ndarray) -> float:
    """Peak-to-trough range, log scale — see the Rust twin for why log."""
    finite = profile[np.isfinite(profile)]
    if finite.size < 2:
        return float("nan")
    return float(np.log1p(finite.max() - finite.min()))


def shape_correlation(profile: np.ndarray, baseline: np.ndarray) -> float:
    """Pearson r between the sample profile and the baseline's mean profile."""
    ok = np.isfinite(profile) & np.isfinite(baseline)
    if ok.sum() < 3:
        return float("nan")
    a, b = profile[ok], baseline[ok]
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def fisher_z(r: float) -> float:
    """``arctanh(r)``. Required, not cosmetic.

    A correlation is bounded at 1.0. Measured on the real cohort the shape
    correlation sits at mean 0.844 with sigma 0.099, so the largest attainable
    positive z is about 1.5 and **302 of 400 anchors could not reach +2**
    however tumour-like the sample. On the Fisher scale the ceiling is gone.
    """
    if not np.isfinite(r):
        return float("nan")
    return float(np.arctanh(np.clip(r, -0.999999, 0.999999)))


def phase_shift(
    profile: np.ndarray, baseline: np.ndarray, max_lag: int = PHASE_MAX_LAG
) -> tuple[float, bool]:
    """Displacement of the sample profile against the baseline, in positions.

    Returns ``(lag, hit_limit)``. A shift is invisible to any per-position
    summary: the same profile one nucleosome along differs at every position
    while being unchanged in shape.

    ``hit_limit`` is ``nrl_at_band_limit`` one level down — a search that ends
    on its own window edge has found a boundary, not a measurement. Measured
    incidence at this window: 1.8% of anchors.
    """
    n = min(profile.size, baseline.size)
    if n < 4 * max_lag:
        return float("nan"), False
    core = baseline[max_lag : n - max_lag]
    best_r, best_lag = -np.inf, 0
    for k in range(-max_lag, max_lag + 1):
        window = profile[max_lag + k : max_lag + k + core.size]
        r = shape_correlation(window, core)
        if np.isfinite(r) and r > best_r:
            best_r, best_lag = r, k
    if not np.isfinite(best_r):
        return float("nan"), False
    return float(best_lag), abs(best_lag) == max_lag


def apply_wps_pon(
    wps_path: Path,
    pon,
    output_base: Optional[Path] = None,
    column: str = "wps_nuc",
    baseline_attr: str = "wps_baseline",
) -> int:
    """Add PON-derived columns to a WPS table. Returns anchors scored.

    ``baseline_attr`` selects which vector baseline to score against:
    ``wps_baseline`` for the genome-wide anchors, ``wps_baseline_panel`` for
    the assay-specific ones. The panel baseline was built and stored by every
    panel-mode PON and read by nothing, so the ``.WPS.panel`` output shipped
    raw -- the anchors closest to the targeted regions were the only ones with
    no comparison to a healthy cohort.

    Emits, per anchor:

    ``{column}_z``                 200-element z vector against the profile
    ``wps_log_amplitude``          + ``_z``
    ``wps_shape_corr``             + ``_z`` (z computed on the Fisher scale)
    ``wps_phase_shift_bp``         raw displacement, deliberately not z-scored
    ``wps_phase_at_search_limit``  the boundary flag

    The raw derived values are emitted beside their z-scores deliberately: a z
    is uninterpretable without a PON, and these three are readable on their own.
    """
    frame = read_exact_table(wps_path)
    if frame is None:
        logger.warning(f"WPS output not found for PON scoring: {wps_path}")
        return 0
    if "region_id" not in frame.columns or column not in frame.columns:
        logger.warning(f"WPS output lacks region_id/{column}: {wps_path}")
        return 0

    vector_baseline = getattr(pon, baseline_attr, None)
    # The shape statistics are fitted on the genome-wide anchor set only. The
    # panel anchors overlap it and number a few hundred, too few to justify a
    # second block, so they borrow it rather than go unscored.
    shape_baseline = getattr(pon, "wps_shape_baseline", None)
    if vector_baseline is None:
        logger.warning(
            f"PON has no {baseline_attr}; WPS keeps its raw profile and gets "
            "no z-scores. Rebuild the PON with build-pon."
        )
        return 0
    if shape_baseline is None:
        logger.info(
            "PON has no wps_shape_baseline (built before 0.9.0); per-position "
            "z vectors are still written, the derived shape z-scores are not."
        )

    z_vectors: list = []
    amps: list[float] = []
    corrs: list[float] = []
    shifts: list[float] = []
    at_limit: list[bool] = []
    z_amp: list[float] = []
    z_corr: list[float] = []
    n_absent = 0
    n_scored = 0
    n_at_limit = 0

    for region_id, raw in zip(frame["region_id"], frame[column]):
        profile = np.asarray(raw, dtype=float)
        mean = vector_baseline.get_baseline_vector(str(region_id), column)

        amplitude = log_amplitude(profile)
        amps.append(amplitude)

        if mean is None:
            n_absent += 1
            z_vectors.append(None)
            corrs.append(float("nan"))
            shifts.append(float("nan"))
            at_limit.append(False)
            z_amp.append(float("nan"))
            z_corr.append(float("nan"))
            continue

        mean = np.asarray(mean, dtype=float)
        z_vectors.append(
            vector_baseline.compute_z_vector(str(region_id), profile, column)
        )
        r = shape_correlation(profile, mean)
        corrs.append(r)
        lag, hit = phase_shift(profile, mean)
        shifts.append(lag)
        at_limit.append(hit)
        n_at_limit += int(hit)
        n_scored += 1

        if shape_baseline is None:
            z_amp.append(float("nan"))
            z_corr.append(float("nan"))
            continue
        observed = {
            "log_amplitude": amplitude,
            # z on the Fisher scale, reported next to the readable correlation
            "shape_corr_fisher": fisher_z(r),
        }
        scores = {
            stat: shape_baseline.compute_zscore(str(region_id), stat, observed[stat])
            for stat in WPS_SHAPE_STATS
        }
        z_amp.append(_nan_if_none(scores["log_amplitude"]))
        z_corr.append(_nan_if_none(scores["shape_corr_fisher"]))

    frame[f"{column}_z"] = z_vectors
    frame["wps_log_amplitude"] = amps
    frame["wps_log_amplitude_z"] = z_amp
    frame["wps_shape_corr"] = corrs
    frame["wps_shape_corr_z"] = z_corr
    frame["wps_phase_shift_bp"] = shifts
    frame["wps_phase_at_search_limit"] = at_limit

    # Parquet, always -- not whatever `--output-format` says.
    #
    # WPS is the one family the output contract pins to Parquet
    # (docs/reference/output-files.md), because the 200-element vectors are
    # unusable as text. The Rust step writes `.WPS.parquet` unconditionally,
    # so honouring `--output-format` here wrote the z-scores to a *second*
    # file, `.WPS.tsv`, and left the parquet as the raw profile.
    #
    # Downstream reads Parquet only (invariant #2). Measured on a real run:
    # `.WPS.tsv` carried all seven PON columns, `.WPS.parquet` carried none,
    # and the consumer would have read the one without them -- WPS scoring
    # present on disk and absent from the product.
    write_table(
        frame,
        output_base or wps_path.with_suffix(""),
        output_format="parquet",
        compress=False,
    )

    logger.info(f"WPS PON: {n_scored}/{len(frame)} anchors scored ({wps_path.name})")
    if n_absent:
        logger.warning(
            f"WPS PON: {n_absent}/{len(frame)} anchors are absent from the "
            "baseline and get no z-score. Anchors backed by fewer than 3 "
            "samples are dropped at build time, so this is expected to be "
            "large for duplex PONs."
        )
    if n_at_limit:
        logger.warning(
            f"WPS PON: {n_at_limit}/{len(frame)} anchors hit the +/-"
            f"{PHASE_MAX_LAG} phase-search window. wps_phase_shift_bp is the "
            "edge of the search there, not a measurement -- see "
            "wps_phase_at_search_limit."
        )
    return n_scored


def _nan_if_none(value: Optional[float]) -> float:
    """An anchor absent from the shape baseline reads the same as an unmeasured one."""
    return float("nan") if value is None else float(value)
