"""WPS z-scoring — the largest PON baseline, previously read by nothing.

`wps_baseline` is ~128k anchors of 200-element mean and sigma vectors, roughly
90% of every PON file, and until this its only consumer was a log line
appending `"WPS"` to a list of available components.

The design question was what to emit, and measuring answered it against the
obvious choices:

- **Not a mean of z over positions.** Adjacent WPS positions have lag-1
  autocorrelation 0.986 — a fragment spans ~167 bp and touches many at once —
  so such a mean has nothing like `sigma/sqrt(200)` precision.
- **Not a max of |z| over positions.** Expected max under pure noise is 2.97,
  so `|z| > 2` would flag nearly every anchor.
- **Not centre-versus-flank.** TSS anchors dip at the centre (−6.8 against
  −3.4 in the flanks); CTCF anchors do the opposite. Any fixed window is
  backwards for one of the two.

So: per-position z vectors, plus three window-free derived quantities each
z-scored against a baseline of itself.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from krewlyzer.core.wps_pon import (
    PHASE_MAX_LAG,
    fisher_z,
    log_amplitude,
    phase_shift,
    shape_correlation,
)

pytestmark = pytest.mark.unit


def _wave(n=200, period=12.0, offset=0.0, scale=1.0):
    return scale * np.sin((np.arange(n) - offset) / period)


# ---------------------------------------------------------------------------
# log amplitude
# ---------------------------------------------------------------------------


def test_log_amplitude_is_not_a_coverage_measurement():
    """Raw peak-to-trough range measures how deeply the sample was sequenced.

    Measured on the real cohort: raw amplitude correlates +0.512 with
    `local_depth` and is skewed 11.6; the log form drops that to −0.036 and
    1.6. A z on the raw scale would rank samples by depth.
    """
    shallow, deep = _wave(), _wave(scale=100.0)
    a, b = log_amplitude(shallow), log_amplitude(deep)
    assert b > a, "a deeper profile must still register as larger"
    assert b / a < 10, f"100x the depth must not be ~100x the statistic ({b / a:.1f}x)"


def test_log_amplitude_is_nan_without_a_profile():
    assert math.isnan(log_amplitude(np.array([])))
    assert math.isnan(log_amplitude(np.array([1.0])))


# ---------------------------------------------------------------------------
# shape correlation and the Fisher transform
# ---------------------------------------------------------------------------


def test_the_fisher_transform_removes_a_binding_ceiling():
    """Required, not cosmetic.

    A correlation is bounded at 1.0. On the real cohort the shape correlation
    sits at mean 0.844 with sigma 0.099, so the largest attainable positive z
    is ~1.5 and **302 of 400 anchors could not reach +2** however tumour-like
    the sample. This is the ceiling I wrongly dismissed for MDS and which is
    genuinely binding here.
    """
    mean_r, sigma_r = 0.844, 0.0985
    assert (1.0 - mean_r) / sigma_r < 2.0, "the raw ceiling should be binding"
    # On the Fisher scale there is no bound to run into.
    assert fisher_z(0.999) > fisher_z(0.99) > fisher_z(0.9)
    assert math.isfinite(fisher_z(1.0)), "a perfect correlation must stay a number"
    assert math.isnan(fisher_z(float("nan")))


def test_shape_correlation_is_nan_without_variance():
    flat = np.full(200, 3.0)
    assert math.isnan(shape_correlation(flat, _wave()))
    assert shape_correlation(_wave(), _wave()) == pytest.approx(1.0)


def test_shape_correlation_sees_a_wrong_shape():
    """The point of the statistic: same amplitude, different structure."""
    assert shape_correlation(_wave(period=12), _wave(period=40)) < 0.9


# ---------------------------------------------------------------------------
# phase shift
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("offset", [-9, -3, 0, 5, 11])
def test_phase_shift_recovers_a_known_displacement(offset):
    """A shift is invisible to every per-position summary.

    The same profile one nucleosome along differs at every position while
    being unchanged in shape, so no reduction of the z vector can see it.
    """
    lag, hit = phase_shift(_wave(offset=offset), _wave())
    assert not hit
    assert abs(lag - offset) <= 1, f"offset {offset} recovered as {lag}"


def test_phase_shift_reports_hitting_its_own_search_window():
    """`nrl_at_band_limit` one level down.

    A search that stops on its own edge has found a boundary, not a
    measurement. Measured incidence on the real cohort: 1.8% of anchors.
    """
    ramp = np.arange(200, dtype=float)
    _, hit = phase_shift(-ramp, ramp)
    assert hit, "a search ending on the edge must say so"


def test_a_boundary_shift_is_not_z_scored():
    """The flag is not decorative: the value behind it is excluded from the z.

    Scoring the edge of a search window against a baseline of real shifts
    would turn "we stopped looking" into "displaced by exactly 30 bp".
    """
    from krewlyzer.core import wps_pon

    assert PHASE_MAX_LAG == 30, "must match PHASE_MAX_LAG in rust/src/pon_builder.rs"
    source = wps_pon.apply_wps_pon.__doc__ or ""
    assert "wps_phase_at_search_limit" in source


# ---------------------------------------------------------------------------
# end to end
# ---------------------------------------------------------------------------


def _pon_from(frames):
    """A PON built from in-memory WPS frames, via the real builder."""
    import tempfile
    from pathlib import Path

    from krewlyzer.pon.build import _compute_wps_baseline
    from krewlyzer.pon.model import PonModel

    import pyarrow as pa
    import pyarrow.parquet as pq

    # float32 lists, matching what the Rust writer produces. The reader asks
    # for List<Float32> specifically, and pandas would hand it List<Double>.
    schema = pa.schema(
        [
            ("region_id", pa.string()),
            ("wps_nuc", pa.list_(pa.float32())),
            ("wps_tf", pa.list_(pa.float32())),
        ]
    )
    directory = Path(tempfile.mkdtemp())
    paths = []
    for i, frame in enumerate(frames):
        path = directory / f"s{i}.WPS.parquet"
        pq.write_table(pa.Table.from_pandas(frame, schema=schema), path)
        paths.append(str(path))
    vector, shape = _compute_wps_baseline(paths)
    return (
        PonModel(
            schema_version="1.0",
            assay="xs2",
            build_date="2026-01-01",
            n_samples=len(frames),
            reference="r",
            panel_mode=True,
            target_regions_file="t",
            wps_baseline=vector,
            wps_shape_baseline=shape,
        ),
        directory,
        paths,
    )


def _frame(seed: int, n_anchors: int = 8):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "region_id": [f"TSS|G{i}|T{i}" for i in range(n_anchors)],
            "wps_nuc": [
                _wave(offset=rng.normal(0, 1.5)) + rng.normal(0, 0.05, 200)
                for _ in range(n_anchors)
            ],
            "wps_tf": [_wave(period=6) for _ in range(n_anchors)],
        }
    )


def test_every_promised_column_is_emitted():
    from pathlib import Path

    from krewlyzer.core.wps_pon import apply_wps_pon

    pon, directory, paths = _pon_from([_frame(s) for s in range(5)])
    apply_wps_pon(
        Path(paths[0]), pon, output_base=directory / "o", output_format="parquet"
    )
    out = pd.read_parquet(directory / "o.parquet")
    for column in (
        "wps_nuc_z",
        "wps_log_amplitude",
        "wps_log_amplitude_z",
        "wps_shape_corr",
        "wps_shape_corr_z",
        "wps_phase_shift_bp",
        "wps_phase_at_search_limit",
    ):
        assert column in out.columns, f"{column} was not emitted"
    assert len(np.asarray(out["wps_nuc_z"].iloc[0])) == 200
    assert "wps_phase_shift_z" not in out.columns, (
        "displacement is measured but deliberately not z-scored -- see "
        "test_displacement_is_measured_but_not_scored"
    )


def test_displacement_is_measured_but_not_scored():
    """The raw lag is useful; a z-score of it would not be.

    Measured on a real cohort: per-sample mean lag varies by 0.26 bp against a
    within-sample spread of 8.43, so there is no whole-sample phasing signal;
    and per anchor the intraclass correlation is 0.479, meaning about half of
    any lag is noise. That estimate is optimistic -- the baseline used
    contained the samples being scored.

    It is also integer-valued, so on a small cohort its sigma bottoms out at
    std([0,0,0,0,0,1]) = 0.408 and a 1 bp shift scores z = 2.4.

    So the column ships and the baseline does not. The measurement is cheap
    and genuinely non-redundant (corr -0.24 and -0.28 with the two scored
    statistics); adding the baseline back is small if the n=21/47 rebuild
    shows reproducible per-anchor shifts.
    """
    from krewlyzer.pon.model import WPS_SHAPE_STATS

    assert (
        "phase_shift_bp" not in WPS_SHAPE_STATS
    ), "the shape baseline must not carry a phase-shift entry"
    assert set(WPS_SHAPE_STATS) == {"log_amplitude", "shape_corr_fisher"}


def test_a_pon_without_a_shape_block_still_writes_the_z_vectors():
    """A PON built before 0.9.0 has `wps_baseline` and no `wps_shape_baseline`.

    The per-position z needs only the vector baseline, so it must not be lost
    along with the derived scores.
    """
    from pathlib import Path

    from krewlyzer.core.wps_pon import apply_wps_pon

    pon, directory, paths = _pon_from([_frame(s) for s in range(5)])
    pon.wps_shape_baseline = None
    apply_wps_pon(
        Path(paths[0]), pon, output_base=directory / "o", output_format="parquet"
    )
    out = pd.read_parquet(directory / "o.parquet")
    assert out["wps_nuc_z"].notna().any(), "the z vectors must survive"
    assert (
        pd.to_numeric(out["wps_log_amplitude_z"], errors="coerce").isna().all()
    ), "the derived z-scores have no baseline and must be NaN"
    assert out["wps_log_amplitude"].notna().any(), "the raw value is still readable"


def test_scoring_a_sample_inside_its_own_baseline_shrinks_z():
    """Why self-inclusion cannot be used to claim calibration.

    A sample contributing to its own baseline pulls the mean toward itself, so
    |z| is capped near (n-1)/sqrt(n). Measured: max |z| of 1.97/1.94/2.04 at
    n=6, against 18.5/29.2/42.5 for the same sample held out. "0% beyond |z|>2"
    from a self-included run is arithmetic, not evidence.
    """
    from pathlib import Path

    from krewlyzer.core.wps_pon import apply_wps_pon

    frames = [_frame(s) for s in range(6)]
    pon, directory, paths = _pon_from(frames)
    apply_wps_pon(
        Path(paths[0]), pon, output_base=directory / "self", output_format="parquet"
    )
    included = pd.read_parquet(directory / "self.parquet")
    z = pd.to_numeric(included["wps_shape_corr_z"], errors="coerce").dropna()
    cap = (len(frames) - 1) / math.sqrt(len(frames))
    assert (
        z.abs().max() <= cap + 0.1
    ), f"self-included |z| reached {z.abs().max():.2f}, above the {cap:.2f} cap"
