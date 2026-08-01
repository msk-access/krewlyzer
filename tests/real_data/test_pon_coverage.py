"""Is every PON baseline actually consumed, and is the result calibrated?

Two questions a synthetic fixture cannot answer.

**Consumption.** `wps_baseline` is ~128k anchors, roughly 90% of every PON
file, and on the 0.8.3 corpus its only consumer was a log line appending
`"WPS"` to a list of available components. Three more blocks were the same.
Reading the code shows a loader; reading a cohort shows whether anything
downstream carries the result.

**Calibration.** Held-out healthy samples should score near N(0,1). At n=21
that only catches gross miscalibration — but gross is exactly what a
fabricated baseline produces: the shipped `wps_background` block would have
given every sample the identical z.

These are expected to *fail* until Phase B wires the missing scorers, so the
consumption test is marked xfail(strict) — when the wiring lands it turns red
and forces the marker off.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytestmark = pytest.mark.integration

#: PON block -> (table it should score, column it should produce).
#:
#: Straight from `pon/build.py`. A block that produces no column anywhere is
#: computed, stored and shipped for nothing.
BLOCK_CONSUMERS = {
    "wps_baseline": (".WPS.parquet", "wps_nuc_z"),
    "wps_background": (".WPS_background.parquet", "nrl_z"),
    "mds_baseline.kmer": (".EndMotif.parquet", "frequency_z"),
    "region_mds_exon": (".MDS.exon.parquet", "mds_z"),
}


@pytest.mark.xfail(
    strict=True,
    reason="Phase B wires these; strict so it fails loudly once they land",
)
@pytest.mark.parametrize("block", sorted(BLOCK_CONSUMERS))
def test_every_pon_block_produces_a_column(corpus_tables, block):
    suffix, column = BLOCK_CONSUMERS[block]
    frames = corpus_tables.get(suffix)
    if not frames:
        pytest.skip(f"{suffix} absent from this cohort")
    carrying = sum(1 for f in frames if column in f.columns)
    assert carrying > 0, (
        f"the `{block}` baseline is built into every PON, but no sample's "
        f"{suffix} carries `{column}`. The block is shipped and read by nothing."
    )


# ---------------------------------------------------------------------------
# calibration
# ---------------------------------------------------------------------------

#: Table and column of every z-score the cohort currently carries.
Z_COLUMNS = [
    (".MDS.gene.parquet", "mds_z"),
    (".TFBS.parquet", "z_score"),
    (".ATAC.parquet", "z_score"),
    (".OCF.offtarget.parquet", "ocf_z"),
]


@pytest.mark.parametrize("suffix,column", Z_COLUMNS)
def test_a_z_score_is_not_a_constant(corpus_tables, suffix, column):
    """The single check that would have caught the fabricated NRL baseline.

    A baseline of `mean=167, std=5` applied to a metric pinned at 150 gives
    every sample `z = -3.4`. Plausible, extreme enough to look meaningful, and
    identical everywhere.
    """
    frames = corpus_tables.get(suffix)
    if not frames:
        pytest.skip(f"{suffix} absent from this cohort")
    means = [
        float(pd.to_numeric(f[column], errors="coerce").mean(skipna=True))
        for f in frames
        if column in f.columns
    ]
    finite = [m for m in means if np.isfinite(m)]
    if len(finite) < 2:
        pytest.skip(f"{suffix}.{column}: fewer than 2 finite sample means")
    assert np.std(finite) > 0, (
        f"{suffix}.{column} is identical across all {len(finite)} samples "
        f"({finite[0]:.4f}). A constant z is a fabricated baseline."
    )


@pytest.mark.parametrize("suffix,column", Z_COLUMNS)
def test_a_z_score_is_not_wildly_miscalibrated(corpus_tables, suffix, column):
    """A weak bound, and deliberately so.

    These are tumour samples, not the healthy cohort the PON was fitted on, so
    a shifted mean is the expected biology, not a defect. What this rejects is
    the scale being wrong by an order of magnitude -- a sigma floored to 0.01,
    or a baseline fitted on different units.
    """
    frames = corpus_tables.get(suffix)
    if not frames:
        pytest.skip(f"{suffix} absent from this cohort")
    values = np.concatenate(
        [
            pd.to_numeric(f[column], errors="coerce").dropna().to_numpy()
            for f in frames
            if column in f.columns
        ]
    )
    values = values[np.isfinite(values)]
    if values.size < 10:
        pytest.skip(f"{suffix}.{column}: too few finite values")
    extreme = float(np.mean(np.abs(values) > 10))
    assert extreme < 0.5, (
        f"{100 * extreme:.0f}% of {suffix}.{column} exceeds |z| > 10. That is a "
        "scale error -- a floored sigma or a mismatched baseline -- not biology."
    )
