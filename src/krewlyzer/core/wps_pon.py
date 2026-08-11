"""Apply the PON to WPS output — the largest baseline, previously unread.

`wps_baseline` is ~128k anchors of 200-element mean and sigma vectors, roughly
90% of every PON file. Until 0.9.0 its only consumer was a log line appending
``"WPS"`` to a list of available components: the baseline was built, stored and
shipped, and nothing downstream carried a single value derived from it.

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

**The arithmetic lives in Rust.** ``rust/src/wps.rs::apply_pon_zscore`` does the
work; this module is the call. It was written in Python first because three of
its decisions changed under measurement and would not have survived being
written in the faster language first -- the ``log1p``, the Fisher transform, and
the absent phase baseline all came from looking at real cohort numbers. Once
settled, it belonged on the Rust side of the boundary in
``.agents/rules/architecture.md``: 89,034 anchors, a +/-30 lag search each, about
5.4M correlations. The Python reference survives as the equivalence oracle in
``tests/unit/test_rust_python_equivalence.py``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger("core.wps_pon")

#: How far the phase search looks, in positions.
#:
#: Mirrors ``PHASE_MAX_LAG`` in ``rust/src/wps.rs``, which is the copy that runs.
#: Asserted equal in ``tests/unit/test_wps_pon.py`` -- a constant kept in two
#: languages drifts unless something fails when it does (invariant #5).
#:
#: The window was chosen by measurement: on the healthy cohort it was tuned
#: against, the search terminates on its own edge for 3.3% of anchors at +/-20
#: and 1.8% at +/-30. That 1.8% is the *matched* figure and is not a promise --
#: scoring a 0.8.3-era WPS output against a 0.9.0 PON measured 12.5%, which is
#: what a poorly-matched sample and baseline look like rather than a defect.
#: `wps_phase_at_search_limit` exists so the reader can tell which they have.
PHASE_MAX_LAG = 30


def apply_wps_pon(
    wps_path: Path,
    pon_parquet_path: Path,
    output_base: Optional[Path] = None,
    column: str = "wps_nuc",
    baseline_table: str = "wps_baseline",
) -> int:
    """Add PON-derived columns to a WPS table. Returns anchors scored.

    ``baseline_table`` selects which vector baseline to score against:
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

    Parquet, always -- not whatever ``--output-format`` says. WPS is the one
    family the output contract pins to Parquet
    (``docs/reference/output-files.md``), because the 200-element vectors are
    unusable as text: 928 MB against 144 MB for the same 89k anchors. Honouring
    the flag here previously wrote the z-scores to a *second* file, ``.WPS.tsv``,
    and left the parquet as the raw profile -- and downstream reads Parquet
    only (invariant #2), so the scoring never reached the product.

    Returns 0 without writing when the PON has no matching vector baseline, so
    the caller keeps the raw table rather than losing it to a half-written one.
    """
    from krewlyzer import _core

    if not wps_path.exists():
        logger.warning(f"WPS output not found for PON scoring: {wps_path}")
        return 0
    if not Path(pon_parquet_path).exists():
        logger.warning(f"PON not found for WPS scoring: {pon_parquet_path}")
        return 0

    # `.parquet` explicitly, not `with_suffix`: `output_base` is a stem, and the
    # in-place case (base == the input's stem) must land on the same file the
    # Rust step wrote, which is the file downstream reads.
    output_path = (output_base or wps_path.with_suffix("")).with_suffix(".parquet")

    n_scored = _core.wps.apply_pon_zscore(
        str(wps_path),
        str(pon_parquet_path),
        str(output_path),
        baseline_table,
        column,
    )

    # Rust returns 0 *and writes nothing* when the PON carries no matching
    # baseline -- a PON built before the block existed, or a panel PON asked
    # for the genome-wide table. Say so here rather than letting a caller
    # discover the missing columns downstream.
    if n_scored == 0:
        logger.warning(
            f"WPS PON: no anchors scored from '{baseline_table}'. "
            f"{wps_path.name} keeps its raw profile and gets no z-scores. "
            "Rebuild the PON with build-pon."
        )
    return n_scored
