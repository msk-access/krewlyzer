"""Per-motif z-scores from the PON's k-mer baseline.

`mds_baseline` carries 625 k-mer means and standard deviations — every 4-mer
over the alphabet ACGTN — and nothing read them. `EndMotif.parquet` shipped
`Motif, Frequency` only, so a shift in one motif was invisible unless it moved
the whole-sample MDS summary enough to notice.

The join is on the motif string, never on position: the baseline has 625 keys
and the output has 256 (ACGT only), so the two are not aligned and assuming
they were would silently pair unrelated motifs.

## The normalisation the join has to correct for

Sample frequencies sum to 1.0 across the 256 ACGT motifs. The baseline's
expectations for those same 256 sum to **0.972** — the missing 2.79% sits in
the N-containing k-mers the output does not report.

Comparing them directly biases every z upward: measured on a real sample the
naive z has median **+0.37** with 20% beyond \\|z\\|>2, which is a normalisation
artefact and not a motif signal. Restricting both sides to the shared motifs
and renormalising brings the median to −0.21.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from krewlyzer.pon.model import zscore_or_nan

from .output_utils import read_exact_table, write_table

logger = logging.getLogger("core.motif_pon")


def apply_motif_pon(
    motif_path: Path,
    pon,
    output_base: Optional[Path] = None,
    output_format: str = "tsv",
    compress: bool = False,
    baseline_attr: str = "mds_baseline",
) -> int:
    """Add ``frequency_z`` to an EndMotif or BreakPointMotif table.

    ``baseline_attr`` names the k-mer block to score against, and there are
    four -- one per table. Until 0.9.0 all four tables were scored against
    ``mds_baseline``, which is an **end-motif, genome-wide** baseline, so three
    of the four measured the offset between two distributions rather than the
    sample's departure from healthy.

    Breakpoint motifs are the sharper case. An end motif is the 4-mer at the
    fragment's 5' terminus -- what the nuclease left. A breakpoint motif spans
    the cut site and includes reference bases *not present in the fragment*.
    On one sample the two frequency vectors correlate 0.696. Measured against
    the end-motif baseline, ``BreakPointMotif`` came out at median |z| 5.85 on
    XS1 and 11.25 on XS2, with 70 and 136 of 256 motifs beyond |z| = 10; a
    correct baseline gives a median near 0.67.

    Returns the number of motifs scored.
    """
    frame = read_exact_table(motif_path)
    if frame is None:
        logger.warning(f"Motif output not found for PON scoring: {motif_path}")
        return 0
    if "Motif" not in frame.columns or "Frequency" not in frame.columns:
        logger.warning(f"Motif output lacks Motif/Frequency columns: {motif_path}")
        return 0

    baseline = getattr(pon, baseline_attr, None)
    expected = getattr(baseline, "kmer_expected", None) if baseline else None
    stds = getattr(baseline, "kmer_std", None) if baseline else None
    if not expected or not stds:
        # No column, not a NaN one -- `test_no_baseline_means_no_column_rather
        # _than_a_fabricated_one` pins this, and the asymmetry is informative:
        # a PON built before 0.9.0 has the end-motif blocks and not the
        # breakpoint ones, so EndMotif gets `frequency_z` and BreakPointMotif
        # visibly does not. That reads as "no baseline for this table", which
        # is exactly what happened.
        logger.info(
            f"PON has no '{baseline_attr}'; {motif_path.name} keeps its raw "
            "frequencies and gets no z-score. Rebuild the PON with build-pon."
        )
        return 0

    motifs = [str(m) for m in frame["Motif"]]
    shared = [m for m in motifs if m in expected]
    if not shared:
        logger.warning(
            f"No motif in {motif_path.name} matches the PON's k-mer baseline "
            f"({len(expected)} keys). Is this PON for the right assay?"
        )
        frame["frequency_z"] = np.nan
        write_table(
            frame,
            output_base or motif_path.with_suffix(""),
            output_format=output_format,
            compress=compress,
        )
        return 0

    # Both sides restricted to the shared motifs, then put on the same total.
    # Scaling the sample by the baseline's partial mass is algebraically the
    # same as renormalising the baseline and its sigma, and keeps the sigma
    # untouched, which is easier to reason about.
    baseline_mass = float(sum(expected[m] for m in shared))
    if not np.isfinite(baseline_mass) or baseline_mass <= 0:
        logger.warning("PON k-mer expectations sum to zero over the shared motifs")
        return 0

    z_scores = []
    n_absent = 0
    for motif, frequency in zip(motifs, frame["Frequency"]):
        if motif not in expected:
            n_absent += 1
            z_scores.append(float("nan"))
            continue
        z_scores.append(
            zscore_or_nan(
                float(frequency) * baseline_mass, expected[motif], stds.get(motif)
            )
        )

    frame["frequency_z"] = z_scores
    write_table(
        frame,
        output_base or motif_path.with_suffix(""),
        output_format=output_format,
        compress=compress,
    )

    n_scored = int(np.isfinite(np.asarray(z_scores, dtype=float)).sum())
    logger.info(
        f"Motif PON: {n_scored}/{len(frame)} motifs scored, baseline mass "
        f"{baseline_mass:.4f} over {len(shared)} shared ({motif_path.name})"
    )
    if n_absent:
        logger.warning(
            f"Motif PON: {n_absent}/{len(frame)} motifs are absent from the "
            f"baseline's {len(expected)} k-mers and get no z-score."
        )
    if baseline_mass < 0.5:
        logger.warning(
            f"Motif PON: the baseline covers only {baseline_mass:.1%} of the "
            "expected mass over the shared motifs. The renormalisation is "
            "doing a lot of work; check the PON matches this motif length."
        )
    return n_scored
