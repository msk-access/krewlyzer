"""What each output table measures, and which way cancer moves it.

Distilled from the output EDA notebook, whose most useful observation was this:

    Direction matters and differs per axis — this is the single most common
    misreading.

That is not hypothetical here. ``region-mds.md`` asserted the *opposite* of its
own threshold table for a year, and nothing caught it because no test and no
report ever stated a direction. So direction is data now, next to the contract
that already declares the columns, and ``describe-output`` renders it.

## Why this file rather than the docs

``docs/reference/output-files.md`` is narrative: what a table is for, how to use
it, worked examples. This is the one-line form a reader needs *while looking at
the file* — and it has to be machine-readable so the report can print it beside
the measured shape.

The risk is the obvious one: two descriptions that drift. Guarded by
``tests/unit/test_meaning_registry.py``, which fails if a contracted table has
no entry here or an entry names a table that no longer exists. Prose depth stays
in the docs; only the one-liner and the direction live here.

## On thresholds

None are recorded. Every numeric band the notebook examined turned out to be a
display default at best, and refuted at worst — the ATAC/TFBS "4–6 bits normal,
>8 tumour" range flags a perfectly healthy N(167,30) distribution as abnormal.
Directions are robust; magnitudes are cohort-specific. Recording a threshold
here would give it an authority it has not earned.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional


@dataclass(frozen=True)
class Meaning:
    """One line of interpretation for an output table."""

    #: What the table quantifies, in a sentence a reader can act on.
    measures: str
    #: Which way the tumour signal moves, or ``None`` where the table is
    #: descriptive rather than directional (identifiers, raw counts, metadata).
    cancer_direction: Optional[str] = None
    #: Anything that would cause a confident misreading. Kept short.
    caveat: Optional[str] = None


_SIZE_DIRECTION = "higher short-fragment fraction"

MEANINGS: Dict[str, Meaning] = {
    # -- fragment size ------------------------------------------------------
    ".FSD.parquet": Meaning(
        "Fragment-size histogram, 5 bp bins over 65–1000 bp, per chromosome arm. "
        "Summed across arms this is the sample's genome-wide size density — the "
        "most direct picture of the whole thesis.",
        "mass shifts from the ~166 bp mononucleosomal mode toward ~145 bp",
    ),
    ".FSD.ontarget.parquet": Meaning(
        "As FSD, restricted to captured regions. Panel capture biases fragment "
        "size, so on-target and off-target are not interchangeable.",
        "as FSD, but compare on-target to on-target only",
    ),
    ".FSR.parquet": Meaning(
        "Short:long ratio per window, with each side already divided by its PON "
        "mean — so the ratio is healthy-relative and cancels batch effects. "
        "`short_long_log2` is the ML-ready signed form, ~0 in healthy plasma.",
        "higher `short_long_log2`",
    ),
    ".FSR.ontarget.parquet": Meaning(
        "As FSR, over captured regions.", "higher `short_long_log2`"
    ),
    ".FSC.parquet": Meaning(
        "GC-corrected coverage split into six non-overlapping size channels per "
        "genomic bin. Ties fragment size to *location*, so focal signal survives "
        "a genome average that would hide it.",
        _SIZE_DIRECTION + " in affected bins",
        "Not PON-normalised — the healthy expectation is not zero. Read the "
        "spread across bins, not the offset.",
    ),
    ".FSC.ontarget.parquet": Meaning("As FSC, over captured regions.", _SIZE_DIRECTION),
    ".FSC.gene.parquet": Meaning(
        "The six size channels aggregated per gene, plus an RPKM-like "
        "`normalized_depth` usable as a copy-number proxy.",
        _SIZE_DIRECTION + " at affected genes",
    ),
    ".FSC.regions.parquet": Meaning(
        "The same channels per exon or capture tile — the finest FSC resolution. "
        "Carries `strand`, `is_e1` and `is_alt_e1` so promoter-proximal regions "
        "can be selected without re-deriving them.",
        _SIZE_DIRECTION + " at affected regions",
    ),
    # -- motif --------------------------------------------------------------
    ".EndMotif.parquet": Meaning(
        "Frequency of each of the 256 4-mers at fragment 5′ ends. Nuclease "
        "cutting preference is sequence-specific, and tumour cfDNA is cut "
        "differently.",
        "a few motifs dominate — the distribution narrows",
    ),
    ".EndMotif.ontarget.parquet": Meaning(
        "As EndMotif, over captured regions.", "distribution narrows"
    ),
    ".BreakPointMotif.parquet": Meaning(
        "4-mer frequencies at the breakpoint rather than the fragment end — the "
        "sequence context spanning the cut site.",
        "distribution narrows",
    ),
    ".BreakPointMotif.ontarget.parquet": Meaning(
        "As BreakPointMotif, over captured regions.", "distribution narrows"
    ),
    ".EndMotif1mer.parquet": Meaning(
        "Single-base composition at fragment ends. DNASE1L3 leaves 5′ overhangs; "
        "filling them in during library prep shifts the terminal base. A cheap "
        "complement to MDS — MDS measures diversity, this measures the chemistry "
        "of the cut end.",
        "higher C-end fraction",
    ),
    ".MDS.parquet": Meaning(
        "Motif Diversity Score: Shannon entropy of the 256 4-mer end motifs, "
        "normalised by log2(256) so it lands in [0, 1].",
        "**lower** MDS — stereotyped cutting",
        "The direction is the opposite of intuition and was documented "
        "backwards for a year. Lower is the abnormal end.",
    ),
    ".MDS.ontarget.parquet": Meaning("As MDS, over captured regions.", "**lower** MDS"),
    ".MDS.gene.parquet": Meaning(
        "MDS per gene, plus `mds_e1` for the first exon specifically — the "
        "promoter's nucleosome-depleted region, where accessibility differences "
        "are largest and early-cancer signal is strongest (Helzer 2025).",
        "**lower** `mds_e1` at a driver locus",
        "`mds_e1` is NaN where the gene has no captured first exon, which on a "
        "targeted panel is most genes. NaN is not zero.",
    ),
    ".MDS.exon.parquet": Meaning(
        "MDS per exon or capture tile — localises aberrant cutting.",
        "**lower** MDS",
    ),
    # -- orientation and accessibility --------------------------------------
    ".OCF.ontarget.parquet": Meaning(
        "Orientation-aware fragmentation at tissue-specific open chromatin. "
        "Fragments ending upstream versus downstream of a region are counted "
        "separately; a tissue actively shedding DNA produces a phased excess.",
        "the shedding tissue rises; a **fall** in T-cell OCF is equally "
        "informative, since it reflects dilution of normal haematopoietic "
        "background by tumour DNA",
    ),
    ".OCF.offtarget.parquet": Meaning(
        "As OCF, outside captured regions — the unbiased view for panel data.",
        "as OCF",
    ),
    ".TFBS.parquet": Meaning(
        "Shannon entropy of the fragment-size distribution at the binding sites "
        "of 808 transcription factors. A proxy for occupancy and local "
        "accessibility, and the highest-dimensional single feature emitted.",
        "shifts either way depending on the factor — interpret per TF, and "
        "prefer z-scores",
        "Absolute entropy bands are refuted: a healthy-like N(167,30) sample "
        "already exceeds the documented 'normal' range. Use z-scores only.",
    ),
    ".TFBS.ontarget.parquet": Meaning(
        "As TFBS, over captured regions.", "per-factor; z-scores only"
    ),
    ".ATAC.parquet": Meaning(
        "The same size-entropy computation over each of 23 TCGA cancer types' "
        "accessible peaks — a panel-compatible accessibility readout, orthogonal "
        "to OCF's orientation signal.",
        "shifts at the peaks of the matching cancer type; z-scores only",
        "Absolute entropy bands are refuted — see TFBS.",
    ),
    ".ATAC.ontarget.parquet": Meaning(
        "As ATAC, over captured regions.", "z-scores only"
    ),
    # -- nucleosome positioning ---------------------------------------------
    ".WPS.parquet": Meaning(
        "Windowed protection score around TSS and CTCF anchors: how strongly "
        "each position is protected from cutting, which traces nucleosome "
        "placement. Per-position vectors, not scalars.",
        "protection at active promoters weakens as positioning is disrupted",
    ),
    ".WPS.panel.parquet": Meaning("As WPS, over the panel's anchors.", "as WPS"),
    ".WPS_background.parquet": Meaning(
        "WPS stacked over Alu elements, giving a genome-wide chromatin quality "
        "readout and the nucleosome repeat length from its periodicity.",
        "NRL drifts from the ~190 bp healthy repeat length",
        "Check `nrl_at_band_limit`: where it is set, `nrl_bp` is the edge of the "
        "search band rather than a measurement — no peak was found.",
    ),
    # -- provenance ----------------------------------------------------------
    ".metadata.parquet": Meaning(
        "Run provenance and fragment totals. Also the **completion marker** — a "
        "sample without this file is dropped from the cohort silently.",
        None,
    ),
}
