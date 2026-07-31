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
from typing import Dict, Optional, Tuple


@dataclass(frozen=True)
class Meaning:
    """How to read one output table.

    Four fields, because a reader opening an unfamiliar file needs four
    different things and conflating them into one paragraph helps nobody:

    ``measures``  the one-liner, for index and summary tables
    ``why``       the biology — why this measurement exists at all
    ``how``       how the number is arrived at
    ``what``      what to actually look at, once you are staring at the numbers

    ``what`` deliberately never contains a threshold. It can say *the right
    tail matters more than the median*; it cannot say *above 1.3 is positive*.
    Every numeric band in this project that was treated as a cut-off turned out
    to be a display default or refuted outright, and a number printed beside a
    column acquires an authority it has not earned.
    """

    #: What the table quantifies, in a sentence a reader can act on.
    measures: str
    #: Which way the tumour signal moves, or ``None`` where the table is
    #: descriptive rather than directional (identifiers, raw counts, metadata).
    cancer_direction: Optional[str] = None
    #: Anything that would cause a confident misreading. Kept short.
    caveat: Optional[str] = None
    #: The biology. Why anyone measures this.
    why: Optional[str] = None
    #: How the number is computed.
    how: Optional[str] = None
    #: What to look at in the values. Never a threshold.
    what: Optional[str] = None
    #: Columns that exist, or mean anything, only when a PON was supplied.
    #: Without one they are absent or uninterpretable, and the report says so
    #: rather than letting a reader take a raw value for a normalised one.
    pon_columns: Tuple[str, ...] = ()


_SIZE_DIRECTION = "higher short-fragment fraction"

MEANINGS: Dict[str, Meaning] = {
    # -- fragment size ------------------------------------------------------
    ".FSD.parquet": Meaning(
        "Fragment-size histogram, 5 bp bins over 65–1000 bp, per chromosome arm. "
        "Summed across arms this is the sample's genome-wide size density — the "
        "most direct picture of the whole thesis.",
        "mass shifts from the ~166 bp mononucleosomal mode toward ~145 bp",
        why=(
            "Nucleosomes protect ~147 bp of DNA from nuclease digestion, so plasma cfDNA arrives in a sharply peaked size distribution. Tumour cells package and shed DNA differently, and the whole field rests on that difference being measurable."
        ),
        how=(
            "Every fragment is binned by length into 5 bp bins from 65 to 1000 bp, counted per chromosome arm."
        ),
        what=(
            "Look at where the mode sits and how heavy the sub-150 bp shoulder is. Summing across arms gives the genome-wide density; comparing arms surfaces regional differences."
        ),
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
        why=(
            "A ratio of short to long fragments is the most direct summary of the size shift, and dividing each side by a healthy baseline first cancels the batch and GC effects that otherwise dominate."
        ),
        how=(
            "Per window: short and long counts are each divided by their PON mean, then the ratio is taken. `short_long_log2` is the signed log of that, ~0 in healthy plasma."
        ),
        what=(
            "The median across windows is the genome-wide summary. The spread and the right tail matter more — focal high-burden regions live in that tail even when the median looks ordinary."
        ),
        pon_columns=("short_norm", "long_norm", "short_long_ratio", "short_long_log2"),
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
        why=(
            "A genome average hides heterogeneity. Tying fragment size to location lets a focal signal survive that would otherwise be diluted across the genome."
        ),
        how=(
            "Fragments are GC-corrected, then counted into six non-overlapping size channels per genomic bin. The channels partition the total, so the six ratios sum to 1."
        ),
        what=(
            "Compare channel ratios between bins rather than reading any single bin. Unlike FSR these are raw counts, not PON-normalised, so the healthy expectation is not zero."
        ),
        pon_columns=(
            "ultra_short_log2",
            "core_short_log2",
            "mono_nucl_log2",
            "di_nucl_log2",
            "long_log2",
            "ultra_long_log2",
        ),
    ),
    ".FSC.ontarget.parquet": Meaning("As FSC, over captured regions.", _SIZE_DIRECTION),
    ".FSC.gene.parquet": Meaning(
        "The six size channels aggregated per gene, plus an RPKM-like "
        "`normalized_depth` usable as a copy-number proxy.",
        _SIZE_DIRECTION + " at affected genes",
        why=(
            "On a targeted panel the gene is the natural unit: it is what the assay was designed around and what a clinician reasons about."
        ),
        how=(
            "The six channels summed over every captured region belonging to a gene, plus `normalized_depth`, an RPKM-like value usable as a copy-number proxy."
        ),
        what=(
            "`normalized_depth` compares coverage between genes in one sample. The channel ratios say whether a gene's fragments are unusually short."
        ),
    ),
    ".FSC.regions.parquet": Meaning(
        "The same channels per exon or capture tile — the finest FSC resolution. "
        "Carries `strand`, `is_e1` and `is_alt_e1` so promoter-proximal regions "
        "can be selected without re-deriving them.",
        _SIZE_DIRECTION + " at affected regions",
        why=(
            "Exon-level resolution localises a signal to part of a gene — a hotspot exon can differ from the gene average."
        ),
        how=(
            "The same six channels per exon or capture tile. Carries `strand`, `is_e1` and `is_alt_e1` from the gene BED so promoter-proximal regions can be selected without re-deriving them."
        ),
        what=(
            "`is_e1` marks the canonical transcript's first exon; `is_alt_e1` marks another basic protein-coding transcript's first exon. Both are real transcription starts. A region with neither is an internal exon."
        ),
    ),
    # -- motif --------------------------------------------------------------
    ".EndMotif.parquet": Meaning(
        "Frequency of each of the 256 4-mers at fragment 5′ ends. Nuclease "
        "cutting preference is sequence-specific, and tumour cfDNA is cut "
        "differently.",
        "a few motifs dominate — the distribution narrows",
        why=(
            "Nucleases cut with sequence preference, so the bases at fragment ends record which enzymes did the cutting. Tumour cfDNA carries a different signature."
        ),
        how=(
            "The 4-mer at each fragment's 5′ end is counted across all 256 possibilities, then normalised to frequencies summing to 1."
        ),
        what=(
            "Read the shape of the ranked distribution rather than any one motif. A broad, flat profile is healthy; mass concentrating on a few motifs is what MDS quantifies as lower diversity."
        ),
    ),
    ".EndMotif.ontarget.parquet": Meaning(
        "As EndMotif, over captured regions.", "distribution narrows"
    ),
    ".BreakPointMotif.parquet": Meaning(
        "4-mer frequencies at the breakpoint rather than the fragment end — the "
        "sequence context spanning the cut site.",
        "distribution narrows",
        why=(
            "The sequence spanning the cut site, rather than the fragment end, captures context the end motif alone misses."
        ),
        how=("As EndMotif, but the 4-mer straddles the breakpoint."),
        what=(
            "Compare its shape against EndMotif; a divergence between them points at the cutting chemistry rather than at fragment selection."
        ),
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
        why=(
            "DNASE1L3 leaves single-stranded 5′ overhangs. Filling those in during library prep writes templated bases at the end, so the single-base composition records the chemistry of the cut — orthogonal to how diverse the 4-mers are."
        ),
        how=("The single 5′ base of each fragment, counted over A/C/G/T."),
        what=(
            "The C fraction is the quantity of interest. It is a cheap complement to MDS: MDS measures diversity, this measures the cut chemistry."
        ),
    ),
    ".MDS.parquet": Meaning(
        "Motif Diversity Score: Shannon entropy of the 256 4-mer end motifs, "
        "normalised by log2(256) so it lands in [0, 1].",
        "**lower** MDS — stereotyped cutting",
        "The direction is the opposite of intuition and was documented "
        "backwards for a year. Lower is the abnormal end.",
        why=(
            "If cutting becomes stereotyped, the *diversity* of end motifs falls even when no single motif stands out. Entropy captures that in one number."
        ),
        how=(
            "Shannon entropy of the 256 4-mer frequencies, divided by log2(256) so the result lands in [0, 1]."
        ),
        what=(
            "Lower is the abnormal direction. Compare against a PON via `mds_z` rather than reading the raw value, which varies with depth and library prep."
        ),
        pon_columns=("mds_z",),
    ),
    ".MDS.ontarget.parquet": Meaning("As MDS, over captured regions.", "**lower** MDS"),
    ".MDS.gene.parquet": Meaning(
        "MDS per gene, plus `mds_e1` for the first exon specifically — the "
        "promoter's nucleosome-depleted region, where accessibility differences "
        "are largest and early-cancer signal is strongest (Helzer 2025).",
        "**lower** `mds_e1` at a driver locus",
        "`mds_e1` is NaN where the gene has no captured first exon, which on a "
        "targeted panel is most genes. NaN is not zero.",
        why=(
            "A global score says cutting is abnormal; per-gene says where. Promoter-proximal first exons sit in nucleosome-depleted regions, where accessibility differences between tumour and normal are largest (Helzer 2025)."
        ),
        how=(
            "The same entropy computation over only the fragments overlapping each gene's regions, plus `mds_e1` for the first exon specifically."
        ),
        what=(
            "`mds_e1` at a known driver is the targeted question. NaN means the gene has no captured first exon — on a hotspot panel that is most genes, and NaN is not zero."
        ),
        pon_columns=("mds_z", "mds_e1_z"),
    ),
    ".MDS.exon.parquet": Meaning(
        "MDS per exon or capture tile — localises aberrant cutting.",
        "**lower** MDS",
        why=("The finest localisation available: which exon, not which gene."),
        how=("Entropy over the fragments overlapping each exon or capture tile."),
        what=(
            "Sparse by nature — an exon with few fragments has an unstable entropy. Read `n_fragments` alongside `mds`."
        ),
    ),
    # -- orientation and accessibility --------------------------------------
    ".OCF.ontarget.parquet": Meaning(
        "Orientation-aware fragmentation at tissue-specific open chromatin. "
        "Fragments ending upstream versus downstream of a region are counted "
        "separately; a tissue actively shedding DNA produces a phased excess.",
        "the shedding tissue rises; a **fall** in T-cell OCF is equally "
        "informative, since it reflects dilution of normal haematopoietic "
        "background by tumour DNA",
        why=(
            "At a tissue-specific open chromatin region, DNA from that tissue fragments in a characteristic phased pattern. The orientation of fragment ends therefore identifies which tissue is shedding."
        ),
        how=(
            "Fragments ending upstream and downstream of each region are counted separately and the asymmetry is normalised per tissue."
        ),
        what=(
            "A rise means that tissue is shedding. A fall in T-cell is equally informative — it reflects dilution of the normal haematopoietic background by tumour DNA."
        ),
        pon_columns=("ocf_z",),
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
        why=(
            "Fragment-size diversity at a transcription factor's binding sites is a proxy for that factor's occupancy and the local accessibility. Cancer rewires TF programmes, and this reads that out without needing accessibility data for the sample."
        ),
        how=(
            "Shannon entropy of the fragment-size distribution at the binding sites of each of 808 factors (GTRD)."
        ),
        what=(
            "Interpretable per factor — lineage factors like HNF4A for hepatocyte or SPI1 for haematopoietic, E2F for proliferation. Use the z-scores; absolute entropy bands are refuted."
        ),
        pon_columns=("z_score",),
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
        why=(
            "Each cancer type has a characteristic set of accessible enhancers. Aberrant nucleosome positioning at those loci changes local fragment-size diversity, giving a panel-compatible accessibility readout orthogonal to OCF's orientation signal."
        ),
        how=(
            "The same entropy computation at each of 23 TCGA cancer types' ATAC-seq peak sets."
        ),
        what=(
            "The relative ranking across types is the signal. Absolute entropy is not interpretable — a healthy N(167,30) distribution already exceeds the documented 'normal' band."
        ),
        pon_columns=("z_score",),
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
        why=(
            "A positioned nucleosome protects the DNA beneath it, so the density of fragment endpoints across a region traces where nucleosomes sit. At a TSS that positioning reflects whether the gene is active."
        ),
        how=(
            "For each position, fragments fully spanning a window score positive and fragments ending inside it score negative, summed across TSS and CTCF anchors. Per-position vectors, not scalars."
        ),
        what=(
            "Reduce the vectors with a mean or max before comparing; consumers use list_avg/list_max. The depth and phasing of the dip at the anchor are the features."
        ),
    ),
    ".WPS.panel.parquet": Meaning("As WPS, over the panel's anchors.", "as WPS"),
    ".WPS_background.parquet": Meaning(
        "WPS stacked over Alu elements, giving a genome-wide chromatin quality "
        "readout and the nucleosome repeat length from its periodicity.",
        "NRL drifts from the ~190 bp healthy repeat length",
        "Check `nrl_at_band_limit`: where it is set, `nrl_bp` is the edge of the "
        "search band rather than a measurement — no peak was found.",
        why=(
            "Stacking WPS over the ~1 million Alu elements averages away locus-specific noise and leaves a genome-wide chromatin quality readout, including the nucleosome repeat length."
        ),
        how=(
            "WPS profiles are stacked over Alu consensus positions, then the repeat length is read from the periodicity of the stacked wave by FFT."
        ),
        what=(
            "Check `nrl_at_band_limit` first. Where it is set, `nrl_bp` is the edge of the search band rather than a measurement — no periodic peak was found, and the value is not a length."
        ),
    ),
    # -- provenance ----------------------------------------------------------
    ".metadata.parquet": Meaning(
        "Run provenance and fragment totals. Also the **completion marker** — a "
        "sample without this file is dropped from the cohort silently.",
        None,
        why=(
            "Provenance: which inputs produced this directory, and how many fragments survived filtering. Also the completion marker the downstream consumer looks for."
        ),
        how=("Written at the end of a run, once every other output has been produced."),
        what=(
            "`total_fragments` bounds everything else — a table computed from few fragments is noisy regardless of what it says. Its absence means the sample is silently dropped downstream."
        ),
    ),
    # -- optional outputs -----------------------------------------------------
    #
    # Not in ``CONTRACT``: each needs an input the pipeline cannot assume, so
    # gating on them would fail every run that legitimately lacks one. They are
    # still described, because a file sitting in the directory needs an entry
    # here whether or not anything downstream reads it -- a section with a shape
    # and no interpretation was the exact complaint that put `describe-output`
    # on `NOT_CONSUMED` in the first place.
    ".mFSD.parquet": Meaning(
        "Fragment sizes carrying the mutant allele versus the reference allele, "
        "per variant. The only output that separates tumour from background "
        "within a single sample rather than against a cohort.",
        "**shorter** ALT than REF fragments",
        "`ALT_Count` can be zero for a variant that was called but is not "
        "observed here. Nothing derived from it — sizes, LLR, KS — is a "
        "measurement in that case, and `ALT_MeanSize` reads 0, not the "
        "~166 bp a plausible-looking value would suggest.",
        why=(
            "Tumour-derived fragments are shorter than the healthy background. Splitting by allele makes that comparison internal: the reference fragments at the same locus are the control, so no panel of normals is needed."
        ),
        how=(
            "Fragments overlapping each variant position are assigned to REF or ALT by the base they carry, and the two size distributions are compared by mean, by log-likelihood ratio against a size model, and by a KS test."
        ),
        what=(
            "Read `ALT_Count` before anything else — a handful of fragments cannot support a distribution comparison, and the KS statistic is computed from as few as two. Then the sign of `ALT_MeanSize - REF_MeanSize`, which is the whole measurement."
        ),
    ),
    ".UXM.parquet": Meaning(
        "Fragment-level methylation: each read classified unmethylated, mixed "
        "or methylated, then summarised per marker region.",
        "the **U** fraction rises at markers for the tissue of origin",
        "Needs a methylation-aware BAM (`--bisulfite-bam`). Absent from an "
        "ordinary run, which is not a failure.",
        why=(
            "Methylation is the sharpest tissue-of-origin signal in cfDNA, and calling it per fragment rather than per CpG keeps the read as the unit of evidence -- a fragment is either from that tissue or it is not."
        ),
        how=(
            "Each read's CpGs are read from its bisulfite tag and the read is classified U, X or M by the fraction methylated; the three counts are then summed per marker region."
        ),
        what=(
            "`U + X + M` sums to 1 by construction, so only two of the three are independent. Compare the U fraction at markers for a candidate tissue against the same markers in a normal."
        ),
    ),
    ".OCF.parquet": Meaning(
        "As OCF, for a whole-genome run -- no capture, so no on/off-target " "split.",
        "as OCF",
        pon_columns=("ocf_z",),
    ),
    ".FSC.regions.e1only.parquet": Meaning(
        "As FSC by region, restricted to first exons.",
        "as FSC",
    ),
    ".fsc_counts.parquet": Meaning(
        "The raw per-bin fragment counts behind FSC, before GC correction and "
        "before conversion to ratios.",
        None,
        why=(
            "FSC ships ratios, which are scale-free and therefore hide how much evidence is behind them. This is the count each ratio was computed from."
        ),
        how=(
            "The same single pass that produces FSC, written out before normalisation."
        ),
        what=(
            "`total` per bin. A ratio from a bin with a handful of fragments is noise with a plausible-looking value, and the ratio alone cannot tell you which bins those are."
        ),
    ),
    ".correction_factors.parquet": Meaning(
        "The GC bias model this run applied: the correction factor fitted for "
        "each fragment-length and GC-content bin.",
        None,
        "Diagnostic output. Factors far from 1 mean the correction is doing a "
        "lot of work, which is a fact about the library, not about the tumour.",
        why=(
            "GC content biases coverage during library preparation and sequencing, and short fragments are affected differently from long ones. Every FSC and FSR ratio is corrected using this model, so the model is part of the result."
        ),
        how=(
            "Observed and expected fragment counts are tabulated per length/GC bin and the correction factor is their ratio; in panel mode the model is fitted on off-target fragments, which are unbiased by capture."
        ),
        what=(
            "How far the factors sit from 1, and whether any bin is fitted from few fragments -- `observed` says. A large correction derived from a thin bin is the case to distrust."
        ),
    ),
}
