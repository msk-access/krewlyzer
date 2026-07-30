"""The output contract krewlyzer owes its consumers.

This is a *declaration*, not a reimplementation of the writers: it records what
a finished output directory must look like so that a downstream reader can rely
on it. The immediate consumer is ``kreview``, which reads **parquet only** and
resolves each sample as ``{results_dir}/{sample_id}/{sample_id}{suffix}``.

Two properties of that consumer shape this file:

1. ``.metadata.parquet`` is a completion marker. A sample without it is dropped
   from the cohort silently -- not warned about, not errored on.
2. Every reader entry point swallows exceptions and yields an empty feature
   dict, so a malformed table degrades to missing columns rather than a crash.

Both mean a defect travels all the way to a model fit without anyone seeing it.
The contract therefore has to be asserted here, upstream, or not at all.

A column that never varies is a failure, not a pass. Four of the five
``WPS_background`` metrics were structurally constant across an entire
production cohort while passing every schema check that existed, which is the
reason ``vary`` is a required part of a column's declaration and
``Vary.NEVER`` costs you a written justification.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Tuple


class Vary(str, Enum):
    """How a column must change to be considered informative."""

    CROSS = "cross"  # must differ between samples
    WITHIN = "within"  # must differ between rows of one sample
    BOTH = "both"
    NEVER = "never"  # legitimately constant -- requires constant_reason


class Kind(str, Enum):
    """Loose type expectation. Deliberately coarse: pandas and pyarrow disagree
    about int64 vs float64 often enough that pinning exact dtypes produces
    false failures without catching real ones."""

    NUMERIC = "numeric"
    STRING = "string"
    LIST = "list"


@dataclass(frozen=True)
class ColumnRule:
    name: str
    kind: Kind = Kind.NUMERIC
    vary: Vary = Vary.CROSS
    required: bool = True
    constant_reason: Optional[str] = None

    def __post_init__(self) -> None:
        # The whole point of the gate. Silencing a degeneracy finding must cost
        # you a sentence explaining why the constant is correct, so that the
        # next person can judge whether it still is.
        if self.vary is Vary.NEVER and not self.constant_reason:
            raise ValueError(
                f"{self.name}: vary=NEVER requires constant_reason explaining "
                "why a constant value is correct for this column"
            )


@dataclass(frozen=True)
class Rows:
    """Row-count expectation."""

    exactly: Optional[int] = None
    at_least: Optional[int] = None

    def check(self, n: int) -> Optional[str]:
        if self.exactly is not None and n != self.exactly:
            return f"expected exactly {self.exactly} row(s), found {n}"
        if self.at_least is not None and n < self.at_least:
            return f"expected at least {self.at_least} row(s), found {n}"
        return None


@dataclass(frozen=True)
class TableRule:
    suffix: str
    columns: Tuple[ColumnRule, ...]
    rows: Rows = field(default_factory=lambda: Rows(at_least=1))
    checks: Tuple[str, ...] = ()  # names resolved in checks.py
    scan_rows: Optional[int] = None
    """Rows the domain checks need. ``None`` means all of them.

    Set it wherever reading the whole table is disproportionate: WPS carries
    200-float vectors over ~15k anchors and is ~120 MB per sample, but its only
    check inspects the first row, and fingerprinting for degeneracy needs a
    bounded slice rather than the lot. At cohort scale (>14k samples) the
    difference is hours of I/O.
    """

    @property
    def family(self) -> str:
        return self.suffix.lstrip(".")


# --------------------------------------------------------------------------
# Shorthand constructors -- the contract below is long enough that spelling
# out ColumnRule(...) for every entry buries the interesting fields.
# --------------------------------------------------------------------------


def metric(name: str, vary: Vary = Vary.CROSS) -> ColumnRule:
    return ColumnRule(name, Kind.NUMERIC, vary)


def label(name: str, reason: str) -> ColumnRule:
    """A key/identifier column. Constant-ness is expected, so it needs a why."""
    return ColumnRule(name, Kind.STRING, Vary.NEVER, constant_reason=reason)


_ID = "identifier column; its values are the join key, not a measurement"
_TISSUE = "fixed tissue/label vocabulary from the bundled atlas"

_ENTROPY_COLS = (
    label("label", _TISSUE),
    metric("count"),
    metric("mean_size"),
    metric("entropy"),
    metric("z_score"),
)

_MOTIF_COLS = (
    label("Motif", "the 256 4-mers are a fixed alphabet"),
    metric("Frequency", Vary.BOTH),
)

_FSC_RATIOS = (
    metric("ultra_short_ratio", Vary.BOTH),
    metric("core_short_ratio", Vary.BOTH),
    metric("mono_nucl_ratio", Vary.BOTH),
    metric("di_nucl_ratio", Vary.BOTH),
    metric("long_ratio", Vary.BOTH),
)

_FSR_COLS = (
    label("region", _ID),
    metric("total_count", Vary.BOTH),
    metric("short_long_ratio", Vary.BOTH),
    metric("short_long_log2", Vary.BOTH),
    metric("short_frac", Vary.BOTH),
    metric("long_frac", Vary.BOTH),
)

# All four are per-position vectors over the anchor window, not scalars --
# consumers reduce them with list_avg/list_max.
_WPS_COLS = (
    label("region_type", "fixed anchor vocabulary (TSS, CTCF)"),
    ColumnRule("wps_nuc", Kind.LIST, Vary.BOTH),
    ColumnRule("wps_tf", Kind.LIST, Vary.BOTH),
    ColumnRule("prot_frac_nuc", Kind.LIST, Vary.BOTH),
    ColumnRule("prot_frac_tf", Kind.LIST, Vary.BOTH),
)


CONTRACT: Tuple[TableRule, ...] = (
    # -- fragment size ------------------------------------------------------
    TableRule(
        ".FSD.parquet",
        (label("region", _ID), metric("total", Vary.BOTH)),
        checks=("fsd_only_size_bins",),
    ),
    TableRule(
        ".FSD.ontarget.parquet",
        (label("region", _ID), metric("total", Vary.BOTH)),
        checks=("fsd_only_size_bins",),
    ),
    TableRule(".FSR.parquet", _FSR_COLS, checks=("fsr_region_format",)),
    TableRule(".FSR.ontarget.parquet", _FSR_COLS),
    TableRule(
        ".FSC.parquet",
        (label("chrom", _ID), metric("total", Vary.BOTH)),
        checks=("chr_prefixed", "fsc_has_log2", "fsc_channels_sum_to_total"),
    ),
    TableRule(
        ".FSC.ontarget.parquet",
        (label("chrom", _ID), metric("total", Vary.BOTH)),
        checks=("chr_prefixed", "fsc_has_log2"),
    ),
    TableRule(
        ".FSC.gene.parquet",
        (label("gene", _ID), *_FSC_RATIOS, metric("normalized_depth", Vary.BOTH)),
        checks=("fsc_gene_ratios_sum_to_one",),
    ),
    TableRule(
        ".FSC.regions.parquet",
        (label("gene", _ID), metric("total", Vary.BOTH), *_FSC_RATIOS),
    ),
    # -- motif --------------------------------------------------------------
    TableRule(
        ".EndMotif.parquet",
        _MOTIF_COLS,
        rows=Rows(exactly=256),
        checks=("acgt_4mers", "frequency_sums_to_one"),
    ),
    TableRule(
        ".EndMotif.ontarget.parquet",
        _MOTIF_COLS,
        rows=Rows(exactly=256),
        checks=("acgt_4mers", "frequency_sums_to_one"),
    ),
    TableRule(
        ".BreakPointMotif.parquet",
        _MOTIF_COLS,
        rows=Rows(exactly=256),
        checks=("acgt_4mers",),
    ),
    TableRule(
        ".BreakPointMotif.ontarget.parquet",
        _MOTIF_COLS,
        rows=Rows(exactly=256),
        checks=("acgt_4mers",),
    ),
    TableRule(
        ".EndMotif1mer.parquet",
        (label("base", "the four DNA bases"), metric("fraction", Vary.BOTH)),
        rows=Rows(exactly=4),
        checks=("acgt_bases", "fraction_sums_to_one"),
    ),
    TableRule(".MDS.parquet", (metric("MDS"),), rows=Rows(exactly=1)),
    TableRule(".MDS.ontarget.parquet", (metric("MDS"),), rows=Rows(exactly=1)),
    TableRule(
        ".MDS.gene.parquet",
        (
            label("gene", _ID),
            metric("mds_mean", Vary.BOTH),
            metric("mds_e1", Vary.BOTH),
            metric("mds_std", Vary.BOTH),
            metric("mds_z", Vary.BOTH),
            metric("mds_e1_z", Vary.BOTH),
        ),
    ),
    TableRule(".MDS.exon.parquet", (label("gene", _ID), metric("mds", Vary.BOTH))),
    # -- regulatory ---------------------------------------------------------
    TableRule(
        ".OCF.ontarget.parquet",
        (
            label("tissue", _TISSUE),
            metric("OCF", Vary.BOTH),
            metric("ocf_z", Vary.BOTH),
        ),
    ),
    TableRule(
        ".OCF.offtarget.parquet",
        (
            label("tissue", _TISSUE),
            metric("OCF", Vary.BOTH),
            metric("ocf_z", Vary.BOTH),
        ),
    ),
    TableRule(".TFBS.parquet", _ENTROPY_COLS),
    TableRule(".TFBS.ontarget.parquet", _ENTROPY_COLS),
    TableRule(".ATAC.parquet", _ENTROPY_COLS),
    TableRule(".ATAC.ontarget.parquet", _ENTROPY_COLS),
    # -- nucleosome ---------------------------------------------------------
    TableRule(
        ".WPS.parquet", _WPS_COLS, checks=("wps_arrays_nonempty",), scan_rows=256
    ),
    TableRule(
        ".WPS.panel.parquet",
        (*_WPS_COLS, metric("local_depth", Vary.BOTH)),
        checks=("wps_arrays_nonempty",),
        scan_rows=256,
    ),
    TableRule(
        ".WPS_background.parquet",
        (
            label("group_id", _ID),
            metric("nrl_bp"),
            metric("nrl_deviation_bp"),
            metric("periodicity_score"),
            metric("adjusted_score"),
            metric("fragment_ratio", Vary.BOTH),
        ),
        checks=("unique_group_id",),
    ),
    # -- completion marker --------------------------------------------------
    TableRule(
        ".metadata.parquet",
        (ColumnRule("sample_id", Kind.STRING, Vary.CROSS), metric("total_fragments")),
        rows=Rows(exactly=1),
    ),
)

COMPLETION_MARKER = ".metadata.parquet"

# Produced by krewlyzer but not read by kreview. Inventoried so that a missing
# one is visible, never gated.
NOT_CONSUMED: Tuple[str, ...] = (
    ".UXM.parquet",
    ".mFSD.parquet",
    ".FSC.regions.e1only.parquet",
    ".fsc_counts.parquet",
    ".correction_factors.parquet",
    ".OCF.parquet",
)
