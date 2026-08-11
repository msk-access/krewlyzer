"""Named domain checks referenced by ``TableRule.checks``.

Each check takes a DataFrame and returns a list of human-readable problems --
empty means it passed. The caller turns them into Findings, so checks stay
free of reporting concerns and are trivial to unit test.

Float tolerance policy: there is no global epsilon, and the tolerance is
derived rather than chosen. For every sum-style check the dominant error is not
float arithmetic but the **decimal precision the value was written at** --
krewlyzer's writers emit ``%.2f`` / ``%.4f`` / ``%.6f`` depending on the table.
A value written at *d* decimals carries up to ``0.5e-d`` of rounding, and
summing *n* of them accumulates *n* times that in the worst case. So each check
infers *d* from the data and scales by the number of terms.

Deriving it this way matters: assuming exact integer conservation for the FSC
channels is wrong, because GC correction makes them weighted floats, and a
hand-picked epsilon would either mask real drift on 6-decimal tables or
manufacture failures on 2-decimal ones.
"""

from __future__ import annotations

import re
from typing import Callable, Dict, List

import numpy as np
import pandas as pd

ACGT = set("ACGT")
_SIZE_BIN = re.compile(r"^\d+-\d+$")
#: The PON log-ratio for one size bin, e.g. `65-69_logR`. Named from the bin, so
#: it cannot be enumerated in the contract the way a fixed column can.
_SIZE_BIN_LOGRATIO = re.compile(r"^\d+-\d+_logR$")
_REGION = re.compile(r"^chr[\dXYMT]+:\d+-\d+$")


_MAX_DECIMALS = 8


def written_decimals(series: pd.Series, sample: int = 500) -> int:
    """Infer the decimal precision the column was written at.

    Capped, because a genuinely full-precision float (1/3) would otherwise
    report 17 and make the tolerance meaninglessly tight.
    """
    values = pd.to_numeric(series, errors="coerce").dropna().iloc[:sample]
    if values.empty:
        return _MAX_DECIMALS
    widest = 0
    for v in values:
        text = repr(float(v))
        if "e" in text or "E" in text:
            return _MAX_DECIMALS
        frac = text.partition(".")[2]
        widest = max(widest, len(frac.rstrip("0")))
        if widest >= _MAX_DECIMALS:
            return _MAX_DECIMALS
    return widest


def quantization_tolerance(n_terms: int, decimals: int) -> float:
    """Worst-case rounding accumulated by summing n values written at d dp."""
    return max(n_terms, 1) * 0.5 * (10.0**-decimals)


def _sum_tolerance(series: pd.Series, n_terms: int) -> float:
    dtype_floor = (
        n_terms * float(np.finfo(np.float32).eps) * 4
        if series.dtype == np.float32
        else 1e-9
    )
    return max(dtype_floor, quantization_tolerance(n_terms, written_decimals(series)))


def _sums_to_one(df: pd.DataFrame, column: str) -> List[str]:
    if column not in df.columns:
        return [f"missing column '{column}'"]
    total = float(df[column].sum())
    tol = _sum_tolerance(df[column], len(df))
    if abs(total - 1.0) > tol:
        return [f"'{column}' sums to {total:.9f}, expected 1.0 +/- {tol:.2e}"]
    return []


def frequency_sums_to_one(df: pd.DataFrame) -> List[str]:
    return _sums_to_one(df, "Frequency")


def fraction_sums_to_one(df: pd.DataFrame) -> List[str]:
    return _sums_to_one(df, "fraction")


def acgt_4mers(df: pd.DataFrame) -> List[str]:
    if "Motif" not in df.columns:
        return ["missing column 'Motif'"]
    motifs = df["Motif"].astype(str)
    bad_len = motifs[motifs.str.len() != 4]
    if not bad_len.empty:
        return [
            f"{len(bad_len)} motif(s) are not 4-mers, e.g. {bad_len.iloc[0]!r}. "
            "Downstream purine/pyrimidine and DNase1L3 features key on a fixed "
            "4-mer alphabet and silently vanish otherwise"
        ]
    non_acgt = motifs[~motifs.apply(lambda m: set(m) <= ACGT)]
    if not non_acgt.empty:
        return [
            f"{len(non_acgt)} motif(s) contain non-ACGT bases, e.g. "
            f"{non_acgt.iloc[0]!r}; consumers compute entropy over the whole "
            "column without filtering"
        ]
    return []


def acgt_bases(df: pd.DataFrame) -> List[str]:
    if "base" not in df.columns:
        return ["missing column 'base'"]
    found = set(df["base"].astype(str))
    if found != ACGT:
        return [f"expected bases {sorted(ACGT)}, found {sorted(found)}"]
    return []


def chr_prefixed(df: pd.DataFrame) -> List[str]:
    if "chrom" not in df.columns:
        return ["missing column 'chrom'"]
    values = df["chrom"].astype(str)
    bare = values[~values.str.startswith("chr")]
    if not bare.empty:
        return [
            f"{len(bare)} row(s) have a bare chromosome name, e.g. "
            f"{bare.iloc[0]!r}; per-chromosome features are keyed on 'chrN' "
            "and silently produce nothing for Ensembl-style names"
        ]
    return []


def fsr_region_format(df: pd.DataFrame) -> List[str]:
    if "region" not in df.columns:
        return ["missing column 'region'"]
    values = df["region"].astype(str)
    bad = values[~values.str.match(_REGION)]
    if not bad.empty:
        return [
            f"{len(bad)} region label(s) do not match 'chrN:start-end', e.g. "
            f"{bad.iloc[0]!r}; the chromosome is parsed back out with a regex"
        ]
    return []


def fsd_only_size_bins(df: pd.DataFrame) -> List[str]:
    """Every numeric column besides the reserved ones is treated as a size bin.

    A stray numeric column (a QC metric, a count) silently becomes a fake
    density feature downstream and corrupts the derived entropy, so the column
    set is part of the contract rather than an implementation detail.
    """
    # PON scoring writes into this same table, and the reserved set was never
    # updated when it started: `{bin}_logR` and `pon_stability` are the
    # documented product of `apply_pon_logratio`, and this check reported all
    # 68 of them as stray on every scored sample. A gate that fires on correct
    # output is a gate nobody reads. It went unnoticed because the synthetic
    # cohort carried no PON columns until 0.9.0, so the check had never once
    # seen the scored shape it runs against in production.
    #
    # Worth confirming against the consumer if the chance arises: the message
    # below asserts that *any* unrecognised numeric column is taken for a size
    # bin. If that is literally true downstream rather than a description of
    # the `\d+-\d+` convention, these columns would need their own table --
    # but they have shipped in this one for several releases.
    reserved = {"region", "total", "chrom", "pon_stability"}
    numeric = set(df.select_dtypes(include="number").columns)
    stray = sorted(
        c
        for c in numeric - reserved
        if not _SIZE_BIN.match(str(c)) and not _SIZE_BIN_LOGRATIO.match(str(c))
    )
    if stray:
        return [
            f"non-size-bin numeric column(s) {stray}: consumers treat every "
            "numeric column that is not region/total/chrom as a fragment-size "
            "bucket"
        ]
    return []


def fsc_has_log2(df: pd.DataFrame) -> List[str]:
    if not [c for c in df.columns if str(c).endswith("_log2")]:
        return ["no '*_log2' columns; consumers select them by suffix"]
    return []


def fsc_channels_sum_to_total(df: pd.DataFrame) -> List[str]:
    """The six channels must partition 'total'.

    Not an integer check: GC correction makes these weighted floats, written
    at the table's decimal precision, so six rounded addends can miss their
    rounded total by six half-ulps of that precision.
    """
    channels = [
        "ultra_short",
        "core_short",
        "mono_nucl",
        "di_nucl",
        "long",
        "ultra_long",
    ]
    missing = [c for c in channels if c not in df.columns]
    if missing or "total" not in df.columns:
        return []  # schema check already reports the absence
    tol = quantization_tolerance(len(channels) + 1, written_decimals(df["total"]))
    diff = (df[channels].sum(axis=1) - df["total"]).abs()
    off = int((diff > tol).sum())
    if off:
        worst = float(diff.max())
        return [
            f"{off} bin(s) where the six channels do not sum to 'total' "
            f"(largest gap {worst:g}, tolerance {tol:g}); channel/total ratios "
            "are only interpretable if the channels partition the total"
        ]
    return []


#: The five ratio columns kreview consumed before ``ultra_long_ratio`` existed.
#:
#: They sum to ``1 - ultra_long_ratio``, **not** to 1. A consumer that selects
#: only these five and renormalises is dividing by the wrong base. Kept as a
#: named list so that relation has somewhere to be written down, and so a
#: change to the column set is a change to this line rather than a silent
#: shift in what the five mean.
KREVIEW_FSC_RATIOS = [
    "ultra_short_ratio",
    "core_short_ratio",
    "mono_nucl_ratio",
    "di_nucl_ratio",
    "long_ratio",
]

FSC_GENE_RATIOS = KREVIEW_FSC_RATIOS + ["ultra_long_ratio"]


def fsc_gene_ratios_sum_to_one(df: pd.DataFrame) -> List[str]:
    """All six gene/region ratios must partition the fragments.

    ``ultra_long_ratio`` is the sixth channel, added when the gene-level bands
    were aligned to the genome-bin bands. Before that the gene path had five
    channels whose ``long`` silently absorbed everything over 400bp, so those
    five summed to 1 -- which is part of why nobody noticed the bands had
    drifted.

    Only one relation is asserted, because there is only one. The statement
    that kreview's original five sum to ``1 - ultra_long_ratio`` is the *same*
    equation rearranged::

        legacy5 - (1 - ultra_long) == (legacy5 + ultra_long) - 1

    so a separate check for it can never fail independently. It is documented
    on :data:`KREVIEW_FSC_RATIOS` and pinned end to end by
    ``tests/invariants/test_fsc_band_alignment.py`` against real output, where
    it guards against the column set changing rather than the arithmetic.
    """
    missing = [c for c in FSC_GENE_RATIOS if c not in df.columns]
    if missing:
        return []  # schema check already reports the absence
    if df.empty:
        return []

    tol = quantization_tolerance(
        len(FSC_GENE_RATIOS), written_decimals(df["ultra_short_ratio"])
    )
    summed = df[FSC_GENE_RATIOS].sum(axis=1)
    off = int(((summed - 1.0).abs() > tol).sum())
    if not off:
        return []

    worst = float((summed - 1.0).abs().max())
    return [
        f"{off} row(s) whose six ratios sum to 1 +/- more than {tol:g} "
        f"(worst deviation {worst:.6f}); the size channels no longer partition "
        "the counted fragments, so every channel/total ratio is against the "
        "wrong base"
    ]


def wps_arrays_nonempty(df: pd.DataFrame) -> List[str]:
    problems = []
    for col in ("wps_nuc", "wps_tf"):
        if col not in df.columns:
            continue
        first = df[col].iloc[0] if len(df) else None
        if first is None or len(first) == 0:
            problems.append(
                f"'{col}' holds empty arrays; consumers aggregate these with "
                "list_avg/list_max and hard-fail on a type mismatch"
            )
    return problems


#: Above this, a z-score is a broken divisor rather than biology.
#:
#: A z of 100 means a sample sits a hundred healthy standard deviations from
#: the cohort. Real cfDNA does not do that; a sigma that should have been
#: absent does. Deliberately loose -- the point is to catch arithmetic that has
#: gone wrong by orders of magnitude, not to police interesting outliers.
#: Measured for reference: real WPS z-scores against the 0.9.0 models reach
#: 6.7, and 34,572 anchors over 6.9M positions produced nothing above 10.
IMPLAUSIBLE_Z = 100.0


def plausible_z_scores(df: pd.DataFrame) -> List[str]:
    """No z-score should be astronomically large.

    This is a divisor check wearing a statistics hat. Every fabricated-sigma
    defect in 0.9.0 -- the hardcoded 167.0/5.0, six sigma floors, a standard
    normal substituted for a missing block -- produced *plausible* numbers and
    survived every schema check. The failure mode this catches is the opposite
    one: a sigma so small the z explodes, which no schema notices either
    because the column is present and finite.

    Covers scalar `*_z` columns and vector ones (`wps_nuc_z` is 200 elements
    per row), since a per-position z is exactly where a near-zero sigma hides.
    """
    problems = []
    for col in [c for c in df.columns if c.endswith("_z") or c == "z_score"]:
        values = df[col].dropna()
        if values.empty:
            continue
        first = values.iloc[0]
        if isinstance(first, (list, tuple, np.ndarray)):
            flat = np.concatenate(
                [np.asarray(v, dtype=float) for v in values if v is not None]
            )
        else:
            flat = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
        finite = flat[np.isfinite(flat)]
        if finite.size == 0:
            continue
        worst = float(np.abs(finite).max())
        if worst > IMPLAUSIBLE_Z:
            n = int((np.abs(finite) > IMPLAUSIBLE_Z).sum())
            problems.append(
                f"'{col}' reaches |z| = {worst:.3g} ({n} value(s) over "
                f"{IMPLAUSIBLE_Z:g}); a z that large is a near-zero sigma, not "
                "a measurement -- check the PON block it was scored against"
            )
    return problems


def no_collided_columns(df: pd.DataFrame) -> List[str]:
    """A `.1`-suffixed or repeated column name is a frame collision on disk.

    pandas renames a duplicate column to `name.1` when it reads one, and writes
    that name straight back out. So the suffix in a shipped file means some
    step appended a set of columns that were already there -- and the reader
    cannot tell which copy is current. `read_csv` returns the first, which is
    the oldest.

    Found in FSD: `apply_pon_logratio` matched its own `65-69_logR` output as
    size bin 65, so every re-run appended another generation. A real 67-bin
    sample went 69 columns to 137 to 273, carrying `_logR`, `_logR.1` and
    `_logR.2` with the correct answer somewhere among them.

    Deliberately generic. Nothing about this check is FSD-specific, and the
    next step to append instead of replace will be caught without anyone
    thinking to look for it.
    """
    problems = []
    names = list(df.columns)

    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:
        problems.append(
            f"duplicate column name(s) {repeated[:5]}; a reader cannot tell "
            "which copy is current, and pandas returns the first"
        )

    # `foo.1`, `foo.12` -- pandas' own de-duplication, written back to disk.
    #
    # The un-suffixed name has to be present too. Without that condition every
    # `foo.N` flags itself, because a name always contributes its own base to
    # the comparison set, and legitimate names like `chr1.2` get reported.
    # Requiring the pair keeps this specific to an actual collision.
    present = set(names)
    suffixed = sorted(
        n
        for n in names
        if "." in n and n.rsplit(".", 1)[1].isdigit() and n.rsplit(".", 1)[0] in present
    )
    if suffixed:
        problems.append(
            f"column(s) {suffixed[:5]} carry a pandas de-duplication suffix; "
            "some step appended columns that were already present rather than "
            "replacing them"
        )
    return problems


def unique_group_id(df: pd.DataFrame) -> List[str]:
    if "group_id" not in df.columns:
        return ["missing column 'group_id'"]
    if not df["group_id"].is_unique:
        dupes = int(len(df) - df["group_id"].nunique())
        return [
            f"{dupes} duplicate group_id value(s); consumers key on it without "
            "aggregating, so duplicates are silently dropped"
        ]
    return []


REGISTRY: Dict[str, Callable[[pd.DataFrame], List[str]]] = {
    "acgt_4mers": acgt_4mers,
    "acgt_bases": acgt_bases,
    "chr_prefixed": chr_prefixed,
    "fraction_sums_to_one": fraction_sums_to_one,
    "frequency_sums_to_one": frequency_sums_to_one,
    "fsc_channels_sum_to_total": fsc_channels_sum_to_total,
    "fsc_gene_ratios_sum_to_one": fsc_gene_ratios_sum_to_one,
    "fsc_has_log2": fsc_has_log2,
    "fsd_only_size_bins": fsd_only_size_bins,
    "fsr_region_format": fsr_region_format,
    "no_collided_columns": no_collided_columns,
    "plausible_z_scores": plausible_z_scores,
    "unique_group_id": unique_group_id,
    "wps_arrays_nonempty": wps_arrays_nonempty,
}
