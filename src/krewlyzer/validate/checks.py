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
    reserved = {"region", "total", "chrom"}
    numeric = set(df.select_dtypes(include="number").columns)
    stray = sorted(c for c in numeric - reserved if not _SIZE_BIN.match(str(c)))
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


def fsc_gene_ratios_sum_to_one(df: pd.DataFrame) -> List[str]:
    ratios = [
        "ultra_short_ratio",
        "core_short_ratio",
        "mono_nucl_ratio",
        "di_nucl_ratio",
        "long_ratio",
    ]
    present = [c for c in ratios if c in df.columns]
    if len(present) != len(ratios):
        return []
    extra = "ultra_long_ratio"
    if extra in df.columns:
        present = present + [extra]
    summed = df[present].sum(axis=1)
    tol = quantization_tolerance(len(present), written_decimals(df[present[0]]))
    off = int(((summed - 1.0).abs() > tol).sum())
    if off:
        worst = float((summed - 1.0).abs().max())
        return [
            f"{off} gene(s) whose ratios sum to 1 +/- more than {tol:g} "
            f"(worst deviation {worst:.6f})"
        ]
    return []


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
    "unique_group_id": unique_group_id,
    "wps_arrays_nonempty": wps_arrays_nonempty,
}
