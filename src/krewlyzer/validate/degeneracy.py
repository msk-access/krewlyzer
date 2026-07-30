"""Detect columns that carry no information.

A metric that cannot vary with the data is worse than a missing metric: it is
present, plausible, and survives every schema check, so it gets modelled. This
module is the reason the gate exists.

Two axes:

* **cross-sample** -- the column is bit-identical in every sample. On a real
  cohort this means the metric is not a function of the input at all.
* **within-sample** -- the column is constant down the rows of every sample,
  for tables where rows are independent observations.

Below ``min_samples`` a cross-sample verdict is reported as SKIP, never PASS.
One sample cannot demonstrate variation, and reporting "pass" there is exactly
how a gate ends up certifying a constant.
"""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .contract import ColumnRule, Kind, Vary
from .findings import Category, Finding, Severity

# Cap the rows used to fingerprint a column. Deterministic head slice: enough
# to distinguish samples in practice, bounded for the 200-float WPS vectors
# across ~15k anchors.
_SIGNATURE_ROWS = 2000


def _reduce_lists(head: pd.Series) -> pd.Series:
    """Collapse array-valued rows to one float each.

    Array columns cannot be hashed or counted directly -- ``Series.nunique``
    raises on the numpy arrays backing the WPS vectors -- and hashing millions
    of floats to compare two samples is wasteful. The mean moves whenever any
    position in the vector moves, which is all either caller needs.
    """
    return pd.Series(
        [float(np.mean(v)) if v is not None and len(v) else np.nan for v in head],
        dtype=float,
    )


def column_signature(series: pd.Series, kind: Kind) -> str:
    """A stable fingerprint of a column's contents."""
    head = series.iloc[:_SIGNATURE_ROWS]
    if kind is Kind.STRING:
        return hashlib.sha1("\x00".join(head.astype(str)).encode()).hexdigest()[:16]
    reduced = (
        _reduce_lists(head)
        if kind is Kind.LIST
        else pd.to_numeric(head, errors="coerce")
    )
    # NaN has no stable byte pattern to hash, so map it to a sentinel that no
    # real measurement produces.
    rounded = np.round(np.nan_to_num(reduced.to_numpy(dtype=float), nan=-9.87e30), 9)
    return hashlib.sha1(rounded.tobytes()).hexdigest()[:16]


def _nunique(series: pd.Series, kind: Kind) -> int:
    """Distinct-value count that also works for array-valued columns."""
    head = series.iloc[:_SIGNATURE_ROWS]
    if kind is Kind.LIST:
        return int(_reduce_lists(head).nunique(dropna=False))
    return int(series.nunique(dropna=False))


def _constant_value(series: pd.Series) -> Optional[float]:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().any() and numeric.nunique(dropna=True) == 1:
        return float(numeric.dropna().iloc[0])
    return None


@dataclass(frozen=True)
class Observation:
    """What one sample contributes about one column.

    Deliberately not the data. Retaining a Series per column per sample is fine
    for a handful of samples and an out-of-memory error at cohort scale, so the
    reduction happens at read time and only this survives.
    """

    sample: str
    signature: str
    n_distinct: int
    n_rows: int
    constant_value: Optional[float]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, raw: Dict[str, Any]) -> "Observation":
        return cls(
            sample=raw["sample"],
            signature=raw["signature"],
            n_distinct=int(raw["n_distinct"]),
            n_rows=int(raw["n_rows"]),
            constant_value=raw.get("constant_value"),
        )


def observe(sample: str, series: pd.Series, kind: Kind) -> Observation:
    return Observation(
        sample=sample,
        signature=column_signature(series, kind),
        n_distinct=_nunique(series, kind),
        n_rows=len(series),
        constant_value=_constant_value(series),
    )


def evaluate(
    suffix: str,
    rule: ColumnRule,
    per_sample: Sequence[Observation],
    min_samples: int,
) -> List[Finding]:
    """Judge one column across the samples that contained it."""
    if rule.vary is Vary.NEVER or not per_sample:
        return []

    findings: List[Finding] = []
    family = suffix.lstrip(".").replace(".parquet", "")
    wants_cross = rule.vary in (Vary.CROSS, Vary.BOTH)
    wants_within = rule.vary in (Vary.WITHIN, Vary.BOTH)

    # Two samples is the floor for the comparison to mean anything: with one,
    # every signature is trivially identical and the check would report the
    # whole cohort as degenerate. Clamp rather than trust the caller's flag.
    effective_min = max(min_samples, 2)

    if wants_cross:
        if len(per_sample) < effective_min:
            findings.append(
                Finding(
                    id=f"{family}.DEGENERACY.{rule.name}",
                    severity=Severity.SKIP,
                    category=Category.DEGENERACY,
                    message=(
                        f"cross-sample variation not evaluated: {len(per_sample)} "
                        f"sample(s) present, {effective_min} required"
                    ),
                    table=suffix,
                    column=rule.name,
                    samples=[o.sample for o in per_sample],
                    evidence={"n_samples": len(per_sample), "min_samples": min_samples},
                )
            )
        else:
            signatures = {o.sample: o.signature for o in per_sample}
            if len(set(signatures.values())) == 1:
                const = per_sample[0].constant_value
                shown = "identical" if const is None else f"identically {const:g}"
                findings.append(
                    Finding(
                        id=f"{family}.DEGENERACY.{rule.name}",
                        severity=Severity.ERROR,
                        category=Category.DEGENERACY,
                        message=(
                            f"'{rule.name}' is {shown} in all "
                            f"{len(per_sample)} samples -- it carries no "
                            "information and cannot be a function of the input"
                        ),
                        table=suffix,
                        column=rule.name,
                        samples=[o.sample for o in per_sample],
                        evidence={
                            "n_samples": len(per_sample),
                            "constant_value": const,
                            "signature": next(iter(signatures.values())),
                        },
                    )
                )

    if wants_within:
        flat = [o.sample for o in per_sample if o.n_distinct <= 1]
        if len(flat) == len(per_sample) and any(o.n_rows > 1 for o in per_sample):
            findings.append(
                Finding(
                    id=f"{family}.FLAT.{rule.name}",
                    severity=Severity.ERROR,
                    category=Category.DEGENERACY,
                    message=(
                        f"'{rule.name}' is constant down the rows of every "
                        f"sample ({len(flat)}/{len(per_sample)})"
                    ),
                    table=suffix,
                    column=rule.name,
                    samples=flat,
                    evidence={"n_samples": len(per_sample)},
                )
            )

    return findings
