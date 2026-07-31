"""Describe what a finished sample directory actually contains.

`validate-output` answers *"is this correct?"*. This answers the question people
ask first: *"what are these 26 files, and what is in them?"*

Everything reported here is either **measured from the file** — row and column
counts, dtypes, a value sample, how much is null — or **read from
``contract.py``**, which already declares what each table owes its consumers.
Nothing about the biology is restated: that lives in
``docs/reference/output-files.md`` and is linked, not copied. A second prose
description of every column is exactly the kind of duplicate this codebase
keeps having to delete.

The consequence worth knowing: a table that gains a column, or stops being
consumed, changes this report automatically. It cannot go stale on its own.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from krewlyzer.core.output_utils import read_table
from krewlyzer.validate.contract import (
    COMPLETION_MARKER,
    CONTRACT,
    NOT_CONSUMED,
    TableRule,
    Vary,
)

#: How many distinct values to show for a column before calling it continuous.
_CATEGORICAL_MAX = 8

#: Rows read per table. The report describes shape and content, not every row,
#: and WPS carries 200-float vectors over ~15k anchors at ~120 MB per sample.
_SCAN_ROWS = 5000

#: Columns whose *values* are identifiers, never measurements.
#:
#: Their example values are redacted. This is not cosmetic: sample directories
#: here are named for the patient, and several tables carry the sample id as a
#: column value -- so a report generated from renamed files still leaked a real
#: identifier through `Sample` and `sample_id` examples. A report intended to be
#: hosted must not carry one.
#:
#: Knowing a column *holds* an identifier is the useful fact; which identifier
#: is not, and is the one thing that must not leave the machine.
_IDENTIFIER_COLUMNS = frozenset(
    {"sample", "sample_id", "sample_name", "patient", "patient_id", "subject", "id"}
)

#: What replaces a redacted value.
REDACTED = "<identifier — redacted>"


def is_identifier_column(name: str) -> bool:
    """Match on the normalised name, so `Sample` and `sample_id` both hit."""
    return name.strip().lower().replace("-", "_") in _IDENTIFIER_COLUMNS


@dataclass
class ColumnFacts:
    """What a column looks like, measured rather than declared."""

    name: str
    dtype: str
    n_null: int
    n_distinct: int
    example: str
    #: Present only for numeric columns with at least one finite value.
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    #: From the contract, when the column is declared there.
    declared_vary: Optional[str] = None
    constant_reason: Optional[str] = None

    @property
    def looks_categorical(self) -> bool:
        return self.n_distinct <= _CATEGORICAL_MAX


@dataclass
class TableFacts:
    """One output file: where it is, how big, and what is in it."""

    suffix: str
    path: Optional[Path]
    consumed: bool
    is_completion_marker: bool
    n_rows: Optional[int] = None
    n_cols: Optional[int] = None
    size_bytes: Optional[int] = None
    truncated: bool = False
    columns: List[ColumnFacts] = field(default_factory=list)
    error: Optional[str] = None

    @property
    def family(self) -> str:
        return self.suffix.lstrip(".").rsplit(".parquet", 1)[0]

    @property
    def present(self) -> bool:
        return self.path is not None and self.error is None


@dataclass
class SampleReport:
    sample_id: str
    directory: Path
    tables: List[TableFacts]

    @property
    def present(self) -> List[TableFacts]:
        return [t for t in self.tables if t.present]

    @property
    def missing(self) -> List[TableFacts]:
        return [t for t in self.tables if not t.present]

    @property
    def has_completion_marker(self) -> bool:
        return any(t.is_completion_marker and t.present for t in self.tables)

    @property
    def total_bytes(self) -> int:
        return sum(t.size_bytes or 0 for t in self.tables)


def human_bytes(n: Optional[int]) -> str:
    if not n:
        return "—"
    size = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} GB"


def _example(series: pd.Series) -> str:
    """A short, readable sample of the values.

    Lists are summarised by length rather than printed: ``wps_nuc`` holds a
    200-float vector per row, and dumping one would bury the report.
    """
    values = series.dropna()
    if values.empty:
        return "—"
    first = values.iloc[0]
    if (
        isinstance(first, (list, tuple))
        or hasattr(first, "__len__")
        and not isinstance(first, (str, bytes))
    ):
        try:
            return f"[{len(first)} values] {first[0]:.4g} …"
        except (TypeError, IndexError, ValueError):
            return f"[{len(first)} values]"
    if isinstance(first, float):
        return "—" if math.isnan(first) else f"{first:.6g}"
    text = str(first)
    return text if len(text) <= 40 else text[:37] + "…"


def _column_facts(
    series: pd.Series, name: str, rule_by_name: Dict[str, Any]
) -> ColumnFacts:
    declared = rule_by_name.get(name)
    try:
        n_distinct = int(series.nunique(dropna=True))
    except TypeError:
        # Unhashable cells -- the list-valued WPS columns. Distinctness is not
        # a meaningful question for them, and guessing would be worse than
        # saying so.
        n_distinct = -1

    facts = ColumnFacts(
        name=name,
        dtype=str(series.dtype),
        n_null=int(series.isna().sum()),
        n_distinct=n_distinct,
        example=REDACTED if is_identifier_column(name) else _example(series),
        declared_vary=declared.vary.value if declared else None,
        constant_reason=declared.constant_reason if declared else None,
    )
    # A numeric range cannot encode an identifier, so it is left alone; only
    # the example value is withheld.
    if pd.api.types.is_numeric_dtype(series):
        finite = pd.to_numeric(series, errors="coerce").dropna()
        if not finite.empty:
            facts.minimum = float(finite.min())
            facts.maximum = float(finite.max())
    return facts


def describe_table(
    directory: Path, sample_id: str, rule: TableRule, consumed: bool
) -> TableFacts:
    base = directory / f"{sample_id}{rule.suffix}"
    facts = TableFacts(
        suffix=rule.suffix,
        path=None,
        consumed=consumed,
        is_completion_marker=rule.suffix == COMPLETION_MARKER,
    )

    df = read_table(base)
    if df is None:
        return facts

    facts.path = base
    try:
        facts.size_bytes = base.stat().st_size
    except OSError:
        pass

    facts.n_rows = len(df)
    facts.n_cols = len(df.columns)
    if len(df) > _SCAN_ROWS:
        facts.truncated = True
        df = df.head(_SCAN_ROWS)

    rule_by_name = {c.name: c for c in rule.columns}
    facts.columns = [_column_facts(df[c], str(c), rule_by_name) for c in df.columns]
    return facts


def describe_sample(directory: Path, sample_id: Optional[str] = None) -> SampleReport:
    """Inspect every contracted output for one sample.

    ``sample_id`` defaults to the directory name, which is the layout the
    downstream consumer assumes: ``{results_dir}/{sample_id}/{sample_id}{suffix}``.
    """
    directory = Path(directory)
    sample_id = sample_id or directory.name
    not_consumed = set(NOT_CONSUMED)

    tables = [
        describe_table(directory, sample_id, rule, rule.suffix not in not_consumed)
        for rule in CONTRACT
    ]
    return SampleReport(sample_id=sample_id, directory=directory, tables=tables)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

_DOCS = "https://msk-access.github.io/krewlyzer/reference/output-files/"


def _vary_note(col: ColumnFacts) -> str:
    if col.declared_vary is None:
        return ""
    if col.declared_vary == Vary.NEVER.value:
        return f"constant by design — {col.constant_reason}"
    return {
        Vary.CROSS.value: "must differ between samples",
        Vary.WITHIN.value: "must differ between rows",
        Vary.BOTH.value: "must differ between rows and between samples",
    }.get(col.declared_vary, col.declared_vary)


def render_markdown(report: SampleReport) -> str:
    """A report someone can read top to bottom and know what they have."""
    out: List[str] = []
    add = out.append

    add(f"# Output report — `{report.sample_id}`\n")
    add(f"`{report.directory}`\n")

    present, missing = report.present, report.missing
    add(
        f"**{len(present)} of {len(report.tables)} contracted tables present**, "
        f"{human_bytes(report.total_bytes)} total.\n"
    )
    if not report.has_completion_marker:
        add(
            f"> **This sample has no `{COMPLETION_MARKER}`.** The downstream "
            "consumer treats that file as its completion marker and drops the "
            "sample from the cohort silently — no warning, no error.\n"
        )

    add("## Contents\n")
    add("| Table | Rows | Cols | Size | Read downstream |")
    add("|---|---:|---:|---:|:---:|")
    for t in present:
        rows = f"{t.n_rows:,}" if t.n_rows is not None else "—"
        add(
            f"| [`{t.family}`](#{t.family.lower().replace('.', '')}) | {rows} | "
            f"{t.n_cols} | {human_bytes(t.size_bytes)} | "
            f"{'yes' if t.consumed else 'no'} |"
        )
    add("")

    if missing:
        add("### Absent\n")
        for t in missing:
            why = "" if t.consumed else " *(not read downstream)*"
            add(f"- `{t.family}`{why}")
        add("")

    add("## Tables\n")
    for t in present:
        add(f"### {t.family}\n")
        meta = [
            f"{t.n_rows:,} rows × {t.n_cols} columns",
            human_bytes(t.size_bytes),
            "read downstream" if t.consumed else "not read downstream",
        ]
        add(" · ".join(meta) + "\n")
        if t.truncated:
            add(
                f"> Column facts sampled from the first {_SCAN_ROWS:,} rows; "
                "counts above are for the whole table.\n"
            )

        add("| Column | Type | Range / example | Distinct | Null | Contract |")
        add("|---|---|---|---:|---:|---|")
        for c in t.columns:
            if c.minimum is not None and c.maximum is not None:
                value = f"{c.minimum:.6g} … {c.maximum:.6g}"
            else:
                value = f"`{c.example}`"
            distinct = "—" if c.n_distinct < 0 else f"{c.n_distinct:,}"
            add(
                f"| `{c.name}` | {c.dtype} | {value} | {distinct} | "
                f"{c.n_null:,} | {_vary_note(c)} |"
            )
        add("")

    add("---\n")
    add(
        f"Shape, ranges and examples are measured from this sample. What each "
        f"column *means* is in the [output reference]({_DOCS}); this report "
        "deliberately does not restate it, so the two cannot disagree.\n"
    )
    return "\n".join(out)
