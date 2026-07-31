"""Describe what a finished sample directory actually contains.

`validate-output` answers *"is this correct?"*. This answers the question people
ask first: *"what are all these files, and what is in them?"*

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
from krewlyzer.validate.meaning import MEANINGS, Meaning
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
    #: One-line interpretation from `meaning.py`. Absent only if a table were
    #: added to the contract without one, which a test forbids.
    meaning: Optional[Meaning] = None

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
    raise AssertionError("unreachable: the GB arm above always returns")


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
    # Ranges are withheld for identifier columns too. The rendering below
    # prefers min/max over `example`, so redacting only the example left a
    # numeric identifier -- an accession, an integer patient key -- printed in
    # full by the very branch that was supposed to hide it. A range over one
    # distinct value *is* the value.
    if pd.api.types.is_numeric_dtype(series) and not is_identifier_column(name):
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
        meaning=MEANINGS.get(rule.suffix),
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
    """Inspect every output krewlyzer can produce for one sample.

    ``sample_id`` defaults to the directory name, which is the layout the
    downstream consumer assumes: ``{results_dir}/{sample_id}/{sample_id}{suffix}``.

    Covers ``CONTRACT`` *and* ``NOT_CONSUMED``. The latter is optional output —
    mFSD needs a variant list, UXM a methylation-aware BAM — and gating on it
    would be wrong, but a report that answers "what is in this folder" has to
    describe a file that is sitting in the folder. Omitting them made the
    report say mFSD was absent while the file was right there.
    """
    directory = Path(directory)
    sample_id = sample_id or directory.name
    not_consumed = set(NOT_CONSUMED)

    tables = [
        describe_table(directory, sample_id, rule, rule.suffix not in not_consumed)
        for rule in CONTRACT
    ]

    # Optional outputs have no TableRule -- nothing declares their columns,
    # because nothing downstream requires them. An empty rule is honest here:
    # every column is still measured, none is checked against a declaration.
    described = {t.suffix for t in tables}
    for suffix in NOT_CONSUMED:
        if suffix in described:
            continue
        tables.append(
            describe_table(
                directory, sample_id, TableRule(suffix, columns=()), consumed=False
            )
        )

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
    gated = sum(1 for t in present if t.consumed)
    add(
        f"**{len(present)} of {len(report.tables)} possible outputs present** "
        f"({gated} gated), {human_bytes(report.total_bytes)} total.\n"
    )
    if not report.has_completion_marker:
        add(
            f"> **This sample has no `{COMPLETION_MARKER}`.** The downstream "
            "consumer treats that file as its completion marker and drops the "
            "sample from the cohort silently — no warning, no error.\n"
        )

    add("## Which way is the tumour signal?\n")
    add(
        "Direction differs per axis, and getting one backwards is the commonest "
        "misreading of this output — `MDS` was documented the wrong way round "
        "for a year. A single elevated metric is a hypothesis; several agreeing "
        "across independent axes is evidence.\n"
    )
    add("| Table | Measures | Cancer direction |")
    add("|---|---|---|")
    for t in present:
        m = t.meaning
        if not m or not m.cancer_direction:
            continue
        add(f"| `{t.family}` | {m.measures.split('.')[0]}. | {m.cancer_direction} |")
    add("")
    add(
        "*No thresholds are given. Every numeric band examined turned out to be "
        "a display default or refuted outright — the documented ATAC/TFBS "
        "entropy range flags a perfectly healthy distribution as abnormal. "
        "Directions are robust; magnitudes are cohort-specific.*\n"
    )

    add("## Contents\n")
    add(
        "*Gated* tables are checked by `krewlyzer validate-output`; a failure "
        "there fails the run. *Inventory* tables are described and never gate — "
        "they need an input the pipeline cannot assume, so an absent one is not "
        "a fault. Neither says whether a downstream tool reads the table: that "
        "is a fact about another repository, and a label nothing keeps in sync "
        "is worse than no label.\n"
    )
    add("| Table | Rows | Cols | Size | Contract |")
    add("|---|---:|---:|---:|:---:|")
    for t in present:
        rows = f"{t.n_rows:,}" if t.n_rows is not None else "—"
        add(
            f"| [`{t.family}`](#{t.family.lower().replace('.', '')}) | {rows} | "
            f"{t.n_cols} | {human_bytes(t.size_bytes)} | "
            f"{'gated' if t.consumed else 'inventory'} |"
        )
    add("")

    if missing:
        add("### Absent\n")
        for t in missing:
            why = "" if t.consumed else " *(inventory — never gates a run)*"
            add(f"- `{t.family}`{why}")
        add("")

    add("## Tables\n")
    for t in present:
        add(f"### {t.family}\n")
        meta = [
            f"{t.n_rows:,} rows × {t.n_cols} columns",
            human_bytes(t.size_bytes),
            "gated" if t.consumed else "inventory",
        ]
        add(" · ".join(meta) + "\n")
        if t.meaning:
            add(t.meaning.measures + "\n")
            if t.meaning.cancer_direction:
                add(f"**Cancer direction:** {t.meaning.cancer_direction}\n")
            if t.meaning.caveat:
                add(f"> ⚠️ {t.meaning.caveat}\n")
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
        "Shape, ranges and examples are measured from this sample. The one-line "
        "interpretation and direction for each table come from the output "
        f"contract, so they cannot drift from it; the longer narrative — how to "
        f"use each table, worked examples — is in the [output reference]({_DOCS}).\n"
    )
    add(
        "Values in identifier columns are redacted. Everything else is this "
        "sample's real data.\n"
    )
    return "\n".join(out)


def render_html_page(report: SampleReport, markdown: str) -> str:
    """Wrap the Markdown in a self-contained page, attachable or hostable as-is.

    Named apart from ``htmlreport.render_html`` on purpose: that one builds
    the full charted report, this one only wraps ``render_markdown`` output.

    Markdown is embedded and rendered client-side only if `markdown` is
    installed; otherwise it is shown in a <pre>, which stays readable. A hard
    dependency on a renderer would make the whole command unavailable in a
    minimal install for the sake of formatting.
    """
    try:
        import markdown as _md  # type: ignore[import-untyped]

        body = _md.markdown(markdown, extensions=["tables", "toc"])
    except ImportError:
        from html import escape

        body = f"<pre>{escape(markdown)}</pre>"

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>krewlyzer output — {report.sample_id}</title>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 60rem; margin: 2rem auto; padding: 0 1rem; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; display: block;
          overflow-x: auto; }}
  th, td {{ border: 1px solid #8884; padding: .35rem .6rem; text-align: left;
           white-space: nowrap; }}
  th {{ background: #8881; }}
  code {{ font-size: .9em; }}
  blockquote {{ border-left: 3px solid #e5a; margin: 1rem 0; padding: .1rem 1rem;
               background: #e5a1; }}
  h3 {{ margin-top: 2rem; border-top: 1px solid #8883; padding-top: 1rem; }}
</style></head><body>
{body}
</body></html>"""
