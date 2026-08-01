"""Rendering for contract-gate results: a table to read, JSON to track."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

from rich.console import Console
from rich.table import Table

from .findings import Category, Finding, Severity
from .gate import Result

_STYLE = {
    Severity.ERROR: "red",
    Severity.WARN: "yellow",
    Severity.SKIP: "dim",
}

_ORDER = {Severity.ERROR: 0, Severity.WARN: 1, Severity.SKIP: 2}


def _sorted(findings: List[Finding]) -> List[Finding]:
    return sorted(findings, key=lambda f: (_ORDER[f.severity], f.table or "", f.id))


#: Shown when a run has degeneracy errors, unless the caller supplies its own.
#:
#: The two gates have different remedies -- an output column can be declared
#: `vary=NEVER` with a written reason, a PON block cannot -- and offering the
#: wrong one is worse than offering none.
DEFAULT_DEGENERACY_NOTE = (
    "a constant column passes every schema check and still reaches a model "
    "fit. Fix the metric, or declare it vary=NEVER with a written reason -- "
    "do not silence the check."
)

PON_DEGENERACY_NOTE = (
    "a baseline that cannot vary with its cohort was not fitted to one. There "
    "is no way to declare this expected -- rebuild the PON with build-pon, and "
    "check the builder is reading the columns it thinks it is."
)


def render(
    result: Result,
    console: Console,
    degeneracy_note: str = DEFAULT_DEGENERACY_NOTE,
) -> None:
    counts = result.counts()
    console.print(
        f"\n[bold]{len(result.samples)} sample(s)[/bold] · "
        f"[red]{counts['error']} error[/red] · "
        f"[yellow]{counts['warn']} warning[/yellow] · "
        f"[dim]{counts['skip']} skipped[/dim]"
    )

    if not result.findings:
        console.print("[green]contract satisfied[/green]")
        return

    table = Table(show_lines=False, header_style="bold")
    table.add_column("severity", width=8)
    table.add_column("category", width=11)
    table.add_column("table", width=30, overflow="fold")
    table.add_column("finding", overflow="fold")

    for f in _sorted(result.findings):
        n = len(f.samples)
        scope = f" [dim]({n} sample{'s' if n != 1 else ''})[/dim]" if n else ""
        table.add_row(
            f"[{_STYLE[f.severity]}]{f.severity.value}[/{_STYLE[f.severity]}]",
            f.category.value,
            (f.table or "").lstrip("."),
            f.message + scope,
        )
    console.print(table)

    degenerate = [
        f
        for f in result.findings
        if f.category is Category.DEGENERACY and f.severity is Severity.ERROR
    ]
    if degenerate:
        console.print(f"\n[bold]Note:[/bold] {degeneracy_note}")


def to_json(result: Result, path: Path, version: str) -> None:
    payload: Dict[str, Any] = {
        "schema_version": "1",
        "krewlyzer_version": version,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_samples": len(result.samples),
        "samples": result.samples,
        "summary": result.counts(),
        "exit_code": result.exit_code,
        "findings": [f.to_dict() for f in _sorted(result.findings)],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
