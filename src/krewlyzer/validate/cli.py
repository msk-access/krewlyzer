"""``krewlyzer validate-output`` -- check a finished output directory."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import typer
from rich.console import Console

from . import report
from .gate import run

console = Console(stderr=True)


def validate_output(
    results_dir: Path = typer.Argument(
        ...,
        exists=True,
        help="Cohort directory laid out as {results_dir}/{sample_id}/, or a "
        "single sample directory.",
    ),
    sample_id: Optional[List[str]] = typer.Option(
        None, "--sample-id", "-s", help="Restrict to these samples (repeatable)."
    ),
    min_samples: int = typer.Option(
        3,
        "--min-samples",
        help="Below this, cross-sample degeneracy is reported SKIP, not PASS: "
        "too few samples to demonstrate that a metric varies.",
    ),
    json_report: Optional[Path] = typer.Option(
        None, "--json-report", help="Write findings as JSON for trend tracking."
    ),
) -> None:
    """Validate outputs against the downstream contract.

    Checks three things, in increasing order of what a schema alone can catch:
    that every consumed table is present and shaped correctly, that domain
    invariants hold (frequencies sum to 1, chromosomes are chr-prefixed,
    channels partition the total), and that no metric is constant across the
    cohort.

    Exit codes: 0 satisfied, 1 contract violation, 2 structural -- the inputs
    are not comparable (missing directory, unreadable Parquet), so a caller can
    retry on 2 and escalate on 1.
    """
    from krewlyzer import __version__

    result = run(results_dir, min_samples=min_samples, only_samples=sample_id)
    report.render(result, console)

    if json_report is not None:
        report.to_json(result, json_report, __version__)
        console.print(f"[dim]wrote {json_report}[/dim]")

    raise typer.Exit(result.exit_code)
