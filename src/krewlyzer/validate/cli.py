"""``krewlyzer validate-output`` -- check a finished output directory."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import typer
from rich.console import Console

from . import report
from .findings import Category, Finding, Severity
from .gate import Fingerprint, Result, evaluate_cohort, run

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
    fingerprint_out: Optional[Path] = typer.Option(
        None,
        "--fingerprint-out",
        help="Write the per-sample fingerprint(s) for a later validate-cohort "
        "pass. Kilobytes per sample; this is what lets a cohort be checked "
        "without re-reading it.",
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

    if fingerprint_out is not None:
        if len(result.fingerprints) == 1:
            result.fingerprints[0].save(fingerprint_out)
        else:
            fingerprint_out.mkdir(parents=True, exist_ok=True)
            for fp in result.fingerprints:
                fp.save(fingerprint_out / f"{fp.sample}.fingerprint.json")
        console.print(f"[dim]wrote fingerprint(s) to {fingerprint_out}[/dim]")

    raise typer.Exit(result.exit_code)


def validate_cohort(
    fingerprints: List[Path] = typer.Argument(
        ...,
        exists=True,
        help="Fingerprint JSON files (or directories of them) written by "
        "validate-output --fingerprint-out.",
    ),
    min_samples: int = typer.Option(3, "--min-samples"),
    json_report: Optional[Path] = typer.Option(None, "--json-report"),
) -> None:
    """Cross-sample half of the gate: find metrics that never vary.

    Degeneracy cannot be judged one sample at a time -- a single sample cannot
    distinguish "this metric is a constant" from "this is its value here" --
    so it is deliberately split out. Scatter validate-output across the cohort,
    gather the fingerprints, run this. It reads kilobytes per sample instead of
    re-reading gigabytes.
    """
    from krewlyzer import __version__

    paths: List[Path] = []
    for entry in fingerprints:
        paths.extend(sorted(entry.glob("*.json")) if entry.is_dir() else [entry])

    loaded = []
    result = Result()
    for path in paths:
        try:
            loaded.append(Fingerprint.load(path))
        except Exception as exc:
            result.findings.append(
                Finding(
                    id="FINGERPRINT.UNREADABLE",
                    severity=Severity.ERROR,
                    category=Category.STRUCTURAL,
                    message=f"could not read {path.name}: {exc}",
                )
            )

    result.fingerprints = loaded
    result.samples = [fp.sample for fp in loaded]
    if not loaded and not result.findings:
        result.findings.append(
            Finding(
                id="FINGERPRINT.NONE",
                severity=Severity.ERROR,
                category=Category.STRUCTURAL,
                message="no fingerprints found",
            )
        )
    result.findings.extend(evaluate_cohort(loaded, min_samples=min_samples))

    report.render(result, console)
    if json_report is not None:
        report.to_json(result, json_report, __version__)
        console.print(f"[dim]wrote {json_report}[/dim]")

    raise typer.Exit(result.exit_code)
