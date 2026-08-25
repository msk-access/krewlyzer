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


def read_expected(path: Path) -> List[str]:
    """Sample ids from a Nextflow samplesheet or a plain one-per-line list.

    Both, because the two are used at different moments: the samplesheet is the
    thing that was actually run and cannot drift from it, while a plain list is
    what you have when the cohort was assembled some other way.

    A CSV *without* a ``sample`` column is an error rather than a fallback to
    "first column". Guessing there would silently reconcile against the wrong
    field and report every sample missing, which reads as a catastrophic run
    failure instead of a malformed input.
    """
    import csv

    text = path.read_text().splitlines()
    if not text:
        raise typer.BadParameter(f"{path} is empty")

    header = text[0]
    if "," in header:
        reader = csv.DictReader(text)
        if reader.fieldnames and "sample" in reader.fieldnames:
            ids = [
                row["sample"].strip()
                for row in reader
                if row.get("sample") and row["sample"].strip()
            ]
            if not ids:
                raise typer.BadParameter(f"{path} has a 'sample' column but no rows")
            return ids
        raise typer.BadParameter(
            f"{path} looks like a CSV but has no 'sample' column "
            f"(found: {', '.join(reader.fieldnames or []) or 'nothing'}). "
            "Pass a plain list of ids instead, one per line."
        )

    ids = [line.strip() for line in text if line.strip() and not line.startswith("#")]
    if not ids:
        raise typer.BadParameter(f"{path} contains no sample ids")
    return ids


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
    expect: Optional[Path] = typer.Option(
        None,
        "--expect",
        exists=True,
        help="The samples that were meant to run: a Nextflow samplesheet (uses "
        "its 'sample' column) or a plain list, one id per line. Without this, "
        "a sample that produced nothing is invisible -- there is no directory "
        "to find, so the cohort is simply smaller than intended.",
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

    expected = read_expected(expect) if expect is not None else None
    if expected is not None:
        console.print(f"[dim]expecting {len(expected)} sample(s) from {expect}[/dim]")

    result = run(
        results_dir,
        min_samples=min_samples,
        only_samples=sample_id,
        expected_samples=expected,
    )
    report.render(result, console)

    if json_report is not None:
        report.to_json(result, json_report, __version__)
        console.print(f"[dim]wrote {json_report}[/dim]")

    if fingerprint_out is not None:
        # One sample writes the file you named; a cohort writes a directory of
        # them, since a single path cannot hold many.
        if len(result.fingerprints) == 1:
            result.fingerprints[0].save(fingerprint_out)
        else:
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


def validate_pon(
    pon_model: List[Path] = typer.Argument(
        ...,
        exists=True,
        help="PON parquet file(s) to check.",
    ),
    json_report: Optional[Path] = typer.Option(
        None, "--json-report", help="Write findings as JSON."
    ),
) -> None:
    """Check a built PON before anything is scored against it.

    `validate-output` gates the results of a run. Nothing gated the reference
    a run measures itself against -- and every PON defect this release fixes
    was sitting in four shipped files, visible in the file and invisible to
    every check.

    The load-bearing assertion is the same invariant as the output gate, one
    level up: a baseline that cannot vary with its cohort was not fitted to
    one. `wps_background` shipped a hardcoded 167.0/5.0 identical across all
    28 groups and all four models; a single check that sigma differs between
    groups would have caught it.

    Exit codes: 0 satisfied, 1 violation, 2 structural -- as validate-output,
    so a workflow can retry on 2 and escalate on 1.
    """
    from krewlyzer import __version__

    from .pon_gate import check_pon, describe, exit_code

    result = Result()
    for path in pon_model:
        provenance = describe(path)
        console.print(f"[bold]{path.name}[/bold]  {provenance or '(no metadata)'}")
        result.findings.extend(check_pon(path))

    report.render(
        result,
        console,
        degeneracy_note=report.PON_DEGENERACY_NOTE,
        unit="model",
        n_checked=len(pon_model),
    )
    if json_report is not None:
        report.to_json(result, json_report, __version__)
        console.print(f"[dim]wrote {json_report}[/dim]")

    raise typer.Exit(exit_code(result.findings))


def stamp_pon(
    pon_model: List[Path] = typer.Argument(
        ..., exists=True, help="PON parquet file(s) to stamp."
    ),
    version: str = typer.Option(
        ...,
        "--version",
        help="The release these models ship with, e.g. 0.9.0.",
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Report what would change and write nothing."
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Stamp even if validate-pon fails. For re-stamping a model that "
        "already passed, not for blessing one that does not.",
    ),
) -> None:
    """Stamp built PONs with the release they ship with.

    A PON is built from `develop`, where the version still reads the previous
    release -- so the model records that, however new the code is. Rather than
    bump the version before a four-hour build, set it here as part of cutting
    the release.

    After this, `krewlyzer_version` means **the release this model is published
    with**, not the code that produced it. That is what a compatibility guard
    needs, and `build_date` still records when it was actually built.

    Only the metadata row changes; every baseline is copied through unchanged
    and the cohort digest is untouched. Re-run `validate-pon` afterwards --
    the file that ships should be the file that was checked.
    """
    from .pon_stamp import stamp_release

    for path in pon_model:
        try:
            previous = stamp_release(path, version, dry_run=dry_run, force=force)
        except ValueError as exc:
            console.print(f"[red]{path.name}: {exc}[/red]")
            raise typer.Exit(2)
        arrow = "would become" if dry_run else "->"
        console.print(f"  {path.name}: {previous or '(unset)'} {arrow} {version}")

    if not dry_run:
        console.print(
            "\n[yellow]Re-run validate-pon on these files[/yellow] -- the file "
            "that ships should be the file that was checked."
        )


def describe_output(
    sample_dir: Path = typer.Argument(
        ...,
        help="One sample directory, e.g. RESULTS/{sample_id}/",
        exists=True,
        file_okay=False,
    ),
    sample_id: Optional[str] = typer.Option(
        None,
        "--sample-id",
        help="Override the sample id (defaults to the directory name)",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help=(
            "Write Markdown here instead of stdout. A .html suffix renders "
            "HTML, which needs the `markdown` package: pip install "
            "'krewlyzer[report]'."
        ),
    ),
) -> None:
    """Describe what a sample's output files contain.

    `validate-output` answers "is this correct?"; this answers the question
    people ask first -- "what are these files, and what is in them?". Shape,
    dtypes, ranges and examples are measured from the sample; which tables are
    gated comes from the output contract. Nothing about the biology is
    restated, so this cannot drift from the reference documentation.
    """
    from .describe import (
        describe_sample,
        markdown_renderer_available,
        render_html_page,
        render_markdown,
    )

    report_obj = describe_sample(sample_dir, sample_id)
    markdown = render_markdown(report_obj)

    if output is None:
        print(markdown)
    elif output.suffix.lower() in (".html", ".htm"):
        # Say it before claiming success. Asking for .html and receiving a page
        # of raw Markdown is confusing precisely because the file *is* HTML --
        # nothing about it looks wrong except the contents.
        if not markdown_renderer_available():
            console.print(
                "[yellow]warning[/yellow] the `markdown` package is not "
                "installed, so the page will contain the Markdown source "
                "rather than rendered HTML.\n"
                "          install it with: pip install 'krewlyzer[report]'"
            )
        output.write_text(render_html_page(report_obj, markdown), encoding="utf-8")
        console.print(f"[green]wrote[/green] {output}")
    else:
        output.write_text(markdown, encoding="utf-8")
        console.print(f"[green]wrote[/green] {output}")

    missing = [t for t in report_obj.missing if t.consumed]
    if missing:
        console.print(
            f"[yellow]{len(missing)} gated table(s) absent[/yellow] "
            "— run validate-output for the contract check"
        )


def report_sample(
    sample_dir: Path = typer.Argument(
        ...,
        help="One sample directory, e.g. RESULTS/{sample_id}/",
        exists=True,
        file_okay=False,
    ),
    output: Path = typer.Option(
        Path("report.html"),
        "--output",
        "-o",
        help="Where to write the HTML report",
    ),
    sample_id: Optional[str] = typer.Option(
        None,
        "--sample-id",
        help="Override the sample id (defaults to the directory name)",
    ),
    z_threshold: float = typer.Option(
        2.0,
        "--z-threshold",
        help="|z| at which an axis is flagged. Conventional, not a clinical "
        "cut-off — the report says so.",
    ),
) -> None:
    """Build a single-sample report: verdict, charts, and the tables behind them.

    For internal use. The output contains one sample's actual measurements, so
    it is generated on demand rather than committed or published; use
    `describe-output` for anything that leaves the machine.
    """
    from .describe import describe_sample
    from .htmlreport import render_html
    from .plots import build_charts
    from .verdict import compute_verdict
    from krewlyzer.core.output_utils import read_table

    described = describe_sample(sample_dir, sample_id)
    resolved_id = described.sample_id

    tables = {
        t.suffix: read_table(sample_dir / f"{resolved_id}{t.suffix}")
        for t in described.tables
    }
    charts = build_charts(tables)
    verdict = compute_verdict(sample_dir, resolved_id, z_threshold)

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_html(described, verdict, charts, tables), encoding="utf-8")

    drawn = sum(1 for c in charts if c.drawn)
    console.print(
        f"[green]wrote[/green] {output}  "
        f"({len(described.present)} tables, {drawn}/{len(charts)} charts, "
        f"{verdict.headline.lower()})"
    )
    undrawn = [c for c in charts if not c.drawn]
    for c in undrawn:
        console.print(f"  [dim]no {c.title.lower()}: {c.reason}[/dim]")
