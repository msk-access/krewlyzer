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
        help="Write Markdown here instead of stdout. A .html suffix renders HTML.",
    ),
) -> None:
    """Describe what a sample's output files contain.

    `validate-output` answers "is this correct?"; this answers the question
    people ask first -- "what are these files, and what is in them?". Shape,
    dtypes, ranges and examples are measured from the sample; which tables are
    read downstream comes from the output contract. Nothing about the biology
    is restated, so this cannot drift from the reference documentation.
    """
    from .describe import describe_sample, render_markdown

    report_obj = describe_sample(sample_dir, sample_id)
    markdown = render_markdown(report_obj)

    if output is None:
        print(markdown)
    elif output.suffix.lower() in (".html", ".htm"):
        output.write_text(_render_html(report_obj, markdown), encoding="utf-8")
        console.print(f"[green]wrote[/green] {output}")
    else:
        output.write_text(markdown, encoding="utf-8")
        console.print(f"[green]wrote[/green] {output}")

    missing = [t for t in report_obj.missing if t.consumed]
    if missing:
        console.print(
            f"[yellow]{len(missing)} table(s) read downstream are absent[/yellow]"
        )


def _render_html(report_obj, markdown: str) -> str:
    """Self-contained HTML, so the report can be attached or hosted as-is.

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
<title>krewlyzer output — {report_obj.sample_id}</title>
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
