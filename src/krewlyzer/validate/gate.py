"""Walk an output directory and judge it against the contract.

Two stages, so this scales past a laptop:

* :func:`check_sample` -- everything decidable from one sample (presence,
  schema, domain invariants), plus a small **fingerprint** of each column.
* :func:`evaluate_cohort` -- the cross-sample degeneracy comparison, fed by
  those fingerprints rather than by re-reading the tables.

The split matters because a sample directory is ~1.5 GB (the WPS table alone is
~120 MB) and a real cohort is tens of thousands of samples. Fingerprints are a
hash and a couple of counts per column, so the gather step reads megabytes
instead of terabytes and the scatter step parallelises with no coordination.
:func:`run` glues both together for the single-machine case.
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import checks as check_registry
from . import degeneracy
from .contract import (
    CONTRACT,
    COMPLETION_MARKER,
    Kind,
    NOT_CONSUMED,
    STAGE_TABLES,
    STAGES,
    TableRule,
)
from ..core.output_provenance import TABLE_FORMATS
from ..core.output_utils import TABLE_EXTENSIONS
from .degeneracy import Observation
from .findings import Category, Finding, Severity

logger = logging.getLogger(__name__)

EXIT_PASS = 0
EXIT_VIOLATION = 1
EXIT_STRUCTURAL = 2

FINGERPRINT_VERSION = "1"


@dataclass
class Fingerprint:
    """What a cohort needs to know about one sample. Kilobytes, not gigabytes."""

    sample: str
    observations: Dict[str, Observation] = field(default_factory=dict)
    """Keyed ``"{suffix}::{column}"``."""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "fingerprint_version": FINGERPRINT_VERSION,
            "sample": self.sample,
            "observations": {k: o.to_dict() for k, o in self.observations.items()},
        }

    @classmethod
    def from_dict(cls, raw: Dict[str, Any]) -> "Fingerprint":
        version = raw.get("fingerprint_version")
        if version != FINGERPRINT_VERSION:
            raise ValueError(
                f"fingerprint version {version!r} != {FINGERPRINT_VERSION!r}; "
                "regenerate the per-sample fingerprints"
            )
        return cls(
            sample=raw["sample"],
            observations={
                k: Observation.from_dict(v) for k, v in raw["observations"].items()
            },
        )

    @classmethod
    def load(cls, path: Path) -> "Fingerprint":
        return cls.from_dict(json.loads(path.read_text()))

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2))


@dataclass
class Result:
    findings: List[Finding] = field(default_factory=list)
    samples: List[str] = field(default_factory=list)
    fingerprints: List[Fingerprint] = field(default_factory=list)

    @property
    def exit_code(self) -> int:
        if any(f.category is Category.STRUCTURAL for f in self.findings):
            return EXIT_STRUCTURAL
        if any(f.severity is Severity.ERROR for f in self.findings):
            return EXIT_VIOLATION
        return EXIT_PASS

    def counts(self) -> Dict[str, int]:
        out = {s.value: 0 for s in Severity}
        for f in self.findings:
            out[f.severity.value] += 1
        return out


def discover_samples(results_dir: Path) -> List[Tuple[str, Path]]:
    """Find ``{results_dir}/{sample_id}/`` directories, as kreview does.

    Falls back to treating ``results_dir`` itself as one sample when it holds
    Parquet files directly, so a single run can be checked without reshuffling
    it into a cohort layout.
    """
    found: List[Tuple[str, Path]] = []
    for child in sorted(p for p in results_dir.iterdir() if p.is_dir()):
        if any(child.glob(f"{child.name}.*.parquet")):
            found.append((child.name, child))
    if found:
        return found

    stems = {p.name.split(".", 1)[0] for p in results_dir.glob("*.parquet")}
    if len(stems) == 1:
        return [(stems.pop(), results_dir)]
    return []


def _kind_matches(series: pd.Series, kind: Kind) -> bool:
    non_null = series.dropna()
    if non_null.empty:
        return True
    first = non_null.iloc[0]
    if kind is Kind.LIST:
        return hasattr(first, "__len__") and not isinstance(first, (str, bytes))
    if kind is Kind.NUMERIC:
        return pd.api.types.is_numeric_dtype(series)
    return not pd.api.types.is_numeric_dtype(series)


def _read_table(path: Path, rule: TableRule) -> Tuple[pd.DataFrame, int]:
    """Read only what the checks need, and report the true row count.

    Full materialisation is what makes the naive version unusable on a real
    cohort. Project to the declared columns, and where ``scan_rows`` says the
    checks do not need the whole table, stop after that many rows -- the row
    count still comes from the footer, so the schema check stays honest.
    """
    handle = pq.ParquetFile(path)
    n_rows = handle.metadata.num_rows
    available = set(handle.schema_arrow.names)
    wanted = [c.name for c in rule.columns if c.name in available]

    # `plausible_z_scores` runs on every table, so its columns have to be read
    # even when the contract does not declare them -- otherwise the check sees
    # an absent column and reports nothing, which is the quietest way for a
    # check to do no work.
    wanted = list(
        dict.fromkeys(
            wanted + [c for c in available if c.endswith("_z") or c == "z_score"]
        )
    )

    # Domain checks may read columns the contract does not declare (FSD size
    # bins, FSC channels, *_log2), so they get the whole row instead.
    columns = None if rule.checks else (wanted or None)

    if rule.scan_rows is None:
        return handle.read(columns=columns).to_pandas(), n_rows

    batches = []
    seen = 0
    for batch in handle.iter_batches(
        batch_size=min(rule.scan_rows, 4096), columns=columns
    ):
        batches.append(batch)
        seen += batch.num_rows
        if seen >= rule.scan_rows:
            break
    if not batches:
        return pd.DataFrame(), n_rows
    return pa.Table.from_batches(batches).to_pandas(), n_rows


def _check_table(
    sample: str,
    rule: TableRule,
    df: pd.DataFrame,
    n_rows: int,
    pon_applied: Optional[bool] = None,
    stage_status: Optional[Dict[str, str]] = None,
) -> Tuple[List[Finding], Dict[str, Observation]]:
    family = rule.family.replace(".parquet", "")
    findings: List[Finding] = []
    observations: Dict[str, Observation] = {}

    row_problem = rule.rows.check(n_rows)
    if row_problem:
        findings.append(
            Finding(
                id=f"{family}.ROWS",
                severity=Severity.ERROR,
                category=Category.SCHEMA,
                message=row_problem,
                table=rule.suffix,
                samples=[sample],
                evidence={"n_rows": n_rows},
            )
        )

    for col in rule.columns:
        if col.name not in df.columns:
            # A PON-derived column is required only when the sample says it was
            # scored. Absent provenance (a pre-0.9.0 directory) is `None`, and
            # is not evidence either way, so it stays quiet.
            status = (stage_status or {}).get(col.requires_stage or "")
            # A z-score stage's column is required where the stage ran, and also where a single
            # tool rewrote its table unscored (`not_requested`) in a sample recorded as scored:
            # downstream reads that sample as PON-scored. (`unavailable`, a PON with no baseline
            # for it, and `failed`, reported as STAGE.FAILED, need no column.)
            rewritten_unscored = status == "not_requested" and pon_applied is True
            expected = (
                col.required
                and (not col.requires_pon or pon_applied is True)
                and (
                    col.requires_stage is None or status == "ran" or rewritten_unscored
                )
            )
            if expected:
                if col.requires_pon:
                    why = (
                        " -- this sample records pon_applied=True, so PON scoring "
                        "ran and its output should be here"
                    )
                elif rewritten_unscored:
                    why = (
                        f" -- this sample records pon_applied=True but "
                        f"{col.requires_stage}_status=not_requested: the table was rewritten "
                        "without the PON (a single tool with --skip-pon?). Rerun it with the "
                        "sample's PON, or rerun run-all"
                    )
                elif col.requires_stage:
                    why = f" -- this sample records {col.requires_stage}_status=ran"
                else:
                    why = ""
                findings.append(
                    Finding(
                        id=f"{family}.MISSING_COLUMN.{col.name}",
                        severity=Severity.ERROR,
                        category=Category.SCHEMA,
                        message=f"required column '{col.name}' is absent{why}",
                        table=rule.suffix,
                        column=col.name,
                        samples=[sample],
                        evidence={"columns_present": sorted(map(str, df.columns))},
                    )
                )
            continue
        if not _kind_matches(df[col.name], col.kind):
            findings.append(
                Finding(
                    id=f"{family}.WRONG_KIND.{col.name}",
                    severity=Severity.ERROR,
                    category=Category.SCHEMA,
                    message=(
                        f"column '{col.name}' should be {col.kind.value}, found "
                        f"dtype {df[col.name].dtype}"
                    ),
                    table=rule.suffix,
                    column=col.name,
                    samples=[sample],
                )
            )
            continue
        observations[col.name] = degeneracy.observe(sample, df[col.name], col.kind)

    # `plausible_z_scores` runs on every table, not per rule.
    #
    # It is a divisor check: any `*_z` column reaching |z| = 100 means a sigma
    # that should have been absent. Listing it on 25 rules would work today and
    # be silently missing from the 26th, and it costs nothing where there is no
    # z column to look at.
    for name in (*rule.checks, "plausible_z_scores", "no_collided_columns"):
        fn = check_registry.REGISTRY.get(name)
        if fn is None:
            raise KeyError(f"unknown check '{name}' referenced by {rule.suffix}")
        for problem in fn(df):
            findings.append(
                Finding(
                    id=f"{family}.{name.upper()}",
                    severity=Severity.ERROR,
                    category=Category.DOMAIN,
                    message=problem,
                    table=rule.suffix,
                    samples=[sample],
                )
            )

    return findings, observations


def _pon_applied(sample_dir: Path, sample: str) -> Optional[bool]:
    """Was this sample actually scored against a PON?

    Read from `{sample}.metadata.parquet`, which `run_features` stamps after
    the PON decision is known. Three-valued on purpose:

    ``True``   scored -- the PON-derived columns must be present
    ``False``  deliberately unscored (--skip-pon, no PON, or a PON the version
               guard refused) -- their absence is correct
    ``None``   unrecorded, so a build older than 0.9.0 wrote this directory and
               there is nothing to check against

    `None` is not treated as `False`: they mean different things, and only the
    first is a reason to stay quiet about missing z-scores.
    """
    marker = sample_dir / f"{sample}{COMPLETION_MARKER}"
    if not marker.exists():
        return None
    try:
        frame = pd.read_parquet(marker, columns=["pon_applied"])
    except Exception:
        # Column absent, or the file is unreadable. Either way, unknown --
        # the marker's own absence is already reported separately.
        return None
    if frame.empty or frame["pon_applied"].isna().all():
        return None
    return bool(frame["pon_applied"].iloc[0])


def _stage_outcomes(sample_dir: Path, sample: str) -> Dict[str, Tuple[str, str]]:
    """``{stage: (status, error_class)}`` from `{sample}.metadata.parquet`, for the stages it
    records. Empty for a directory written before outcomes were recorded, which then behaves
    exactly as it did."""
    marker = sample_dir / f"{sample}{COMPLETION_MARKER}"
    if not marker.exists():
        return {}
    try:
        present = set(pq.read_schema(marker).names)
        columns = [
            f"{stage}_{part}"
            for stage in STAGES
            for part in ("status", "error")
            if f"{stage}_{part}" in present
        ]
        frame = pd.read_parquet(marker, columns=columns) if columns else None
    except Exception:
        return {}  # an unreadable marker is reported where the contract reads it
    if frame is None or frame.empty:
        return {}

    def text(column: str) -> str:
        value = frame[column].iloc[0] if column in frame.columns else None
        return value if isinstance(value, str) else ""

    found = {}
    for stage in STAGES:
        status = text(f"{stage}_status")
        if status:
            found[stage] = (status, text(f"{stage}_error"))
    return found


#: Which stage writes each table (`contract.STAGE_TABLES`, inverted).
_STAGE_OF = {
    suffix: stage for stage, suffixes in STAGE_TABLES.items() for suffix in suffixes
}
_CONSUMED = {rule.suffix for rule in CONTRACT}


def _stage_findings(
    sample: str,
    sample_dir: Path,
    outcomes: Dict[str, Tuple[str, str]],
) -> List[Finding]:
    """What the recorded outcomes say, before any table is read.

    A failed stage is one finding with its reason, at ERROR when anything downstream reads one
    of its tables and WARN otherwise -- the severity its absent tables had, so for absent tables
    `--strict-validation` fails the runs it failed before. It fails more runs than before in
    three cases, each deliberate: a failed stage whose tables are only columns (the MDS
    z-scores, which nothing reported), and a consumed table **present** beside a stage that
    failed (`STAGE.STALE_TABLE`: a previous run's table, or a partial one) or that did not run
    in this run (`STAGE.UNEXPECTED_TABLE`: left from another configuration). Downstream reads a
    consumed table as this run's either way. For a table nothing reads, both are WARN.
    """
    findings: List[Finding] = []
    for stage, (status, error) in outcomes.items():
        tables = STAGE_TABLES.get(stage, ()) or tuple(
            rule.suffix
            for rule in CONTRACT
            if any(col.requires_stage == stage for col in rule.columns)
        )
        consumed = any(t in _CONSUMED for t in tables)
        severity = Severity.ERROR if consumed else Severity.WARN
        names = ", ".join(t.lstrip(".").replace(".parquet", "") for t in tables)
        if status == "failed":
            findings.append(
                Finding(
                    id="STAGE.FAILED",
                    severity=severity,
                    category=Category.COMPLETION,
                    message=(
                        f"{stage} was run and failed ({error or 'unrecorded'}); run-all "
                        f"carried on without it, so {names} is missing or incomplete. The "
                        "run's log has the full error"
                    ),
                    table=tables[0] if len(tables) == 1 else None,
                    samples=[sample],
                    evidence={"stage": stage, "error": error, "tables": list(tables)},
                )
            )
        elif status == "unavailable" and not consumed:
            # Its absent tables were an inventory WARN with no reason; this is the reason.
            findings.append(
                Finding(
                    id="STAGE.UNAVAILABLE",
                    severity=Severity.WARN,
                    category=Category.MISSING,
                    message=(
                        f"{stage} was asked for, but an input it needs was missing; "
                        f"{names} was not produced"
                    ),
                    table=tables[0] if len(tables) == 1 else None,
                    samples=[sample],
                    evidence={"stage": stage},
                )
            )
        if stage not in STAGE_TABLES:
            continue  # a column stage: its table belongs to another stage
        for suffix in tables:
            if not (sample_dir / f"{sample}{suffix}").exists():
                continue
            name = suffix.lstrip(".").replace(".parquet", "")
            if status == "failed":
                findings.append(
                    Finding(
                        id="STAGE.STALE_TABLE",
                        severity=(
                            Severity.ERROR if suffix in _CONSUMED else Severity.WARN
                        ),
                        category=Category.COMPLETION,
                        message=(
                            f"{name} is present although {stage} failed in this run: it is "
                            "left from an earlier run, or partial, and would be read as this "
                            "run's"
                        ),
                        table=suffix,
                        samples=[sample],
                        evidence={"stage": stage, "error": error},
                    )
                )
            elif status in ("not_requested", "no_variants", "unavailable"):
                findings.append(
                    Finding(
                        id="STAGE.UNEXPECTED_TABLE",
                        severity=(
                            Severity.ERROR if suffix in _CONSUMED else Severity.WARN
                        ),
                        category=Category.COMPLETION,
                        message=(
                            f"{name} is present although {stage} was {status.replace('_', ' ')} "
                            "in this run: left from an earlier configuration of this directory, "
                            "and would be read as this run's"
                        ),
                        table=suffix,
                        samples=[sample],
                        evidence={"stage": stage, "status": status},
                    )
                )
    return findings


#: Every table the contract or a stage names, and the record, without its extension: the files
#: `FORMAT.MIXED` looks at. BED, JSON, logs and intermediates are not tables anything reads.
_TABLE_BASES = tuple(
    sorted(
        {rule.suffix[: -len(".parquet")] for rule in CONTRACT}
        | {s[: -len(".parquet")] for tables in STAGE_TABLES.values() for s in tables}
        | {COMPLETION_MARKER[: -len(".parquet")]}
    )
)


def _format_findings(sample: str, sample_dir: Path) -> List[Finding]:
    """Tables in an extension outside the record's `table_formats` (`FORMAT.MIXED`).

    The configuration guard refuses a rerun in another format before it writes, except under
    Nextflow: each task starts in an empty folder, and `publishDir` copies a rerun's files over an
    `--outdir` that may hold an earlier format. A table in another extension means a later run
    wrote here in another format, so the Parquet downstream reads may be the earlier run's. A
    record written before `table_formats` existed cannot tell a `both` folder from a mixed one;
    when it holds tables in more than one extension that is said (`FORMAT.UNJUDGED`), not passed.
    """
    marker = sample_dir / f"{sample}{COMPLETION_MARKER}"
    if not marker.exists():
        return []  # COMPLETION.MARKER_ABSENT says so
    try:
        if TABLE_FORMATS in pq.read_schema(marker).names:
            value = pq.read_table(marker, columns=[TABLE_FORMATS]).column(0)[0].as_py()
        else:
            value = None
    except Exception as exc:  # the contract loop reports an unreadable marker
        logger.debug("%s: table_formats unreadable (%s)", marker.name, exc)
        return []
    present: Dict[str, List[str]] = defaultdict(list)
    for base in _TABLE_BASES:
        for ext in TABLE_EXTENSIONS:
            name = f"{sample}{base}{ext}"
            if (sample_dir / name).exists():
                present[ext].append(name)
    if value:
        recorded = value.split(",")
        stray = sorted(
            n for ext, names in present.items() if ext not in recorded for n in names
        )
        if not stray:
            return []
        return [
            Finding(
                id="FORMAT.MIXED",
                severity=Severity.ERROR,
                category=Category.COMPLETION,
                message=(
                    f"{len(stray)} table(s) in another format than this sample's record "
                    f"({value}): a later run wrote here with another --output-format, so the "
                    "Parquet tables downstream reads may be an earlier run's. Rerun the sample "
                    "into a fresh output folder"
                ),
                samples=[sample],
                evidence={"table_formats": value, "files": stray},
            )
        ]
    if len(present) > 1:
        return [
            Finding(
                id="FORMAT.UNJUDGED",
                severity=Severity.SKIP,
                category=Category.COMPLETION,
                message=(
                    f"tables in {', '.join(sorted(present))}, and a record from before "
                    "table_formats was recorded: a --output-format both run and two runs in "
                    "different formats look the same, so this folder could not be checked"
                ),
                samples=[sample],
                evidence={"extensions": sorted(present)},
            )
        ]
    return []


def _absence_explained(suffix: str, outcomes: Dict[str, Tuple[str, str]]) -> bool:
    """Whether a recorded outcome already accounts for this table being absent.

    A failed stage's tables are reported once, as `STAGE.FAILED`. A stage that was not asked for
    (or had no variants, or no input) explains an absent table nothing downstream reads; a
    consumed table keeps its ABSENT ERROR, since downstream reads it either way.
    """
    stage = _STAGE_OF.get(suffix)
    if stage is None or stage not in outcomes:
        return False
    status = outcomes[stage][0]
    if status == "failed":
        return True
    return suffix not in _CONSUMED and status in (
        "not_requested",
        "no_variants",
        "unavailable",
    )


def check_sample(sample: str, sample_dir: Path) -> Tuple[List[Finding], Fingerprint]:
    """Everything decidable from one sample, plus its fingerprint.

    This is the scatter half: it never needs another sample, so a workflow can
    fan it out per sample with no coordination.
    """
    logger.debug("validating %s in %s", sample, sample_dir)
    findings: List[Finding] = []
    fingerprint = Fingerprint(sample=sample)

    marker = sample_dir / f"{sample}{COMPLETION_MARKER}"
    if not marker.exists():
        findings.append(
            Finding(
                id="COMPLETION.MARKER_ABSENT",
                severity=Severity.ERROR,
                category=Category.COMPLETION,
                message=(
                    f"{marker.name} is absent. Consumers use it as the "
                    "completion marker, so this sample is dropped from the "
                    "cohort silently -- no warning, no error"
                ),
                table=COMPLETION_MARKER,
                samples=[sample],
            )
        )

    # Read once, not per table: it comes from the marker and cannot change
    # between tables of the same sample.
    pon_applied = _pon_applied(sample_dir, sample)
    if pon_applied is None:
        logger.debug(
            "%s records no pon_applied; PON-derived columns will not be "
            "required. A build older than 0.9.0 wrote this directory.",
            sample,
        )

    # Also once: what each optional stage did. Empty for an older build's directory.
    outcomes = _stage_outcomes(sample_dir, sample)
    findings.extend(_stage_findings(sample, sample_dir, outcomes))
    findings.extend(_format_findings(sample, sample_dir))
    stage_status = {stage: status for stage, (status, _) in outcomes.items()}

    for rule in CONTRACT:
        path = sample_dir / f"{sample}{rule.suffix}"
        if not path.exists():
            if _absence_explained(rule.suffix, outcomes):
                continue
            findings.append(
                Finding(
                    id=f"{rule.family}.ABSENT".replace(".parquet", ""),
                    severity=Severity.ERROR,
                    category=Category.MISSING,
                    message=f"{rule.suffix.lstrip('.')} is absent",
                    table=rule.suffix,
                    samples=[sample],
                )
            )
            continue
        try:
            df, n_rows = _read_table(path, rule)
        except Exception as exc:  # unreadable == not comparable
            findings.append(
                Finding(
                    id="INPUT.UNREADABLE",
                    severity=Severity.ERROR,
                    category=Category.STRUCTURAL,
                    message=f"could not read {path.name}: {exc}",
                    table=rule.suffix,
                    samples=[sample],
                )
            )
            continue

        table_findings, observations = _check_table(
            sample,
            rule,
            df,
            n_rows,
            pon_applied=pon_applied,
            stage_status=stage_status,
        )
        findings.extend(table_findings)
        for name, observation in observations.items():
            fingerprint.observations[f"{rule.suffix}::{name}"] = observation
        logger.debug(
            "  %s: %d rows, %d column(s) fingerprinted, %d finding(s)",
            rule.suffix.lstrip("."),
            n_rows,
            len(observations),
            len(table_findings),
        )
        del df  # the fingerprint is all the cohort stage needs

    for suffix in NOT_CONSUMED:
        if not (sample_dir / f"{sample}{suffix}").exists():
            if _absence_explained(suffix, outcomes):
                continue
            findings.append(
                Finding(
                    id=f"{suffix.lstrip('.').replace('.parquet', '')}.ABSENT",
                    severity=Severity.WARN,
                    category=Category.MISSING,
                    message=(
                        f"{suffix.lstrip('.')} is absent (not read by "
                        "kreview; reported for inventory only)"
                    ),
                    table=suffix,
                    samples=[sample],
                )
            )

    return findings, fingerprint


def evaluate_cohort(
    fingerprints: Sequence[Fingerprint], min_samples: int = 3
) -> List[Finding]:
    """The gather half: the checks that need more than one sample."""
    collected: Dict[Tuple[str, str], List[Observation]] = defaultdict(list)
    for fp in fingerprints:
        for key, observation in fp.observations.items():
            suffix, _, column = key.partition("::")
            collected[(suffix, column)].append(observation)

    findings: List[Finding] = []
    for rule in CONTRACT:
        for col in rule.columns:
            findings.extend(
                degeneracy.evaluate(
                    rule.suffix,
                    col,
                    collected.get((rule.suffix, col.name), []),
                    min_samples,
                )
            )
    return findings


def reconcile_expected(
    results_dir: Path,
    expected: Sequence[str],
    discovered: Sequence[str],
) -> List[Finding]:
    """Compare the samples you meant to run against the ones that exist.

    Everything else in this module validates the samples it *finds*. Nothing
    can find a sample that produced nothing, and there are two ways for that to
    happen -- both of which leave a cohort quietly smaller than intended:

    * the job never ran, or died before creating its directory;
    * the directory exists but holds no Parquet, so ``discover_samples`` skips
      it. A job killed between mkdir and the first write looks exactly like a
      sample that was never submitted.

    Only the third case -- Parquet present, completion marker absent -- is
    caught today, by ``check_sample``. At a handful of samples the other two are
    obvious. At 16,000 they are invisible, and the consumer reads a short cohort
    without ever knowing it was short.

    The two are reported separately because the remedies differ: an empty
    directory means the job started and its logs exist, while no directory
    usually means it was never submitted.

    ``UNEXPECTED`` is a WARN rather than an ERROR because the usual cause is an
    identifier that does not round-trip -- a suffix in the sheet that the
    pipeline strips, say. That corrupts nothing, but it makes every other count
    here wrong, so it must be visible.
    """
    findings: List[Finding] = []
    want = list(dict.fromkeys(expected))  # de-duplicated, order preserved
    have = set(discovered)

    missing: List[str] = []
    empty: List[str] = []
    for sample in want:
        if sample in have:
            continue
        d = results_dir / sample
        (empty if d.is_dir() else missing).append(sample)

    if missing:
        findings.append(
            Finding(
                id="EXPECTED.NO_OUTPUT_DIRECTORY",
                severity=Severity.ERROR,
                category=Category.COMPLETION,
                message=(
                    f"{len(missing)} of {len(want)} expected sample(s) have no "
                    "output directory at all. The job never ran or died before "
                    "writing anything, and nothing downstream can tell the "
                    "cohort is short"
                ),
                samples=missing,
                evidence={"expected": len(want), "missing": len(missing)},
            )
        )

    if empty:
        findings.append(
            Finding(
                id="EXPECTED.DIRECTORY_HAS_NO_TABLES",
                severity=Severity.ERROR,
                category=Category.COMPLETION,
                message=(
                    f"{len(empty)} expected sample(s) have a directory but no "
                    "Parquet in it, so sample discovery skips them entirely -- "
                    "indistinguishable from never having been submitted"
                ),
                samples=empty,
                evidence={"expected": len(want), "empty": len(empty)},
            )
        )

    unexpected = sorted(have - set(want))
    if unexpected:
        findings.append(
            Finding(
                id="EXPECTED.UNEXPECTED_SAMPLE",
                severity=Severity.WARN,
                category=Category.COMPLETION,
                message=(
                    f"{len(unexpected)} sample(s) produced output but are not in "
                    "the expected list. Usually an identifier that does not "
                    "round-trip between the samplesheet and the pipeline, which "
                    "makes the missing count above unreliable"
                ),
                samples=unexpected,
                evidence={"unexpected": len(unexpected)},
            )
        )

    return findings


def run(
    results_dir: Path,
    min_samples: int = 3,
    only_samples: Optional[Sequence[str]] = None,
    expected_samples: Optional[Sequence[str]] = None,
) -> Result:
    """Both stages on one machine. Nextflow calls the halves separately."""
    result = Result()

    if not results_dir.is_dir():
        result.findings.append(
            Finding(
                id="INPUT.NOT_A_DIRECTORY",
                severity=Severity.ERROR,
                category=Category.STRUCTURAL,
                message=f"{results_dir} is not a directory",
            )
        )
        return result

    samples = discover_samples(results_dir)

    # Before any filtering, and before the empty-cohort exit below: "expected
    # 16,552, found 0" is precisely the case this reconciliation exists for, and
    # returning early on INPUT.NO_SAMPLES would report the least useful half of
    # it. Reconciled against everything discovered, not against --sample-id,
    # which narrows what is checked rather than what was meant to run.
    if expected_samples is not None:
        result.findings.extend(
            reconcile_expected(results_dir, expected_samples, [s for s, _ in samples])
        )

    if only_samples:
        wanted = set(only_samples)
        samples = [(s, p) for s, p in samples if s in wanted]
    if not samples:
        result.findings.append(
            Finding(
                id="INPUT.NO_SAMPLES",
                severity=Severity.ERROR,
                category=Category.STRUCTURAL,
                message=(
                    f"no sample directories found under {results_dir}; expected "
                    "{results_dir}/{sample_id}/{sample_id}.*.parquet"
                ),
            )
        )
        return result

    result.samples = [s for s, _ in samples]
    total = len(samples)
    logger.info("validating %d sample(s) under %s", total, results_dir)
    for index, (sample, sample_dir) in enumerate(samples, start=1):
        findings, fingerprint = check_sample(sample, sample_dir)
        result.findings.extend(findings)
        result.fingerprints.append(fingerprint)
        # Progress matters here: a cohort run reads gigabytes per sample, and a
        # silent hour is indistinguishable from a hang.
        logger.info("[%d/%d] %s: %d finding(s)", index, total, sample, len(findings))

    logger.info(
        "comparing %d fingerprint(s) across the cohort", len(result.fingerprints)
    )
    result.findings.extend(evaluate_cohort(result.fingerprints, min_samples))
    result.findings = _merge(result.findings)
    logger.info(
        "validation complete: %s",
        ", ".join(f"{k}={v}" for k, v in result.counts().items()),
    )
    return result


def _merge(findings: List[Finding]) -> List[Finding]:
    """Collapse the same problem seen in many samples into one finding.

    Per-sample rows make a cohort report unreadable -- 26 samples missing the
    same table is one defect, not 26. The sample list carries the scope.
    """
    merged: Dict[Tuple, Finding] = {}
    for f in findings:
        key = (f.id, f.severity, f.category, f.table, f.column, f.message)
        existing = merged.get(key)
        if existing is None:
            merged[key] = f
        else:
            existing.samples.extend(f.samples)
            for k, v in f.evidence.items():
                existing.evidence.setdefault(k, v)
    return list(merged.values())
