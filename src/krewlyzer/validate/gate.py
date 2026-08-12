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
from .contract import CONTRACT, COMPLETION_MARKER, Kind, NOT_CONSUMED, TableRule
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
            expected = col.required and (not col.requires_pon or pon_applied is True)
            if expected:
                why = (
                    " -- this sample records pon_applied=True, so PON scoring "
                    "ran and its output should be here"
                    if col.requires_pon
                    else ""
                )
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

    for rule in CONTRACT:
        path = sample_dir / f"{sample}{rule.suffix}"
        if not path.exists():
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
            sample, rule, df, n_rows, pon_applied=pon_applied
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
