"""Walk an output directory and judge it against the contract."""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from . import checks as check_registry
from . import degeneracy
from .contract import CONTRACT, COMPLETION_MARKER, Kind, NOT_CONSUMED, TableRule
from .findings import Category, Finding, Severity

logger = logging.getLogger(__name__)

EXIT_PASS = 0
EXIT_VIOLATION = 1
EXIT_STRUCTURAL = 2


@dataclass
class Result:
    findings: List[Finding] = field(default_factory=list)
    samples: List[str] = field(default_factory=list)

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
    if kind is Kind.ANY:
        return True
    non_null = series.dropna()
    if non_null.empty:
        return True
    first = non_null.iloc[0]
    if kind is Kind.LIST:
        return hasattr(first, "__len__") and not isinstance(first, (str, bytes))
    if kind is Kind.NUMERIC:
        return pd.api.types.is_numeric_dtype(series)
    return not pd.api.types.is_numeric_dtype(series)


def _check_table(
    sample: str, rule: TableRule, df: pd.DataFrame
) -> Tuple[List[Finding], Dict[str, pd.Series]]:
    family = rule.family.replace(".parquet", "")
    findings: List[Finding] = []
    columns_for_degeneracy: Dict[str, pd.Series] = {}

    row_problem = rule.rows.check(len(df))
    if row_problem:
        findings.append(
            Finding(
                id=f"{family}.ROWS",
                severity=Severity.ERROR,
                category=Category.SCHEMA,
                message=row_problem,
                table=rule.suffix,
                samples=[sample],
                evidence={"n_rows": len(df)},
            )
        )

    for col in rule.columns:
        if col.name not in df.columns:
            if col.required:
                findings.append(
                    Finding(
                        id=f"{family}.MISSING_COLUMN.{col.name}",
                        severity=Severity.ERROR,
                        category=Category.SCHEMA,
                        message=f"required column '{col.name}' is absent",
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
        columns_for_degeneracy[col.name] = df[col.name]

    for name in rule.checks:
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

    return findings, columns_for_degeneracy


def run(
    results_dir: Path,
    min_samples: int = 3,
    only_samples: Optional[Sequence[str]] = None,
) -> Result:
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

    # column values keyed by (suffix, column) -> [(sample, series), ...]
    collected: Dict[Tuple[str, str], List[Tuple[str, pd.Series]]] = defaultdict(list)

    for sample, sample_dir in samples:
        marker = sample_dir / f"{sample}{COMPLETION_MARKER}"
        if not marker.exists():
            result.findings.append(
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

        for rule in CONTRACT:
            path = sample_dir / f"{sample}{rule.suffix}"
            if not path.exists():
                result.findings.append(
                    Finding(
                        id=f"{rule.family}.ABSENT".replace(".parquet", ""),
                        severity=Severity.WARN if rule.optional else Severity.ERROR,
                        category=Category.MISSING,
                        message=f"{rule.suffix.lstrip('.')} is absent",
                        table=rule.suffix,
                        samples=[sample],
                    )
                )
                continue
            try:
                df = pd.read_parquet(path)
            except Exception as exc:  # unreadable == not comparable
                result.findings.append(
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

            table_findings, columns = _check_table(sample, rule, df)
            result.findings.extend(table_findings)
            for name, series in columns.items():
                collected[(rule.suffix, name)].append((sample, series))

        for suffix in NOT_CONSUMED:
            path = sample_dir / f"{sample}{suffix}"
            if not path.exists():
                result.findings.append(
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

    for rule in CONTRACT:
        for col in rule.columns:
            per_sample = collected.get((rule.suffix, col.name), [])
            result.findings.extend(
                degeneracy.evaluate(rule.suffix, col, per_sample, min_samples)
            )

    result.findings = _merge(result.findings)
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
