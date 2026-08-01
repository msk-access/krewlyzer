"""Check a built PON before anything is scored against it.

`validate-output` gates the results of a run. Nothing gated the *reference* a
run measures itself against — and every defect this release fixes was in that
reference, sitting in four shipped files for months:

- `wps_background` held a hardcoded `167.0 / 5.0 / 0.0 / 1.0`, identical across
  all 28 groups and all four models, from cohorts of 21 and 47 samples.
- Six σ floors turned "no spread measured" into a divisor.
- Four blocks were built, shipped, and read by nothing.
- No model records what it was built from, so none can be reproduced.

Each is visible in the file itself. None was visible to any check.

The load-bearing one is degeneracy, and it is the same invariant the output
gate enforces one level up: **a baseline that cannot vary with its cohort is
worse than a missing one**, because it is present, plausible, and produces
z-scores that survive every schema check. `wps_background` is the proof — a
single assertion that σ differs between groups would have caught it in March.

Exit codes match `validate-output`: 0 satisfied, 1 violation, 2 structural.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from .findings import Category, Finding, Severity

logger = logging.getLogger("validate.pon")

#: Blocks whose entries carry a mean/σ pair, and the column prefix to check.
#:
#: Read from the parquet rather than the model so a block nothing parses yet
#: is still checked — the point is to see what is *in the file*.
SCALAR_BLOCKS = {
    "wps_background": ("group_id", ("nrl", "periodicity")),
    "region_mds": ("gene", ("mds", "mds_e1")),
    "region_mds_exon": ("gene", ("mds",)),
    "fsc_gene_baseline": ("gene", ("depth",)),
    "fsc_region_baseline": ("region_id", ("depth",)),
    "ocf_baseline": ("region_id", ("ocf",)),
    "tfbs_baseline": ("label", ("entropy",)),
    "atac_baseline": ("label", ("entropy",)),
    "wps_shape_baseline": ("region_id", ("log_amplitude", "shape_corr_fisher")),
}

#: Below this, a baseline entry is an anecdote. Matches MIN_SAMPLES_PER_KEY in
#: the builder; asserted here because the two are enforced in different places.
MIN_SAMPLES = 3


def _std_columns(frame: pd.DataFrame, prefixes) -> List[str]:
    return [f"{p}_std" for p in prefixes if f"{p}_std" in frame.columns]


def check_pon(path: Path) -> List[Finding]:
    """Everything wrong with a PON that the file itself can reveal."""
    findings: List[Finding] = []

    try:
        table = pd.read_parquet(path)
    except Exception as exc:
        return [
            Finding(
                id="PON.UNREADABLE",
                severity=Severity.ERROR,
                category=Category.STRUCTURAL,
                message=f"could not read {path.name}: {exc}",
            )
        ]

    if "table" not in table.columns:
        return [
            Finding(
                id="PON.NOT_A_PON",
                severity=Severity.ERROR,
                category=Category.STRUCTURAL,
                message=f"{path.name} has no `table` column; is it a PON?",
            )
        ]

    findings.extend(_check_provenance(table, path))
    findings.extend(_check_blocks(table))
    return findings


def _check_provenance(table: pd.DataFrame, path: Path) -> List[Finding]:
    """A model that cannot be reproduced cannot be audited."""
    found: List[Finding] = []
    meta = table[table["table"] == "metadata"]
    if meta.empty:
        return [
            Finding(
                id="PON.NO_METADATA",
                severity=Severity.ERROR,
                category=Category.STRUCTURAL,
                message=f"{path.name} has no metadata row",
            )
        ]
    row = meta.iloc[0]

    if not str(row.get("krewlyzer_version", "") or "").strip():
        found.append(
            Finding(
                id="PON.NO_VERSION",
                severity=Severity.ERROR,
                category=Category.SCHEMA,
                message=(
                    "no krewlyzer_version recorded. 0.9.0 changes what every "
                    "feature means, so a model that cannot say which version "
                    "built it cannot be checked for compatibility. Rebuild "
                    "with build-pon."
                ),
            )
        )
    if not str(row.get("cohort_digest", "") or "").strip():
        found.append(
            Finding(
                id="PON.NO_COHORT",
                severity=Severity.WARN,
                category=Category.SCHEMA,
                message=(
                    "no cohort_digest recorded, so this model cannot be "
                    "reproduced or compared against another build"
                ),
            )
        )
    n = row.get("n_samples")
    if pd.isna(n) or float(n) < MIN_SAMPLES:
        found.append(
            Finding(
                id="PON.TOO_FEW_SAMPLES",
                severity=Severity.ERROR,
                category=Category.DEGENERACY,
                message=f"built from {n} samples; below the floor of {MIN_SAMPLES}",
                evidence={"n_samples": None if pd.isna(n) else float(n)},
            )
        )
    return found


def _check_blocks(table: pd.DataFrame) -> List[Finding]:
    """Per-block: is this fitted, and can it vary?"""
    found: List[Finding] = []
    for name, (key, prefixes) in SCALAR_BLOCKS.items():
        block = table[table["table"] == name]
        if block.empty:
            continue
        for column in _std_columns(block, prefixes):
            values = pd.to_numeric(block[column], errors="coerce")
            finite = values[np.isfinite(values)]

            if finite.empty:
                found.append(
                    Finding(
                        id="PON.BLOCK_UNFITTED",
                        severity=Severity.WARN,
                        category=Category.DEGENERACY,
                        message=(
                            f"{name}.{column} has no finite value in "
                            f"{len(block)} entries: nothing can be z-scored "
                            "against it"
                        ),
                        table=name,
                        column=column,
                    )
                )
                continue

            # The wps_background signature. One sigma repeated across every
            # entry is what a hardcoded default looks like from the outside,
            # and it is indistinguishable from a fitted one by any other check.
            if len(finite) > 1 and finite.nunique() == 1:
                found.append(
                    Finding(
                        id="PON.BLOCK_DEGENERATE",
                        severity=Severity.ERROR,
                        category=Category.DEGENERACY,
                        message=(
                            f"{name}.{column} is {finite.iloc[0]!r} for all "
                            f"{len(finite)} entries. A baseline that cannot "
                            "vary with its cohort was not fitted to one."
                        ),
                        table=name,
                        column=column,
                        evidence={"value": float(finite.iloc[0]), "n": len(finite)},
                    )
                )

            if (finite <= 0).any():
                found.append(
                    Finding(
                        id="PON.NONPOSITIVE_SIGMA",
                        severity=Severity.ERROR,
                        category=Category.DOMAIN,
                        message=(
                            f"{name}.{column} has {(finite <= 0).sum()} "
                            "non-positive value(s); a z-score divided by one "
                            "is infinite, not conservative"
                        ),
                        table=name,
                        column=column,
                    )
                )

        if "n_samples" in block.columns:
            counts = pd.to_numeric(block["n_samples"], errors="coerce").dropna()
            thin = int((counts < MIN_SAMPLES).sum())
            if thin:
                found.append(
                    Finding(
                        id="PON.THIN_ENTRIES",
                        severity=Severity.ERROR,
                        category=Category.DEGENERACY,
                        message=(
                            f"{name}: {thin}/{len(counts)} entries are backed "
                            f"by fewer than {MIN_SAMPLES} samples"
                        ),
                        table=name,
                        evidence={"thin": thin, "total": len(counts)},
                    )
                )
    return found


def exit_code(findings: List[Finding]) -> int:
    """2 for structural, 1 for a violation, 0 otherwise — as validate-output."""
    if any(f.category is Category.STRUCTURAL for f in findings):
        return 2
    if any(f.severity is Severity.ERROR for f in findings):
        return 1
    return 0


def describe(path: Path) -> Optional[str]:
    """One line of provenance for a log, or None if there is none to give."""
    try:
        meta = pd.read_parquet(path).query("table == 'metadata'")
    except Exception:
        return None
    if meta.empty:
        return None
    row = meta.iloc[0]
    version = str(row.get("krewlyzer_version", "") or "unknown")
    digest = str(row.get("cohort_digest", "") or "unrecorded")
    label = str(row.get("cohort_label", "") or "")
    suffix = f" [{label}]" if label else ""
    return f"built by krewlyzer {version}, cohort {digest}{suffix}"
