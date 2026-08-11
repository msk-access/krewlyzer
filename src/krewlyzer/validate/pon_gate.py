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

#: Blocks whose σ lives in per-bin columns rather than one ``*_std`` column.
#:
#: ``fsd_baseline`` stores one row per (arm, bin) and ``gc_bias`` one row per GC
#: bin, so neither fits `SCALAR_BLOCKS`. Both were therefore unchecked -- and
#: both are applied: FSD supplies the log-ratios, and `gc_bias` reaches FSC and
#: FSR through ``PonModel.get_mean`` / ``get_variance``.
VECTOR_BLOCKS = {
    "fsd_baseline": ("arm", ("std",)),
    "fsd_baseline_ontarget": ("arm", ("std",)),
    "gc_bias": (
        "gc_bin",
        ("short_std", "intermediate_std", "long_std"),
    ),
    "gc_bias_ontarget": (
        "gc_bin",
        ("short_std", "intermediate_std", "long_std"),
    ),
}

#: Blocks every PON must carry, whatever it was built from.
#:
#: All of these are computed from the fragment BED, so no input choice explains
#: their absence. Missing means something failed quietly.
CORE_BLOCKS = (
    "gc_bias",
    "fsd_baseline",
    "wps_baseline",
    "wps_shape_baseline",
    "wps_background",
    "ocf_baseline",
    "tfbs_baseline",
    "atac_baseline",
)

#: Blocks that need BAM/CRAM input, so a BED-built PON legitimately lacks them.
#:
#: Warned rather than failed: the metadata records `panel_mode` but not whether
#: the cohort was BAMs or fragment BEDs, so the gate cannot tell "not asked for"
#: from "went wrong". Saying so is better than guessing either way.
BAM_ONLY_BLOCKS = (
    "mds_baseline",
    # Breakpoint 4-mers are a different distribution from end 4-mers, so they
    # need their own block; scoring BreakPointMotif against `mds_baseline` gave
    # median |z| 5.85 where a fitted baseline gives ~0.67. Every PON built
    # before 0.9.0 lacks it, and a warning is the right level: those models are
    # not broken, they simply predate the block, and the scorer emits no
    # `frequency_z` rather than a wrong one.
    "breakpoint_motif_baseline",
    "region_mds",
    "region_mds_exon",
)

#: Additionally required when ``panel_mode`` is set.
PANEL_BLOCKS = (
    "gc_bias_ontarget",
    "fsd_baseline_ontarget",
    "wps_baseline_panel",
    "ocf_baseline_ontarget",
    "ocf_baseline_offtarget",
    "fsc_gene_baseline",
    "fsc_region_baseline",
    "tfbs_baseline_ontarget",
    "atac_baseline_ontarget",
)

#: Panel blocks that also need BAM/CRAM input.
PANEL_BAM_ONLY_BLOCKS = (
    "mds_baseline_ontarget",
    "breakpoint_motif_baseline_ontarget",
)

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
    findings.extend(_check_required_blocks(table))
    findings.extend(_check_blocks(table))
    return findings


def _check_required_blocks(table: pd.DataFrame) -> List[Finding]:
    """Is everything that should be here, here?

    Every other check reads the blocks *present in the file* and skips the rest,
    so a block that vanished entirely is invisible to all of them -- it looks
    exactly like one that was never expected. Demonstrated on a real model:
    deleting `region_mds`, or `fsc_gene_baseline`, or almost every baseline at
    once, produced an identical finding list each time.

    That is not hypothetical either. `region_mds` really did vanish, because a
    Rust reader was handed a format it could not parse, returned nothing, and
    the builder logged one warning and carried on. Build, gate and downstream
    all reported success on a PON with no per-gene MDS in it.

    This is the packing list: the check that notices a box never arrived.
    """
    found: List[Finding] = []
    present = set(table["table"].dropna().unique())

    meta = table[table["table"] == "metadata"]
    panel_mode = False
    if not meta.empty:
        raw = meta.iloc[0].get("panel_mode")
        panel_mode = bool(raw) and str(raw).lower() not in ("false", "0", "nan")

    expected_error = list(CORE_BLOCKS) + (list(PANEL_BLOCKS) if panel_mode else [])
    expected_warn = list(BAM_ONLY_BLOCKS) + (
        list(PANEL_BAM_ONLY_BLOCKS) if panel_mode else []
    )

    missing = [b for b in expected_error if b not in present]
    if missing:
        found.append(
            Finding(
                id="PON.BLOCK_MISSING",
                severity=Severity.ERROR,
                # DEGENERACY, not STRUCTURAL. Structural means "is this even a
                # PON" and exits 2; a model that parses fine but is short a
                # block is a contract violation like any other, and exits 1.
                # The bundled 0.8.x models are exactly that -- they predate
                # `wps_shape_baseline` -- and they should fail the same way
                # their fabricated sigma already makes them fail.
                category=Category.DEGENERACY,
                message=(
                    f"absent from the model: {', '.join(missing)}. These are "
                    "built from the fragment BED, so no input choice explains "
                    "them being gone -- something failed quietly during the "
                    "build. Nothing else here can see a block that is not in "
                    "the file."
                ),
                evidence={"missing": missing, "panel_mode": panel_mode},
            )
        )

    missing_bam = [b for b in expected_warn if b not in present]
    if missing_bam:
        # `input_kind` decides the severity. A fragment-BED cohort has no
        # k-mers and no per-gene MDS, so their absence is expected; a BAM
        # cohort should have produced them, and their absence is the exact
        # failure that let `region_mds` disappear behind one warning line.
        #
        # Models built before the field exists record nothing, and stay a
        # warning -- guessing either way would be worse than saying so.
        kind = ""
        if not meta.empty:
            kind = str(meta.iloc[0].get("input_kind", "") or "").strip().lower()
        from_bams = kind in ("bam", "mixed", "outputs")

        found.append(
            Finding(
                id="PON.BLOCK_MISSING_BAM_ONLY",
                severity=Severity.ERROR if from_bams else Severity.WARN,
                category=Category.DEGENERACY if from_bams else Category.SCHEMA,
                message=(
                    f"absent from the model: {', '.join(missing_bam)}. "
                    + (
                        f"The cohort was {kind!r}, which should have produced "
                        "them -- check the build log for region-mds warnings."
                        if from_bams
                        else "These need BAM/CRAM input. This model does not "
                        "record what its cohort was made of (built before "
                        "`input_kind`), so a fragment-BED cohort and a failed "
                        "build look the same here."
                    )
                ),
                evidence={"missing": missing_bam, "input_kind": kind},
            )
        )
    return found


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
    """Per-block: is this fitted, and can it vary?

    `VECTOR_BLOCKS` join `SCALAR_BLOCKS` here. They were left out originally
    because their σ sits in per-bin columns rather than one ``*_std`` -- but
    both are applied, so being awkward to check is not a reason not to. A
    degenerate `gc_bias` is the worst of the lot: it feeds every FSC log-ratio
    through ``PonModel.get_mean``.
    """
    found: List[Finding] = []
    checked = {
        **{n: (k, p, False) for n, (k, p) in SCALAR_BLOCKS.items()},
        **{n: (k, p, True) for n, (k, p) in VECTOR_BLOCKS.items()},
    }
    for name, (key, prefixes, is_vector) in checked.items():
        block = table[table["table"] == name]
        if block.empty:
            continue
        # Vector blocks name their σ columns outright; scalar ones take a
        # prefix and have `_std` appended.
        columns = (
            [c for c in prefixes if c in block.columns]
            if is_vector
            else _std_columns(block, prefixes)
        )
        for column in columns:
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
