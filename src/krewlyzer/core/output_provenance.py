"""Record, in the product itself, what produced a sample and what it was scored against.

`{sample}.metadata.parquet` is the downstream consumer's **completion marker**:
a sample without it is dropped from the cohort silently (invariant #2). It
recorded the genome, the assay, the filters and whether GC correction ran — and
nothing at all about the code version or the PON.

The version *was* recorded, in `{sample}.features.json`. Downstream reads
Parquet only and never that file, so the one fact needed to tell two runs apart
lived precisely where the product ignores it.

That became load-bearing in 0.9.0. This release changed what several columns
*mean* — every FSD `_logR` was previously scored against the wrong size bin,
`pon_stability` was wrong by a median of 4709%, and WPS gained seven columns.
Every one of those defects produced output that was present, finite and
plausible. Given a results directory, there was no way to answer *was this
produced before or after the fix?* except by re-running it.

Five fields answer it:

``krewlyzer_version``      which code wrote this sample
``pon_applied``            whether any feature was z-scored at all
``pon_model``              which model, by **basename** — never the full path,
                           which would pin the operator's home directory into a
                           shipped table
``pon_cohort_digest``      which healthy cohort that model was fitted to; the
                           salted digest from `pon/provenance.py`, so it answers
                           "the same cohort?" without naming anyone
``pon_krewlyzer_version``  the release the model was stamped with

`pon_applied` is the one a reader acts on. A `--skip-pon` run and a run with no
PON are both legitimately unscored, and look identical to a reader counting
columns. (A PON that was asked for and cannot be used stops the run instead:
`pon_integration.UnusablePonError`.) Recording
the distinction is also what lets `validate-output` require the PON-derived
columns when — and only when — a PON was actually used.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, Optional

from .output_utils import read_exact_table, resolve_table_path, write_table

logger = logging.getLogger("core.output_provenance")

#: The columns this module adds. Named once so the gate and the tests can both
#: refer to them without restating the list (invariant #5, one level down).
PROVENANCE_COLUMNS = (
    "krewlyzer_version",
    "pon_applied",
    "pon_model",
    "pon_cohort_digest",
    "pon_krewlyzer_version",
)


def provenance_fields(
    pon: Any = None,
    pon_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """The provenance to record for one sample.

    ``pon`` is the loaded model and ``pon_path`` the file it came from; pass
    neither for an unscored run. Both are required for ``pon_applied`` to be
    True: a path without a loaded model is not a scored run.
    """
    from krewlyzer import __version__

    applied = pon is not None and pon_path is not None
    return {
        "krewlyzer_version": __version__,
        "pon_applied": applied,
        # Basename only. The full path would put the operator's home directory
        # into a table that ships, which `scripts/check_phi_guard.sh` exists to
        # keep out of this project.
        "pon_model": Path(pon_path).name if (applied and pon_path) else "",
        "pon_cohort_digest": (
            str(getattr(pon, "cohort_digest", "") or "") if applied else ""
        ),
        "pon_krewlyzer_version": (
            str(getattr(pon, "krewlyzer_version", "") or "") if applied else ""
        ),
    }


def record_path(
    output_dir: Path, sample: str, output_format: str = "tsv", compress: bool = False
) -> Optional[Path]:
    """The metadata file a run of this format reads its record from: its own, else whatever exists.

    One rule for the stamp and the second-PON check, so the two never read different copies (a
    TSV run beside an earlier run's `.parquet`, or the reverse).
    """
    base = output_dir / f"{sample}.metadata"
    own = Path(
        f"{base}.parquet"
        if output_format in ("parquet", "both")
        else f"{base}.tsv.gz" if compress else f"{base}.tsv"
    )
    return own if own.exists() else resolve_table_path(base)


def _flag(value: Any) -> bool:
    """A recorded boolean, with an empty cell false: `bool(NaN)` is True, so an empty
    `pon_applied` read back from TSV counted as applied."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    return bool(value)


def _as_text(value: Any) -> str:
    """A recorded field as text, with None and NaN both empty: an empty field reads back from
    TSV as NaN, which `str()` turns into "nan" and made the same PON look like another.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value)


def refuse_a_second_pon(
    output_dir: Path,
    sample: str,
    pon: Any,
    pon_path: Path,
    output_format: str = "tsv",
    compress: bool = False,
) -> None:
    """Refuse to score a sample against a PON other than the one its record names.

    `pon_applied` is one record per sample. A single tool scoring against another PON would
    leave tables from two PONs under a record that names one, and nothing would flag it.
    """
    from .pon_integration import UnusablePonError

    existing = record_path(output_dir, sample, output_format, compress)
    # As text: a digit-only cohort digest such as '0123…' read back from TSV as an integer.
    text = {"pon_model": str, "pon_cohort_digest": str}
    frame = read_exact_table(existing, dtype=text) if existing is not None else None
    if frame is None or frame.empty or "pon_applied" not in frame.columns:
        return
    if not _flag(frame["pon_applied"].iloc[0]):
        return
    ours = provenance_fields(pon, pon_path)
    recorded = {
        key: frame[key].iloc[0] if key in frame.columns else None
        for key in ("pon_model", "pon_cohort_digest")
    }
    if any(_as_text(recorded[k]) != _as_text(ours[k]) for k in recorded):
        raise UnusablePonError(
            f"{sample} was scored against {recorded['pon_model']} "
            f"(cohort {recorded['pon_cohort_digest']}); refusing to score it against "
            f"{pon_path.name} as well. Re-run run-all to rescore the whole sample, or use "
            "the same PON."
        )


def record_pon_used(
    output_dir: Path,
    sample: str,
    pon: Any,
    pon_path: Optional[Path],
    output_format: str = "tsv",
    compress: bool = False,
) -> bool:
    """A standalone scorer's entry in the sample's record, after it has scored.

    `motif`, `region-entropy` and `region-mds` do not go through `run_features`, so without this
    they checked the record (`refuse_a_second_pon`) but never wrote it: after `extract`, one could
    score with PON A and a later tool with PON B, unnoticed. Raise-only, like every single tool.
    """
    return stamp_metadata(
        output_dir,
        sample,
        pon=pon,
        pon_path=pon_path if pon is not None else None,
        output_format=output_format,
        compress=compress,
        owns_record=False,
    )


def stamp_metadata(
    output_dir: Path,
    sample: str,
    pon: Any = None,
    pon_path: Optional[Path] = None,
    output_format: str = "tsv",
    compress: bool = False,
    gc_factors: Optional[dict] = None,
    owns_record: bool = True,
) -> bool:
    """Add the provenance columns to `{sample}.metadata.*`. True if written.

    `gc_factors` maps "off_target"/"on_target" to the `GcFactorsChoice` the run used; it records
    which kind of factor table each population was corrected with, or why none was (redesign
    item 0: production ran panel WPS, TFBS and ATAC uncorrected and nothing said so).

    A no-op when there is no metadata table, which is the normal case for a
    single-feature CLI run: those write one feature table and no completion
    marker. Logged at debug rather than warned, because it is not a fault.

    `owns_record` is True for run-all and the PON build, which rewrite the whole record. A single
    tool (False) fills it in where none exists, and otherwise may only raise ``pon_applied``, so
    that flag does not depend on the order single tools ran in. Which PON is recorded is kept
    consistent by `refuse_a_second_pon`, which stops a single tool before any work.

    Read-modify-write rather than appended at creation: the metadata table is
    written during extraction, before any PON is loaded, so the PON fields are
    not knowable yet at that point.
    """
    base = output_dir / f"{sample}.metadata"
    # Read the file this run's format writes (`record_path`), then write all of them. Resolving
    # `.tsv` first and then reading Parquet-first meant a TSV run beside an earlier run's
    # `.parquet` read the stale row and wrote it back over this run's.
    existing = record_path(output_dir, sample, output_format, compress)
    if existing is None:
        logger.debug(
            f"No metadata table for {sample}; skipping provenance stamp. "
            "Expected for a single-feature run, which writes no completion marker."
        )
        return False

    frame = read_exact_table(existing)
    if frame is None or frame.empty:
        logger.warning(
            f"{existing.name} could not be read; provenance not recorded. "
            "Downstream cannot tell which build produced this sample."
        )
        return False

    fields = provenance_fields(pon, pon_path)
    recorded = "pon_applied" in frame.columns
    if not owns_record and recorded:
        # A single tool keeps the sample's record (run-all's, or an earlier tool's): rewriting
        # pon_applied and gc_factors_* from one tool's run made validate-output stop checking
        # every other table. It may only *raise* pon_applied, so any order of single tools ends
        # with the same record: a table one of them left unscored is then flagged, not hidden.
        if not (fields["pon_applied"] and not _flag(frame["pon_applied"].iloc[0])):
            logger.info(
                f"{existing.name} already records how this sample was run; a single tool keeps it"
            )
            return False
        gc_factors = (
            None  # raising the PON record says nothing about this tool's GC tables
        )

    for name, value in fields.items():
        frame[name] = value
    for population, choice in (gc_factors or {}).items():
        tag = population.replace("_", "")  # off_target -> offtarget
        frame[f"gc_factors_{tag}"] = choice.record()
        frame[f"gc_correction_applied_{tag}"] = bool(choice.found)

    write_table(frame, base, output_format=output_format, compress=compress)

    if fields["pon_applied"]:
        logger.info(
            f"Provenance: krewlyzer {fields['krewlyzer_version']}, scored against "
            f"{fields['pon_model']} (cohort {fields['pon_cohort_digest'] or 'unrecorded'})"
        )
    else:
        logger.info(
            f"Provenance: krewlyzer {fields['krewlyzer_version']}, no PON applied — "
            "every z-score column is absent by design, not by failure."
        )
    return True
