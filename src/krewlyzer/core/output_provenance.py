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

`pon_applied` is the one a reader acts on. A `--skip-pon` run, a run with no
PON, and a run whose PON the version guard refused are all legitimately
unscored, and all three look identical to a reader counting columns. Recording
the distinction is also what lets `validate-output` require the PON-derived
columns when — and only when — a PON was actually used.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

from .output_utils import read_table, resolve_table_path, write_table

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
    True — a path without a loaded model is exactly the refused-PON case, and
    it is not a scored run.
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


def stamp_metadata(
    output_dir: Path,
    sample: str,
    pon: Any = None,
    pon_path: Optional[Path] = None,
    output_format: str = "tsv",
    compress: bool = False,
) -> bool:
    """Add the provenance columns to `{sample}.metadata.*`. True if written.

    A no-op when there is no metadata table, which is the normal case for a
    single-feature CLI run: those write one feature table and no completion
    marker. Logged at debug rather than warned, because it is not a fault.

    Read-modify-write rather than appended at creation: the metadata table is
    written during extraction, before any PON is loaded, so the PON fields are
    not knowable yet at that point.
    """
    base = output_dir / f"{sample}.metadata"
    existing = resolve_table_path(base)
    if existing is None:
        logger.debug(
            f"No metadata table for {sample}; skipping provenance stamp. "
            "Expected for a single-feature run, which writes no completion marker."
        )
        return False

    frame = read_table(existing)
    if frame is None or frame.empty:
        logger.warning(
            f"{existing.name} could not be read; provenance not recorded. "
            "Downstream cannot tell which build produced this sample."
        )
        return False

    fields = provenance_fields(pon, pon_path)
    for name, value in fields.items():
        frame[name] = value

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
