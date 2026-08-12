"""Stamp a built PON with the release it ships with.

A PON is built from `develop`, where `pyproject.toml` still reads the previous
version — so the model records `0.8.3` however new the code is. Bumping the
version before the rebuild would fix that, at the cost of putting a release
number on unreleased code and blocking a four-hour build on a release chore.

Restamping at release time is the other order: build from develop, then set the
field to the tag as part of cutting the release.

## What the field means afterwards

``krewlyzer_version`` becomes **the release this model is published with**, not
the code that produced it. That is the definition a compatibility guard needs —
"can this model be scored against by that release" — and it is the one worth
recording, since a PON's compatibility tracks feature semantics, which change
at releases.

``build_date`` is left alone, so the two together still say when it was built.

Only the metadata row is touched. Every baseline is copied through byte for
byte, and the cohort digest — a hash of the sample IDs — is unaffected.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger("validate.pon_stamp")


def stamp_release(
    path: Path, version: str, dry_run: bool = False, force: bool = False
) -> str:
    """Set ``krewlyzer_version`` on a PON. Returns the previous value.

    **Refuses to stamp a model that fails ``validate-pon``.** Without that,
    this command would be the shortest path to laundering: run it on one of
    the four models carrying the fabricated ``167.0 / 5.0`` baseline and it
    claims 0.9.0 compatibility, which is precisely what the guard exists to
    prevent. A stamp is an assertion that the model is fit to score against,
    so it has to be earned.

    ``force`` exists for the one legitimate exception — re-stamping a model
    the gate already passed, e.g. correcting a typo in the version — and says
    what it is doing.

    Raises ``ValueError`` if the file has no metadata row, rather than writing
    one: a PON without metadata is malformed, and inventing the row hides that.
    """
    if not force:
        from .pon_gate import check_pon, exit_code

        findings = check_pon(path)
        # NO_VERSION is expected here -- it is the thing being fixed.
        blocking = [
            f
            for f in findings
            if f.id != "PON.NO_VERSION" and f.severity.value == "error"
        ]
        if blocking:
            raise ValueError(
                f"{path.name} fails validate-pon and will not be stamped:\n"
                + "\n".join(f"  {f.id}: {f.message}" for f in blocking[:4])
                + "\nStamping it would assert a compatibility it does not have."
            )
        del exit_code  # only imported alongside check_pon for clarity

    table = pd.read_parquet(path)
    if "table" not in table.columns:
        raise ValueError(f"{path.name} has no `table` column; is it a PON?")

    meta_rows = table.index[table["table"] == "metadata"]
    if len(meta_rows) == 0:
        raise ValueError(f"{path.name} has no metadata row to stamp")
    if len(meta_rows) > 1:
        raise ValueError(f"{path.name} has {len(meta_rows)} metadata rows, expected 1")

    row = meta_rows[0]
    if "krewlyzer_version" not in table.columns:
        # Built before provenance existed. Adding the column is right; what
        # would be wrong is doing it without the gate check above.
        table["krewlyzer_version"] = ""
    previous = str(table.at[row, "krewlyzer_version"] or "")

    if dry_run:
        logger.info(f"{path.name}: {previous or '(unset)'} -> {version} (dry run)")
        return previous

    table.at[row, "krewlyzer_version"] = version
    # Rewrite whole: parquet is columnar and immutable in place, so there is no
    # partial update. Every other row is unchanged by construction.
    table.to_parquet(path, index=False)
    logger.info(f"{path.name}: krewlyzer_version {previous or '(unset)'} -> {version}")
    return previous
