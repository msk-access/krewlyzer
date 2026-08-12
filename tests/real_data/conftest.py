"""A real scored cohort, supplied locally. Never in CI, never committed.

Every defect the PON audit found needed a cohort to see:

| defect                              | what a single fixture could show |
|-------------------------------------|----------------------------------|
| `wps_background` fabricated          | nothing — one sample cannot reveal a constant |
| FSD unnormalised under parquet       | nothing — the fixture PON matches no arm |
| 28.8% single-sample WPS anchors      | nothing — needs per-anchor n across samples |
| `wps_baseline` consumed by nothing   | nothing — the code reads fine |

So this directory exists, and it is deliberately awkward to run: point
``KREWLYZER_TEST_CORPUS`` at a directory of ``{sample_id}/`` subdirectories and
these tests run; leave it unset and they skip. That keeps ``pytest tests/``
green everywhere, keeps CI unaware, and keeps patient data off disk in this
repo.

**PHI.** Sample directories are named for the patient (invariant #4). Nothing
here may put an identifier in a test name, a parameter id, an assertion
message, or a written artifact. Samples are addressed positionally, and
``sample_label`` gives a stable non-reversible handle when one is needed for
output.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

#: Points at a cohort directory laid out as {corpus}/{sample_id}/{sample_id}.*
ENV_VAR = "KREWLYZER_TEST_CORPUS"

#: Below this, cross-sample statements are not worth making -- the same floor
#: `validate-cohort` uses, and for the same reason.
MIN_SAMPLES = 3


def sample_label(directory: Path) -> str:
    """A stable, non-reversible handle for use in messages and artifacts.

    Not the directory name: that is the patient identifier. A truncated hash is
    enough to say "the same sample as the one in the previous line" without
    saying which sample that is.
    """
    return "sample-" + hashlib.sha256(directory.name.encode()).hexdigest()[:8]


@pytest.fixture(scope="session")
def corpus() -> list[Path]:
    """Every sample directory in the cohort, sorted for reproducible ordering."""
    root = os.environ.get(ENV_VAR)
    if not root:
        pytest.skip(f"{ENV_VAR} is not set (local-only gate; see this conftest)")
    path = Path(root).expanduser()
    if not path.is_dir():
        pytest.skip(f"{ENV_VAR}={root} is not a directory")

    samples = sorted(
        d for d in path.iterdir() if d.is_dir() and any(d.glob("*.parquet"))
    )
    if len(samples) < MIN_SAMPLES:
        pytest.skip(
            f"{ENV_VAR} holds {len(samples)} scored sample(s); "
            f"cross-sample checks need at least {MIN_SAMPLES}"
        )
    return samples


@pytest.fixture(scope="session")
def corpus_tables(corpus):
    """``{suffix: [DataFrame, ...]}`` for every table present in the cohort.

    Read once per session: the WPS tables alone are ~44 MB each, and re-reading
    them per test would make this gate slow enough that nobody runs it.
    """
    from krewlyzer.core.output_utils import read_table

    tables: dict[str, list] = {}
    for directory in corpus:
        sample_id = directory.name
        for path in sorted(directory.glob(f"{sample_id}*.parquet")):
            suffix = path.name[len(sample_id) :]
            if ".sync." in suffix:
                continue  # size-resolved companions, not feature tables
            frame = read_table(path)
            if frame is not None:
                tables.setdefault(suffix, []).append(frame)
    return tables
