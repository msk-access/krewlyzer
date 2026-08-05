"""What a PON was built from, recorded without recording who.

The four models in this repository carry `n_samples` and nothing else. There
is no way to tell which samples produced them, whether two were built from the
same cohort, or whether a rebuild reproduced the last one. `n_samples: 21` is a
single integer with nothing behind it — and the build script that produced it
has a stale header comment claiming 47, so even the integer is not
corroborated anywhere.

That is a problem for a file whose whole job is to be the reference every
z-score is measured against.

## Why not just store the sample list

Sample directories are named for the patient (invariant #4). A PON ships in
this repository, in the Docker image and on PyPI, so it is the last place a
list of identifiers may appear.

A **salted hash** gives the useful half without the dangerous half: two builds
from the same cohort produce the same digest, a build from a different cohort
produces a different one, and the digest reveals nothing about who is in it.
The salt is fixed and public — its purpose is domain separation, so a digest
here cannot be compared against one computed elsewhere, not secrecy. Identifier
spaces are small enough that an unsalted hash of a known ID list is reversible
by enumeration; this is why the salt is not optional.

What you can answer with the digest: *is this the cohort I think it is?* What
you cannot: *who is in it?* That is the intended split.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable, Optional

#: Domain separator. Public by design, and stable: changing it changes every
#: digest and silently breaks "same cohort?" comparisons against older models.
COHORT_SALT = b"krewlyzer-pon-cohort-v1"

#: How much of the digest to keep. 16 hex characters is 64 bits -- ample to
#: distinguish cohorts, short enough to read in a log line.
DIGEST_CHARS = 16


def cohort_digest(sample_paths: Iterable[Path]) -> str:
    """A stable, non-reversible fingerprint of the cohort.

    Derived from the sample *stems*, sorted and de-duplicated, so the digest
    survives the cohort being moved between filesystems or re-run from a
    different working directory — which a path-based hash would not, making it
    useless for the one question it exists to answer.
    """
    stems = sorted({Path(p).name.split(".")[0] for p in sample_paths})
    if not stems:
        return ""
    digest = hashlib.sha256(COHORT_SALT)
    for stem in stems:
        digest.update(b"\x00")
        digest.update(stem.encode("utf-8"))
    return digest.hexdigest()[:DIGEST_CHARS]


def build_provenance(
    sample_paths: Iterable[Path],
    krewlyzer_version: str,
    cohort_label: Optional[str] = None,
    input_kind: str = "",
) -> dict:
    """The provenance fields to write into a PON's metadata row.

    ``krewlyzer_version`` is the one that matters most: 0.9.0 changes what
    every feature *means*, so a PON built by an earlier version is not merely
    old, it is measuring something else. Recording it is what lets the loader
    refuse one (Phase C3) instead of silently producing wrong z-scores.

    ``input_kind`` records what the cohort was made of -- ``"bam"``,
    ``"bed"``, ``"mixed"`` or ``"outputs"``. Without it the gate cannot tell a
    block that was never asked for from one that failed: ``mds_baseline`` and
    ``region_mds`` need a BAM, so their absence is legitimate for a fragment-BED
    cohort and a defect for a BAM one. It stayed a warning for both until this
    field existed. Empty for models built before it did.
    """
    return {
        "krewlyzer_version": krewlyzer_version,
        "cohort_digest": cohort_digest(sample_paths),
        "cohort_label": cohort_label or "",
        "input_kind": input_kind,
    }
