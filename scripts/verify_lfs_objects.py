#!/usr/bin/env python3
"""Assert every Git LFS object is present and is the object git recorded.

An LFS pointer stores the sha256 of the file it stands for, so hashing what is
on disk and comparing it to the pointer's oid is an exact identity check --
not a heuristic. That single property makes this script serve two callers:

  1. After `git lfs pull`, it answers "did the objects actually hydrate?".
     Absence is not loud on its own: `tests/conftest.py` marks 30 tests
     `requires_data`, so a checkout without data reports them as SKIPPED, and a
     skip count does not fail a build. A partial fetch is worse -- the paths
     exist as ~130-byte pointer files, so existence checks pass and the failure
     surfaces later as "Corrupt footer" from inside a reader.

  2. After recovering the data from the published container image, it answers
     "is the image's copy the same data this branch expects?". The image is
     built at release time, so it can lag a branch that changed the data. Here
     that is not a risk to reason about: the hashes either match, in which case
     the bytes are provably identical to what LFS would have supplied, or they
     do not and this exits non-zero.

Exit 0 when every object checks out, 1 otherwise.
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

# A pointer file is a three-line text stub; the smallest real object here is
# orders of magnitude larger. Used only to label the failure, never to pass it.
POINTER_SENTINEL = b"version https://git-lfs"
CHUNK = 1 << 20


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def tracked_objects() -> list[tuple[str, Path]]:
    """(oid, path) for every LFS-tracked file, read from the pointers.

    `git lfs ls-files -l` reads the pointer files rather than the objects, so
    it works before anything is hydrated -- which is the whole point, since
    that is exactly when this script needs to run.
    """
    result = subprocess.run(
        ["git", "lfs", "ls-files", "-l"],
        capture_output=True,
        text=True,
        check=True,
    )
    rows = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        # "<oid> <*|-> <path>" -- the marker is hydration state, which we do
        # not trust; the hash below is the authority.
        oid, _marker, path = line.split(None, 2)
        rows.append((oid, Path(path)))
    return rows


def main() -> int:
    objects = tracked_objects()
    if not objects:
        print("no LFS-tracked files found -- is this a git checkout?", file=sys.stderr)
        return 1

    problems: list[str] = []
    for oid, path in objects:
        if not path.is_file():
            problems.append(f"missing  {path}")
            continue
        actual = sha256_of(path)
        if actual == oid:
            continue
        head = path.open("rb").read(len(POINTER_SENTINEL))
        if head == POINTER_SENTINEL:
            problems.append(
                f"pointer  {path} ({path.stat().st_size} bytes, never hydrated)"
            )
        else:
            problems.append(
                f"mismatch {path}\n           expected {oid}\n           found    {actual}"
            )

    if problems:
        print(
            f"{len(problems)} of {len(objects)} LFS objects are wrong:", file=sys.stderr
        )
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        print(
            "\nTests that read src/krewlyzer/data would skip or fail obscurely.\n"
            "Check the LFS budget, and that the pull or the image fallback ran.",
            file=sys.stderr,
        )
        return 1

    print(f"all {len(objects)} LFS objects present and matching their recorded sha256")
    return 0


if __name__ == "__main__":
    sys.exit(main())
