#!/usr/bin/env python3
"""Regenerate ``krewlyzer_all_docs.md`` from ``docs/``.

Why this is generated
---------------------
The file is a single-document concatenation of the whole ``docs/`` tree, used
where one file is easier to hand around than a site -- an LLM context dump, an
offline read, a grep target.

It was maintained by hand, and had drifted exactly the way this repo's
documentation always drifts: four documents missing (including
``docs/reference/output-files.md``, which describes every output table), two
listed that had been deleted, and content predating several output changes.

A hand-maintained copy of another file is a second source of truth, and the
whole point of `.agents/rules/output-contract.md` is that this repo cannot
afford those. Generating it makes staleness impossible; a test makes it
detectable if someone forgets to run this.

Usage
-----
    python scripts/build_all_docs.py            # rewrite the file
    python scripts/build_all_docs.py --check    # exit 1 if it is stale
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Optional

_ROOT = Path(__file__).resolve().parents[1]
DOCS_DIR = _ROOT / "docs"
OUTPUT = _ROOT / "krewlyzer_all_docs.md"

#: Written at the top so nobody edits the file by hand. Anything above the
#: first ``---`` separator is this banner.
BANNER = """<!--
GENERATED FILE -- DO NOT EDIT.

Concatenation of every Markdown file under docs/, produced by
scripts/build_all_docs.py. Edit the source document instead, then run:

    python scripts/build_all_docs.py

tests/unit/test_all_docs_generated.py fails if this file is out of date.
-->
"""


def source_files(docs_dir: Optional[Path] = None) -> List[Path]:
    """Every Markdown file under ``docs/``, in a stable order.

    ``docs_dir`` defaults to :data:`DOCS_DIR` at *call* time, not at definition
    time. Binding it as a default argument would freeze the module-level value
    at import, so overriding ``DOCS_DIR`` -- which the tests do to check the
    empty-tree guard -- would have no effect and the guard would look broken.

    Sorted by POSIX path so the output is reproducible across filesystems --
    an unsorted walk would make the file churn between machines and turn every
    regeneration into an unreviewable diff.
    """
    root = DOCS_DIR if docs_dir is None else docs_dir
    return sorted(root.rglob("*.md"), key=lambda p: p.relative_to(_ROOT).as_posix())


def render(paths: Iterable[Path]) -> str:
    """Build the concatenated document.

    Section format is unchanged from the hand-written file, so this is a
    content update rather than a format change:

        ---

        # FILE: docs/whatever.md

        <contents>
    """
    parts = [BANNER]
    for path in paths:
        rel = path.relative_to(_ROOT).as_posix()
        body = path.read_text(encoding="utf-8").rstrip("\n")
        parts.append(f"\n---\n\n# FILE: {rel}\n\n{body}\n")
    return "".join(parts)


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--check",
        action="store_true",
        help="Do not write; exit 1 if the file is missing or out of date",
    )
    args = ap.parse_args(argv)

    if not DOCS_DIR.is_dir():
        print(f"error: {DOCS_DIR} does not exist", file=sys.stderr)
        return 2

    paths = source_files()
    if not paths:
        # An empty docs/ would otherwise silently truncate the file to a banner.
        print(f"error: no Markdown files under {DOCS_DIR}", file=sys.stderr)
        return 2

    rendered = render(paths)

    if args.check:
        if not OUTPUT.exists():
            print(f"{OUTPUT.name} is missing; run scripts/build_all_docs.py")
            return 1
        if OUTPUT.read_text(encoding="utf-8") != rendered:
            print(
                f"{OUTPUT.name} is out of date; run scripts/build_all_docs.py",
                file=sys.stderr,
            )
            return 1
        print(f"{OUTPUT.name} is up to date ({len(paths)} documents)")
        return 0

    OUTPUT.write_text(rendered, encoding="utf-8")
    print(
        f"wrote {OUTPUT.name}: {len(paths)} documents, {len(rendered.splitlines())} lines"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
