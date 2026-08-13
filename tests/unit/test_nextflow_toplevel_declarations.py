"""Nextflow 26 forbids top-level statements in a script.

The third thing 26.04 rejected, and the first outside the config:

    subworkflows/local/input_check/main.nf:16:1: Statements cannot be mixed
    with script declarations -- move statements into a process, workflow,
    or function

`def get_index = { bam -> ... }` at column 1 is a *statement*, and 26 wants
only declarations at the top level of a `.nf`. A plain `def name(args) { }` is
a declaration and is accepted by both 25 and 26.

Only top-level assignments are a problem. The same file has three closures
inside its `workflow` block -- `resolvePon`, `resolveTargets`,
`resolveWpsAnchors` -- and those are fine, which is why this checks indentation
rather than banning the construct outright.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_NF = sorted((Path(__file__).resolve().parents[2] / "nextflow").rglob("*.nf"))


@pytest.mark.parametrize("path", _NF, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_no_top_level_closure_assignment(path: Path):
    offenders = [
        (i, line.rstrip())
        for i, line in enumerate(path.read_text().splitlines(), start=1)
        if re.match(r"^(def\s+)?[a-zA-Z_][A-Za-z0-9_]*\s*=\s*\{", line)
    ]
    assert not offenders, (
        f"{path.name} assigns a closure at the top level: "
        + "; ".join(f"line {i}: {t[:48]}" for i, t in offenders)
        + ". Nextflow 26 will not compile the script -- 'Statements cannot be "
        "mixed with script declarations'. Declare it as a function instead, "
        "`def name(args) { ... }`, which both 25 and 26 accept. Closures "
        "*inside* a process or workflow block are unaffected."
    )
