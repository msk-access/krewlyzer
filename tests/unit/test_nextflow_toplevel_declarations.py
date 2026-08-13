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
def test_a_closure_that_is_called_should_be_a_function(path: Path):
    """One rule covering both ways 26 rejected this.

    My first version of this test checked indentation, on the belief that only
    *top-level* closures were a problem and workflow-scoped ones were fine.
    They are not: 26 parses a workflow-scoped closure but cannot resolve it
    from inside a nested `.map { }`, giving "`resolvePon` is not defined".
    Placement was not the only thing tightened; scoping was too.

    So the rule is about usage, not position. If a closure is *invoked* like a
    function anywhere in the file, it should be declared as one -- that works
    at any scope in both versions. A closure genuinely used as a value, passed
    to an operator, never matches `name(` and is left alone.
    """
    text = path.read_text()
    offenders = []
    for i, line in enumerate(text.splitlines(), start=1):
        m = re.match(r"\s*(?:def\s+)?([a-zA-Z_][A-Za-z0-9_]*)\s*=\s*\{", line)
        if m and re.search(rf"\b{re.escape(m.group(1))}\s*\(", text):
            offenders.append((i, m.group(1)))

    assert not offenders, (
        f"{path.name} assigns closures that are then called as functions: "
        + "; ".join(f"line {i}: {n}" for i, n in offenders)
        + ". Nextflow 26 rejects this -- at the top level as 'Statements cannot "
        "be mixed with script declarations', and at workflow scope as "
        "'`name` is not defined' when called from a nested closure. Declare "
        "them as `def name(args) { ... }`, which resolves at any scope in both "
        "25 and 26."
    )
