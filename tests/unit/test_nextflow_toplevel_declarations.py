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


def _uses_task_var_outside_closure(line: str) -> bool:
    """Is `meta`/`task` referenced outside every closure on this line?

    A scanner, not a regex, because the two brace kinds nest: a closure may
    contain `${...}` and a `${...}` is itself braces. Two regexes tried before
    this each failed in one direction -- one deleted `${meta.id}` as though it
    were a closure and passed everything, the other could not span a closure
    containing an interpolation and failed the correct code.

    `${` opens an interpolation; a bare `{` opens a closure. A reference is
    safe only while at least one closure is open, because that defers it to
    task time.
    """
    # A STACK, not a counter. Counting only closure-opens while decrementing on
    # every close lets a `${...}` inside a closure pop the closure -- which made
    # the corrected line look like an offender. Each brace must remember which
    # kind it was.
    stack: list[bool] = []  # True = closure, False = ${...} interpolation
    i = 0
    while i < len(line):
        if line[i] == "{":
            stack.append(i == 0 or line[i - 1] != "$")
        elif line[i] == "}":
            if stack:
                stack.pop()
        elif not any(stack) and line.startswith(("meta", "task"), i):
            before = line[i - 1] if i else " "
            after = line[i + 4] if i + 4 < len(line) else " "
            if not before.isalnum() and not after.isalnum() and after != "_":
                return True
        i += 1
    return False


@pytest.mark.parametrize("path", _NF, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_directives_defer_task_variables_to_a_closure(path: Path):
    """A directive string is evaluated when the process is DEFINED.

    `publishDir "${params.outdir}/${meta.id}"` reads `meta` before any task
    exists. 26 refuses the module -- "No such variable: meta" -- while 25
    tolerated it. `path: { ... }` defers evaluation to task time and both
    accept it.

    `tag` is exempt: Nextflow evaluates it lazily by design, and all fourteen
    modules here use `tag "$meta.id"` successfully under both versions.
    """
    offenders = []
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped.startswith(("publishDir", "storeDir")):
            continue
        if _uses_task_var_outside_closure(stripped):
            offenders.append((i, stripped[:60]))

    assert not offenders, (
        f"{path.name} interpolates a task variable in a directive string: "
        + "; ".join(f"line {i}: {t}" for i, t in offenders)
        + ". It is evaluated at process-definition time, before any task "
        "exists, so 26 fails to load the module. Use `path: { ... }`."
    )
