"""A skip flag must stop the asset resolver, not just add a CLI flag.

`--skip_pon` and `--skip_target_regions` reach the modules as `--skip-pon` and
`--skip-target-regions`. The PON and targets reach them separately, through
`resolvePon` / `resolveTargets` in the input_check subworkflow, which return
the samplesheet's file, the assay's bundled file, or `params.pon_model` /
`params.targets`. Those resolvers never looked at the skip flags, so a panel
sample with `--skip_pon` got both `--pon-model X` and `--skip-pon`, and
`asset_resolution.resolve_pon_model` raises ValueError on that pair (likewise
for targets).

The task fails, and under a profile that retries and then ignores errors the
run still reports complete, with those tables missing. That is a silent
failure, so this is checked here.

The fix is at the source: each resolver returns `[]` when its skip flag is set,
before it consults any other source. Every module fed by a resolver then builds
no `--pon-model` / `--target-regions` argument, because `[]` is falsy.

Pure text analysis: Nextflow and Groovy are not installed locally or in CI.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[2]
_INPUT_CHECK = _ROOT / "nextflow" / "subworkflows" / "local" / "input_check" / "main.nf"

pytestmark = [
    pytest.mark.unit,
    pytest.mark.skipif(not _INPUT_CHECK.is_file(), reason="nextflow/ not present"),
]


def _function_body(text: str, name: str) -> str:
    """The text between the braces of a top-level `def name(...) { ... }`.

    Brace-matched rather than regex-matched: the bodies hold map literals and
    GString `${...}` interpolations, so a lazy `\\{.*?\\}` stops early.
    """
    m = re.search(rf"^def\s+{re.escape(name)}\s*\([^)]*\)\s*\{{", text, re.MULTILINE)
    assert m, f"{name}() is not declared as a top-level function in {_INPUT_CHECK.name}"
    depth, start = 1, m.end()
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start:i]
    raise AssertionError(f"unbalanced braces in {name}()")


@pytest.mark.parametrize(
    "resolver, flag",
    [("resolvePon", "skip_pon"), ("resolveTargets", "skip_target_regions")],
)
def test_skip_flag_returns_empty_before_any_other_source(resolver: str, flag: str):
    """The skip guard must be the resolver's first return.

    A guard placed after `if (row.pon) return file(row.pon)` would still hand a
    samplesheet PON to a run that asked to skip it, so presence alone is not
    enough: the guard must come before every other `return`.
    """
    body = _function_body(_INPUT_CHECK.read_text(), resolver)

    guard = re.search(rf"if\s*\(\s*params\.{flag}\s*\)\s*return\s*\[\s*\]", body)
    assert guard, (
        f"{resolver}() never returns [] for params.{flag}. The module then gets "
        f"the resolved file *and* the skip flag, and asset_resolution raises "
        f"ValueError on that pair."
    )

    first_return = re.search(r"\breturn\b", body)
    assert first_return and first_return.start() >= guard.start(), (
        f"{resolver}() can return a file before it checks params.{flag}: "
        f"{body[first_return.start():].splitlines()[0].strip()!r} comes first."
    )
