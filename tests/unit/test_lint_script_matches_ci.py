"""`scripts/lint.sh` must run what the CI lint job runs.

The script exists so a developer can reproduce CI's lint result before pushing.
That only helps if it is faithful: a script that misses a check, or runs it
with laxer flags, produces confident green locally and a failure in CI --
worse than having no script, because it is trusted.

So the flags are asserted rather than maintained by hand. Both of these have
already happened in this repository:

* mypy 2.3.1 locally against a 1.19.1 pin, which does not report
  `[import-untyped]`, so #71 passed locally and failed CI;
* running clippy without `-A clippy::too_many_arguments` reports 12 errors CI
  does not, which sends you fixing things that are not broken.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

try:  # 3.11+
    import tomllib
except ModuleNotFoundError:  # 3.10, which is what CI runs
    import tomli as tomllib  # type: ignore[no-redef]

_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT = _ROOT / "scripts" / "lint.sh"
_WORKFLOW = _ROOT / ".github" / "workflows" / "test.yml"


def _strip_comment_lines(text: str) -> str:
    """Drop whole-line comments, so a flag named only in prose cannot pass.

    The checks below are substring tests. Without this, documenting a flag is
    indistinguishable from passing it: adding `--tests` to the clippy
    invocation with a comment explaining why left the comment alone sufficient
    to satisfy the check with the flag itself deleted from the command --
    verified by deleting it. A guard that a comment can satisfy is not a guard.

    Only whole-line comments are removed, not trailing ones. None of the
    mirrored tokens contain `#`, so this cannot hide a real flag.
    """
    return "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("#")
    )


def _script(*, strip_comments: bool = True) -> str:
    """The script's text, comments stripped by default.

    Stripping is right for a check that something *must appear* -- a comment
    cannot satisfy a requirement -- and wrong for a check that something *must
    not appear*, because a comment can still violate a prohibition. A version
    restated in a comment drifts exactly like a hardcoded one, for the next
    person who reads it. So negative checks pass `strip_comments=False`.
    """
    text = _SCRIPT.read_text(encoding="utf-8")
    return _strip_comment_lines(text) if strip_comments else text


def _lint_job_commands() -> str:
    """The shell of every `run:` step in the workflow's lint job.

    Read as text on purpose. Parsing the YAML would give the same strings by a
    longer route, and the thing being compared is the shell either way.
    """
    text = _WORKFLOW.read_text(encoding="utf-8")
    start = text.index("\n  lint:")
    return _strip_comment_lines(text[start:])


#: Tokens that must appear in the script for each command CI runs. Paths and
#: flags only -- not ordering, and not the uvx wrapper the script adds.
_MIRRORED = {
    "black": ["--check", "src/krewlyzer/", "tests/"],
    "ruff": ["check", "src/krewlyzer/", "tests/"],
    "mypy": ["src/krewlyzer/", "--ignore-missing-imports", "--no-error-summary"],
    "clippy": [
        "--manifest-path rust/Cargo.toml",
        # Load-bearing: without --tests clippy lints the lib target only and
        # never sees test code, where this repo's invariants are asserted.
        # Listed here so removing it from either the workflow or the script
        # fails this mirror check instead of silently narrowing the gate.
        "--tests",
        "-D warnings",
        "-A clippy::too_many_arguments",
        "-A clippy::type_complexity",
    ],
    "check_output_format": ["scripts/check_output_format.py"],
}


def test_the_lint_script_exists_and_is_executable():
    assert _SCRIPT.is_file(), "scripts/lint.sh is referenced by AGENTS.md"
    assert _SCRIPT.stat().st_mode & 0o111, "scripts/lint.sh must be executable"


@pytest.mark.parametrize("tool", sorted(_MIRRORED))
def test_script_runs_the_same_checks_as_ci(tool):
    """Every flag CI passes must appear in the script."""
    script = _script()
    missing = [token for token in _MIRRORED[tool] if token not in script]
    assert not missing, (
        f"scripts/lint.sh omits {missing} for {tool}, so it does not reproduce "
        "the CI lint job"
    )


@pytest.mark.parametrize("tool", sorted(_MIRRORED))
def test_the_mirrored_flags_are_the_ones_ci_actually_uses(tool):
    """And the list above must describe CI, not a stale memory of it.

    Without this the previous test degenerates: someone relaxes a flag in the
    workflow, `_MIRRORED` keeps the old one, and both tests still pass while
    the script and CI have quietly diverged.
    """
    commands = _lint_job_commands()
    missing = [token for token in _MIRRORED[tool] if token not in commands]
    assert not missing, (
        f"_MIRRORED[{tool!r}] lists {missing}, which the lint job no longer "
        "runs -- update this table and scripts/lint.sh together"
    )


def test_every_pinned_dev_tool_is_run_by_the_script():
    """A pinned tool nobody runs locally is a CI-only surprise waiting to happen."""
    with (_ROOT / "pyproject.toml").open("rb") as handle:
        dev = tomllib.load(handle)["project"]["optional-dependencies"]["dev"]

    pinned = {
        m.group(1)
        for spec in dev
        if (m := re.fullmatch(r"([A-Za-z0-9_.-]+)==[\w.]+", spec.strip()))
    }
    assert pinned, "expected == pins in [dev]"

    script = _script()
    unrun = [tool for tool in sorted(pinned) if f"pin {tool}" not in script]
    assert not unrun, (
        f"pinned in [dev] but never run by scripts/lint.sh: {unrun}. "
        "Either run it there or drop the pin."
    )


def test_the_script_reads_versions_from_pyproject():
    """Versions must not be restated in the script, or they will drift.

    The failure this whole change is about is two sources of truth for a tool
    version disagreeing; a hardcoded pin here would recreate it one level down.
    """
    assert "pyproject.toml" in _script()
    # Raw text, deliberately unlike every positive check in this file: see
    # `_script`. Stripping comments here let `# CI runs black@26.1 here` pass,
    # where the guard before comment-stripping correctly refused it.
    hardcoded = re.findall(r"(black|ruff|mypy)@\d+\.\d+", _script(strip_comments=False))
    assert not hardcoded, (
        f"scripts/lint.sh hardcodes versions {hardcoded}; read them from "
        "pyproject.toml with the `pin` helper instead"
    )
