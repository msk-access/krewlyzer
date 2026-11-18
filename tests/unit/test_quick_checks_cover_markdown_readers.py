"""Every test that reads markdown must run in quick-checks.yml.

test.yml skips the expensive suite when a change touches only documentation.
That is safe only because quick-checks.yml, which has no path filter, runs every
test that could be broken by such a change. The list of those tests used to be a
comment saying "add yours here" -- a rule enforced by remembering it. This makes
it a test.

"Reads markdown" means *constructs a path to it*, not mentions it. Found by
walking the AST, so docstrings and failure messages do not count: one invariant
test names `fsc.md` in its docstring and never opens it, and treating that as a
read would demand a maturin build in a job that exists to avoid one.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_ROOT = Path(__file__).resolve().parents[2]
_TESTS = _ROOT / "tests"
_WORKFLOW = _ROOT / ".github" / "workflows" / "quick-checks.yml"
_SELF = Path(__file__).resolve()

#: Readers deliberately not run by quick-checks.yml, and what covers them instead.
#: An entry here is a claim that something else runs this test on every change
#: that could break it; write down what.
_COVERED_ELSEWHERE = {
    "tests/unit/test_all_docs_generated.py": (
        "docs.yml runs `build_all_docs.py --check`, and its `paths` cover every "
        "input this test has: docs/**, krewlyzer_all_docs.md, the script."
    ),
}

_CALLS = {"Path", "PurePath", "open", "glob", "rglob", "read_text", "joinpath"}


def _names_markdown(value: str) -> bool:
    return value == "docs" or value.endswith(".md") or "*.md" in value


def _docstring_ids(tree: ast.AST) -> set[int]:
    ids = set()
    for node in ast.walk(tree):
        if isinstance(
            node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            body = getattr(node, "body", [])
            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
            ):
                ids.add(id(body[0].value))
    return ids


def _reads_markdown(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    skip = _docstring_ids(tree)

    def is_md(node: ast.AST) -> bool:
        return (
            isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and id(node) not in skip
            and _names_markdown(node.value)
        )

    for node in ast.walk(tree):
        # `root / "docs" / "x.md"`
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            if is_md(node.left) or is_md(node.right):
                return True
        # `Path("README.md")`, `docs.rglob("*.md")`, `open("x.md")`
        if isinstance(node, ast.Call):
            fn = node.func
            name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
            if name in _CALLS and any(is_md(a) for a in node.args):
                return True
    return False


def _readers() -> list[str]:
    return sorted(
        str(p.relative_to(_ROOT))
        for p in _TESTS.rglob("*.py")
        if p.resolve() != _SELF and _reads_markdown(p)
    )


def _quick_checks_tests() -> set[str]:
    """Test files the workflow actually runs -- whole-line comments excluded.

    A positive check, so comments must not count: naming a file in a comment is
    not running it. The first version of this guard counted tests/conftest.py as
    listed because the workflow's comment about pysam mentions it -- the same
    mistake as test_lint_script_matches_ci once made with `--tests`.
    """
    lines = _WORKFLOW.read_text(encoding="utf-8").splitlines()
    code = "\n".join(line for line in lines if not line.lstrip().startswith("#"))
    return set(re.findall(r"tests/[\w/]+\.py", code))


def test_the_sweep_still_finds_the_known_readers():
    """Guards the detector, not the workflow.

    If a refactor of these tests changed how they build paths, the sweep could
    quietly find nothing -- and then the check below would pass vacuously.
    """
    found = set(_readers())
    for known in (
        "tests/unit/test_cli_documented.py",
        "tests/unit/test_nextflow_params.py",
        "tests/unit/test_all_docs_generated.py",
    ):
        assert known in found, (
            f"the markdown-reader sweep no longer detects {known}; the detector "
            "has drifted from how tests build paths, so it cannot be trusted"
        )


def test_every_markdown_reader_runs_on_markdown_only_changes():
    listed = _quick_checks_tests()
    missing = [r for r in _readers() if r not in listed and r not in _COVERED_ELSEWHERE]
    assert not missing, (
        f"{missing} read markdown but are not run by quick-checks.yml. test.yml "
        "skips markdown-only changes, so nothing would run them on the change "
        "most likely to break them. Add them to the `Run the text-only tests` "
        "step, or to _COVERED_ELSEWHERE with what covers them instead."
    )


def test_exceptions_still_describe_real_readers():
    """A stale exception is a hole someone forgot to close."""
    readers = set(_readers())
    stale = [e for e in _COVERED_ELSEWHERE if e not in readers]
    assert (
        not stale
    ), f"_COVERED_ELSEWHERE lists files that no longer read markdown: {stale}"


def test_every_listed_test_exists():
    """A renamed test file would otherwise drop out of the job silently --
    pytest given a missing path errors, but only if the step ever runs."""
    missing = sorted(t for t in _quick_checks_tests() if not (_ROOT / t).is_file())
    assert not missing, f"quick-checks.yml runs tests that do not exist: {missing}"
