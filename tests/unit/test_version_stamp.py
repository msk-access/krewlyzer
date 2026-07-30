"""The version krewlyzer stamps into output must be the version it is.

``FeatureSerializer`` defaulted to a hardcoded ``"0.8.3"`` in two signatures
and ``wrapper.py`` passed a third copy at the call site. The release guide did
have ``sed`` commands for all three, so this was not shipping wrong today --
it was four hand-written regexes over five files, driven by an ``OLD_VERSION``
variable the releaser types in, with nothing checking afterwards that they all
matched. A mistyped ``OLD_VERSION``, or a release cut by hand, silently leaves
some copies behind, and the result is a plausible version string that passes
every check there is.

Now there is one definition and the copies import it. These tests pin that:
the stamp equals the package version, no module restates a version literal,
and the files the release guide *does* edit agree with each other.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import List, Tuple

import pytest

import krewlyzer
from krewlyzer.core.feature_serializer import FeatureSerializer

#: Anchored on this file, not on ``krewlyzer.__file__``.
#:
#: An editable install puts the package inside the checkout, so
#: ``Path(krewlyzer.__file__).parents[1]`` is the repo root locally -- and is
#: ``.venv/lib/python3.10`` in CI, where the package is installed properly.
#: The first version of this test used that and passed locally while failing
#: on the runner.
_REPO = Path(__file__).resolve().parents[2]
_SRC = _REPO / "src" / "krewlyzer"

#: Where a version literal is legitimate: the single source of truth itself.
_ALLOWED = {"__init__.py"}

_VERSION_LITERAL = re.compile(r"^\d+\.\d+\.\d+([-.a-z0-9]*)$")


def test_the_serializer_default_is_the_package_version():
    assert FeatureSerializer("S1").version == krewlyzer.__version__


def test_from_outputs_stamps_the_package_version(tmp_path):
    """The classmethod carries its own default, so it needs its own check."""
    serializer = FeatureSerializer.from_outputs("S1", tmp_path)
    assert serializer.version == krewlyzer.__version__


def test_the_stamp_reaches_the_serialized_payload():
    """A correct attribute that never gets written would still be a bug."""
    serializer = FeatureSerializer("S1")
    payload = serializer.to_dict()
    assert payload["krewlyzer_version"] == krewlyzer.__version__


def _hardcoded_version_literals() -> List[Tuple[str, int, str]]:
    """Every version-shaped string literal in the package, with its location.

    An AST walk rather than a grep: a version in a docstring or a comment about
    historical behaviour ("files written by krewlyzer <= 0.8.3") is fine and
    common in this codebase, while one in executable code is the defect.
    """
    found = []
    for path in sorted(_SRC.rglob("*.py")):
        if path.name in _ALLOWED:
            continue
        tree = ast.parse(path.read_text(), filename=str(path))
        # Drop docstrings; they are the documented-history case, not code.
        docstrings = {
            id(node.body[0].value)
            for node in ast.walk(tree)
            if isinstance(
                node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
            )
            and node.body
            and isinstance(node.body[0], ast.Expr)
            and isinstance(node.body[0].value, ast.Constant)
            and isinstance(node.body[0].value.value, str)
        }
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and id(node) not in docstrings
                and _VERSION_LITERAL.match(node.value)
            ):
                found.append(
                    (
                        str(path.relative_to(_SRC)),
                        node.lineno,
                        node.value,
                    )
                )
    return found


def test_no_module_restates_the_version():
    """One source of truth, enforced.

    A second copy is not wrong on the day it is written -- it is wrong one
    release later, silently, in output that has already shipped.
    """
    offenders = _hardcoded_version_literals()
    assert not offenders, (
        "version literal(s) outside krewlyzer/__init__.py:\n"
        + "\n".join(f"  {f}:{line} -> {value!r}" for f, line, value in offenders)
        + "\nImport __version__ instead; the release bump only edits "
        "__init__.py, so a copy here falls a release behind."
    )


@pytest.mark.parametrize(
    "path,pattern",
    [
        ("src/krewlyzer/__init__.py", r'__version__ = "([^"]+)"'),
        ("pyproject.toml", r'^version = "([^"]+)"'),
        ("rust/Cargo.toml", r'^version = "([^"]+)"'),
    ],
)
def test_the_declared_versions_agree(path, pattern):
    """The files the release guide's sed list actually touches.

    Catches a partial bump -- the guide runs four separate sed commands and
    nothing verifies they all landed.
    """
    text = (_REPO / path).read_text()
    match = re.search(pattern, text, re.MULTILINE)
    assert match, f"no version found in {path} with {pattern!r}"
    assert match.group(1) == krewlyzer.__version__, (
        f"{path} declares {match.group(1)!r} but the package reports "
        f"{krewlyzer.__version__!r}; a version bump was applied unevenly"
    )
