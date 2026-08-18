"""Every optional import must be installable, and `[all]` must mean all.

This exists because of a bug it would have caught. `validate/describe.py` did
`import markdown` behind a `try/except ImportError`, and `markdown` appeared in
no dependency group at all -- so no install of krewlyzer could satisfy it, and
`describe-output -o page.html` wrote a valid HTML document whose body was raw
Markdown in a `<pre>`. It worked only where an unrelated package (mkdocs) had
dragged the dependency in, which is precisely the condition under which nobody
notices.

The shape of that defect is general: an optional import is a promise that the
feature is available *if you install the right thing*, and nothing checked that
the right thing existed. `psutil` was in the same state and turned up when this
test was written.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

try:  # 3.11+
    import tomllib
except ModuleNotFoundError:  # 3.10, which is what CI runs
    import tomli as tomllib  # type: ignore[no-redef]

_ROOT = Path(__file__).resolve().parents[2]
_SRC = _ROOT / "src" / "krewlyzer"

#: Optional imports that are deliberately not pip-installable.
#:
#: `krewlyzer._core` is the compiled pyo3 extension, built by maturin and never
#: resolved from an index. Add to this only with a reason: the default answer
#: for an optional import is that it belongs in an extra.
_NOT_ON_PYPI = {"krewlyzer"}

#: Import name -> distribution name, where they differ. None needed yet; the
#: three current optional imports all match. Kept so the next one that does not
#: has an obvious home rather than a special case bolted into the assertion.
_DISTRIBUTION_NAME: dict[str, str] = {}

#: Groups that exist for contributors, not for users of the tool. `[all]` is
#: not expected to cover these -- see the comment above `all` in pyproject.
_CONTRIBUTOR_EXTRAS = {"docs", "test", "dev"}


def _extras() -> dict[str, list[str]]:
    with (_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["project"]["optional-dependencies"]


def _requirement_name(spec: str) -> str:
    """`plotly>=5.18.0` -> `plotly`; `krewlyzer[report]` -> `krewlyzer`."""
    return re.split(r"[\[<>=!;~\s]", spec, maxsplit=1)[0].strip().lower()


def _optional_imports() -> dict[str, Path]:
    """Top-level modules imported inside a `try:` guarded by `except ImportError`.

    Walks the AST rather than grepping: an `import` inside a try block is the
    structure that matters, and a regex over source cannot tell that from an
    import that merely appears near the word ImportError.
    """
    found: dict[str, Path] = {}
    for path in sorted(_SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Try):
                continue
            catches_import_error = any(
                (isinstance(h.type, ast.Name) and h.type.id == "ImportError")
                or (
                    isinstance(h.type, ast.Tuple)
                    and any(
                        isinstance(e, ast.Name) and e.id == "ImportError"
                        for e in h.type.elts
                    )
                )
                for h in node.handlers
            )
            if not catches_import_error:
                continue
            for inner in ast.walk(node):
                if isinstance(inner, ast.Import):
                    for alias in inner.names:
                        found.setdefault(alias.name.split(".")[0], path)
                elif (
                    isinstance(inner, ast.ImportFrom)
                    and inner.module
                    and not inner.level
                ):
                    found.setdefault(inner.module.split(".")[0], path)
    return found


def test_every_optional_import_is_declared_in_some_extra():
    """An import nothing can install is a feature nobody can turn on."""
    declared = {
        _requirement_name(spec) for specs in _extras().values() for spec in specs
    }
    undeclared = {
        module: path
        for module, path in _optional_imports().items()
        if module not in _NOT_ON_PYPI
        and _DISTRIBUTION_NAME.get(module, module).lower() not in declared
    }
    assert not undeclared, "optional imports that no extra provides: " + ", ".join(
        f"{m} ({p.relative_to(_ROOT)})" for m, p in sorted(undeclared.items())
    )


def test_all_covers_every_runtime_extra():
    """`[all]` has to keep meaning "all", including extras added later.

    It self-references (`krewlyzer[report]`) so the contents cannot drift, but
    a *new* runtime extra would not add itself. This is the part that has to be
    checked rather than arranged.
    """
    extras = _extras()
    referenced: set[str] = set()
    for spec in extras["all"]:
        match = re.fullmatch(r"krewlyzer\[([a-z0-9_,-]+)\]", spec.strip())
        if match:
            referenced.update(match.group(1).split(","))

    runtime = set(extras) - _CONTRIBUTOR_EXTRAS - {"all"}
    missing = runtime - referenced
    assert not missing, (
        f"[all] does not include runtime extra(s): {sorted(missing)}. "
        "Add krewlyzer[<name>] to `all` in pyproject.toml, or move the group "
        "into _CONTRIBUTOR_EXTRAS here if it is tooling rather than a feature."
    )


def test_all_does_not_drag_in_contributor_tooling():
    """`pip install 'krewlyzer[all]'` must not deliver mypy and mkdocs.

    The counter-pressure to the test above: "all" is easy to widen until a
    pipeline image carries a linter.
    """
    declared = {_requirement_name(s) for s in _extras()["all"]}
    tooling = {
        _requirement_name(s)
        for name in _CONTRIBUTOR_EXTRAS
        for s in _extras()[name]
        # plotly is in `test` only because the report path is tested there.
        if _requirement_name(s) not in {"plotly", "tomli"}
    }
    assert not (declared & tooling), sorted(declared & tooling)


@pytest.mark.parametrize("module", sorted(_optional_imports()))
def test_optional_imports_are_never_module_level(module):
    """An optional module must not also be imported at module scope.

    One unconditional top-level import turns the whole install into a hard
    dependency and makes the careful fallback decoration -- the package fails
    at import time, long before the code path that was supposed to degrade.

    Only module scope counts. An import nested inside a function, a `try`, or
    an `if` is deferred and conditional by construction: `render_html_page`
    reaches its `import markdown` only after `markdown_renderer_available()`
    says so, which is as good a guard as a try/except and reads better.
    """
    if module in _NOT_ON_PYPI:
        pytest.skip("compiled extension, not resolved from an index")
    offenders = []
    for path in sorted(_SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in tree.body:  # module scope only
            names = []
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                names = [node.module.split(".")[0]]
            if module in names:
                offenders.append(f"{path.relative_to(_ROOT)}:{node.lineno}")
    assert not offenders, (
        f"{module} is optional elsewhere but imported at module scope in "
        f"{offenders} -- that makes it a hard dependency"
    )
