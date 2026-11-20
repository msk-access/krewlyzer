"""Only the resolver may spell a GC factor-table filename (redesign item 0).

Three call sites each built `{sample}.correction_factors[.ontarget].tsv` by hand and tested
whether it existed; production writes `.tsv.gz` and `.parquet`, so they found nothing and ran
uncorrected. `core.gc_assets` (`factor_table_base`, `factor_table_path`, `resolve_gc_factors`)
is now the one place a factor filename is built. This fails if any other module builds one
again. Docstrings and comments are ignored; the allowlist names modules that only *mention* the
suffix for reading QC output, never to choose a table for correction.
"""

import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[2] / "src" / "krewlyzer"
ALLOWED = {
    "core/gc_assets.py": "the resolver itself",
    "core/feature_serializer.py": "features.json only; resolves any format via _resolve_output_path",
    "validate/plots.py": "QC plot of the Parquet table run-all wrote",
    "validate/meaning.py": "documentation of an output suffix",
    "validate/contract.py": "output-contract suffix list",
}


def _docstring_nodes(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(
            node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ):
            body = getattr(node, "body", [])
            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
            ):
                out.add(id(body[0].value))
    return out


def _offenders():
    found = []
    for path in sorted(PACKAGE.rglob("*.py")):
        rel = path.relative_to(PACKAGE).as_posix()
        if rel in ALLOWED:
            continue
        tree = ast.parse(path.read_text())
        docs = _docstring_nodes(tree)
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and id(node) not in docs
            ):
                if ".correction_factors" in node.value:
                    found.append(f"{rel}:{node.lineno}: {node.value!r}")
    return found


def test_no_module_but_the_resolver_builds_a_factor_filename():
    offenders = _offenders()
    assert not offenders, (
        "Build GC factor paths with core.gc_assets (factor_table_path / resolve_gc_factors):\n"
        + "\n".join(offenders)
    )


def test_the_check_would_catch_the_old_lookup(tmp_path):
    # The shape of the defect, in a scratch package: the check must see it.
    pkg = tmp_path / "krewlyzer" / "core"
    pkg.mkdir(parents=True)
    (pkg / "unified_processor.py").write_text(
        "def f(output_dir, sample_name):\n"
        '    return output_dir / f"{sample_name}.correction_factors.tsv"\n'
    )
    global PACKAGE
    saved, PACKAGE = PACKAGE, tmp_path / "krewlyzer"
    try:
        assert _offenders() == [
            "core/unified_processor.py:2: '.correction_factors.tsv'"
        ]
    finally:
        PACKAGE = saved
