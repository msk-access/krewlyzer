"""After writing a file, read *that* file — not whatever resolves.

`read_table` is deliberately parquet-first: given `x.FSD.tsv` it prefers
`x.FSD.parquet`. That is correct when asking "what was produced for this
output", and wrong immediately after writing a specific file.

The Rust backends write plain TSV and Python reads it back to honour
`--output-format`. Where a `.parquet` sibling from an earlier run sat in the
directory, `read_table` preferred it and the freshly computed values were
discarded — silently, with every column keeping its name.

It has cost real data twice:

- **FSD** lost its log-ratios this way (`c92ed86`): computed, logged as
  "41 arms normalized", then overwritten by the pre-normalisation table.
- **Region entropy** would have emitted a *previous* run's z-scores on any
  re-run into the same directory — which is exactly what a Nextflow retry or
  `-resume` produces.

`read_exact_table` is the counterpart, and this module pins both the helper's
behaviour and the structural rule that no read-after-write uses the resolving
reader.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.core.output_utils import read_exact_table, read_table

_SRC = Path(__file__).resolve().parents[2] / "src" / "krewlyzer"

pytestmark = pytest.mark.invariant


# ---------------------------------------------------------------------------
# the helper
# ---------------------------------------------------------------------------


def test_the_resolving_reader_prefers_a_stale_parquet(tmp_path):
    """Not a defect in read_table -- the documented behaviour this exists for.

    Pinned so the two readers cannot quietly converge and make
    read_exact_table redundant-looking.
    """
    pd.DataFrame({"v": [1]}).to_csv(tmp_path / "t.tsv", sep="\t", index=False)
    pd.DataFrame({"v": [999]}).to_parquet(tmp_path / "t.parquet")
    assert read_table(tmp_path / "t.tsv")["v"].iloc[0] == 999


def test_the_exact_reader_ignores_the_sibling(tmp_path):
    pd.DataFrame({"v": [1]}).to_csv(tmp_path / "t.tsv", sep="\t", index=False)
    pd.DataFrame({"v": [999]}).to_parquet(tmp_path / "t.parquet")
    assert read_exact_table(tmp_path / "t.tsv")["v"].iloc[0] == 1
    assert read_exact_table(tmp_path / "t.parquet")["v"].iloc[0] == 999


def test_the_exact_reader_returns_none_for_a_missing_file(tmp_path):
    """Matches read_table, so callers need no special case."""
    assert read_exact_table(tmp_path / "nope.tsv") is None


def test_the_exact_reader_skips_comment_footers(tmp_path):
    """Motif output carries a `#` metadata footer (bf4404e)."""
    (tmp_path / "t.tsv").write_text("a\tb\n1\t2\n# written by krewlyzer\n")
    frame = read_exact_table(tmp_path / "t.tsv")
    assert len(frame) == 1 and frame["a"].iloc[0] == 1


# ---------------------------------------------------------------------------
# the structural rule
# ---------------------------------------------------------------------------

_WRITERS = {
    "write_table",
    "to_csv",
    "to_parquet",
    "apply_pon_logratio",
    "apply_pon_zscore",
    "run_region_mds",
}
_RESOLVING_READERS = {"read_table", "load_entropy_tsv"}

#: Reads that follow a write but legitimately target a *different*, earlier
#: file rather than the one just written. Each needs a reason.
_ALLOWED = {
    # The degradation branch re-reads the raw entropy input, not the output.
    ("core/region_entropy_processor.py", "process_region_entropy"),
    # A 700-line orchestrator: its reads target other features' outputs
    # entirely, not anything it wrote in the same call.
    ("core/unified_processor.py", "run_features"),
}


def _read_after_write() -> list[tuple[str, str, list[str]]]:
    findings = []
    for path in sorted(_SRC.rglob("*.py")):
        tree = ast.parse(path.read_text())
        for fn in [
            n
            for n in ast.walk(tree)
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        ]:
            seq = sorted(
                (
                    c.lineno,
                    (
                        c.func.attr
                        if isinstance(c.func, ast.Attribute)
                        else getattr(c.func, "id", None)
                    ),
                )
                for c in ast.walk(fn)
                if isinstance(c, ast.Call)
            )
            first_write = next((ln for ln, nm in seq if nm in _WRITERS), None)
            if first_write is None:
                continue
            after = [
                f"{nm}@{ln}"
                for ln, nm in seq
                if nm in _RESOLVING_READERS and ln > first_write
            ]
            if after:
                findings.append((path.relative_to(_SRC).as_posix(), fn.name, after))
    return findings


def test_no_unreviewed_read_after_write():
    """A resolving read after a write is a defect until someone says otherwise.

    Adding a name to `_ALLOWED` requires writing down why the read targets a
    different file than the write. That is the whole guard: the pattern is easy
    to reintroduce and impossible to spot in review.
    """
    unreviewed = [
        (f, fn, after)
        for f, fn, after in _read_after_write()
        if (f, fn) not in _ALLOWED
    ]
    assert not unreviewed, (
        "resolving read after a write, not in the reviewed list:\n"
        + "\n".join(f"  {f}::{fn} -> {after}" for f, fn, after in unreviewed)
        + "\nUse read_exact_table, or add to _ALLOWED with a reason."
    )


def test_the_allowlist_has_no_dead_entries():
    """A stale allowlist entry hides the next real one."""
    live = {(f, fn) for f, fn, _ in _read_after_write()}
    dead = sorted(_ALLOWED - live)
    assert not dead, f"_ALLOWED names sites that no longer read after writing: {dead}"
