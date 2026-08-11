"""OCF must survive `--output-format`, and it did not.

Under `--output-format parquet` **every OCF table was silently discarded** —
all six of them: `OCF`, `OCF.sync`, `OCF.{on,off}target`,
`OCF.{on,off}target.sync`. Measured on a real XS1 plasma BAM: 44 logical tables
under `both`, 38 under `parquet`, and exactly the OCF six missing.

The temp directory OCF writes into is an **internal intermediate**: Python
moves the files to their final names and converts them, because
`_core.ocf.apply_pon_zscore` is a line-oriented TSV reader and
`_write_ocf_output` produces the requested format from its result. But
`pipeline.rs` dispatched that intermediate write on the user's
`--output-format`, so under parquet it wrote `all.ocf.parquet` while
`unified_processor.py` looked for `all.ocf.tsv`. Nothing matched, nothing
moved, and the temp directory was deleted. An absent file is indistinguishable
from "OCF was not requested", so there was no warning.

The FSC counts three lines above in `pipeline.rs` already carried the right
rule — *"an internal intermediate; keep as TSV"* — and OCF did not follow it.

This is the third instance of one defect: a hardcoded `.tsv` filename tested
with `.exists()` against a writer that honours `--output-format`. The first was
FSD (see `test_fsd_pon_format_parity.py`), the second FSD's on-target branch,
fifteen lines below the fix for the first.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_ROOT = Path(__file__).resolve().parents[2]


def test_the_ocf_intermediate_is_written_as_tsv():
    """`pipeline.rs` must not dispatch the temp write on the user's format.

    Source-level because the alternative is a full pipeline run, and this pins
    the one line whose change lost six tables. The Python mover on the other
    side of this contract reads those exact `.tsv` names.
    """
    source = (_ROOT / "rust/src/pipeline.rs").read_text()
    ocf_block = source[source.index("if let Some(c) = &self.ocf {") :]
    ocf_block = ocf_block[: ocf_block.index("Ok(())")]

    assert 'write_output_format(&p, "tsv", false)' in ocf_block, (
        "the OCF intermediate is no longer forced to TSV. If it dispatches on "
        "`output_format` again, a parquet run writes all.ocf.parquet, the "
        "mover in unified_processor.py finds no all.ocf.tsv, and all six OCF "
        "tables are silently dropped."
    )
    assert "write_output_format(&p, output_format, compress)" not in ocf_block


def test_the_mover_reports_an_empty_intermediate_rather_than_shrugging():
    """A missing intermediate must be loud.

    Each `if rust_*.exists():` treats absence as "this split was not
    produced" — true for on/off-target in WGS mode, catastrophic for the
    primary table. Without a check up front the two are indistinguishable,
    which is why the loss went unnoticed through a whole release.
    """
    source = (_ROOT / "src/krewlyzer/core/unified_processor.py").read_text()
    assert (
        'ocf_tmp_dir.glob("all.ocf*")' in source
    ), "the guard that reports an empty OCF temp directory is gone"
    marker = source.index('ocf_tmp_dir.glob("all.ocf*")')
    assert (
        "logger.error" in source[marker - 400 : marker + 400]
    ), "an OCF run that produced nothing must be an error, not silence"


def test_the_fsd_ontarget_branch_resolves_the_written_file():
    """The same defect, one branch below its own fix.

    `unified_processor.py` guarded the on-target FSD table with `.exists()` on
    a hardcoded `.FSD.ontarget.tsv`, fifteen lines below the comment explaining
    why the genome-wide branch had stopped doing exactly that. Under parquet
    the on-target table shipped unscored — 69 columns against the genome-wide
    table's 137, and zero `_logR`.
    """
    source = (_ROOT / "src/krewlyzer/core/unified_processor.py").read_text()
    block = source[source.index("# On-target FSD: use fsd_baseline_ontarget") :][:1400]

    assert (
        "resolve_table_path(out_fsd_on)" in block
    ), "the on-target FSD branch is testing a hardcoded .tsv again"
    assert not re.search(r"\bout_fsd_on\.exists\(\)", block), (
        "`.exists()` on the .tsv is back; under --output-format parquet that "
        "file never exists and on-target FSD goes unscored"
    )


def test_python_owns_the_ocf_format_conversion():
    """Why the intermediate can be TSV: the whole chain downstream is.

    If `apply_ocf_python_pon` or `_write_ocf_output` ever stop converting, the
    TSV-only intermediate becomes the shipped output and the reasoning in
    `pipeline.rs` no longer holds.
    """
    from krewlyzer.core import ocf_processor

    source = Path(ocf_processor.__file__).read_text()
    for fn in ("_write_ocf_output", "apply_ocf_python_pon"):
        body = source[source.index(f"def {fn}(") :]
        body = body[: body.index("\ndef ", 1)] if "\ndef " in body[1:] else body
        assert "write_table(" in body, f"{fn} no longer writes the final format"
        assert "cleanup_intermediate_tsv(" in body, (
            f"{fn} leaves its intermediate .tsv behind. A parquet run then "
            "ships two files for one table with no way to tell which is "
            "current -- the on-target path did exactly this."
        )
