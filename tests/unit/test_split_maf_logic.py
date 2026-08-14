"""Exercise SPLIT_MAF's embedded Python without Nextflow.

Nothing in CI runs the pipeline -- no stub run, no integration test -- so a
Groovy change is checked by review alone. The splitting logic is the part that
could silently produce wrong output, and it is plain Python, so it can be
lifted out of the module and run directly.

What matters most here is the empty case. A sample with no matching rows that
receives *no file* is dropped by the caller's re-pairing, RUNALL never runs for
it, and the cohort is short with nothing to show for it -- the same silent loss
`validate-output --expect` exists to catch, reintroduced one stage earlier.
"""

from __future__ import annotations

import re
import subprocess
import textwrap
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_MODULE = (
    Path(__file__).resolve().parents[2]
    / "nextflow/modules/local/krewlyzer/split_maf/main.nf"
)
_HEADER = "Hugo_Symbol\tChromosome\tStart_Position\tTumor_Sample_Barcode"


def _script(sample_ids, single, maf_name):
    """Lift the script block out of the module and bind its Groovy templates."""
    text = _MODULE.read_text()
    body = text.split('"""', 1)[1].split('"""', 1)[0]
    # Dedent BEFORE substituting, mirroring Nextflow.
    #
    # NOTE: Nextflow actually does the opposite -- it interpolates first and
    # dedents after. That difference is not cosmetic: it is why this file could
    # not catch the outage. The module used to inject 16,552 id lines at zero
    # indent, which left `stripIndent()` no common prefix, and the task died on
    # `    import os`. Here the substitution happened after the dedent, so the
    # test saw well-formed Python and passed. The ids are a staged file now, so
    # nothing multi-line is interpolated and the orders cannot diverge again --
    # `test_module_interpolations.py` asserts that.
    body = textwrap.dedent(body)
    body = body.replace("${single_flag}", "true" if single else "false")
    body = body.replace("${expected_n}", str(len(sample_ids)))
    body = body.replace("${ids}", "ids.txt")
    body = body.replace("${maf}", maf_name)
    body = body.replace("${task.process}", "SPLIT_MAF")
    # Groovy escapes these for its own interpolation; Python needs them raw.
    return body.replace("\\\\n", "\\n").replace("\\\\t", "\\t")


def _run(tmp_path, sample_ids, rows, single=False):
    maf = tmp_path / "cohort.maf"
    maf.write_text(_HEADER + "\n" + "".join(rows))
    (tmp_path / "ids.txt").write_text("\n".join(sample_ids) + "\n")
    script = tmp_path / "split.py"
    script.write_text(_script(sample_ids, single, "cohort.maf"))
    proc = subprocess.run(
        [sys.executable, "split.py"], cwd=tmp_path, capture_output=True, text=True
    )
    return proc, tmp_path / "split"


def _row(sample, gene="TP53"):
    return f"{gene}\tchr17\t7577120\t{sample}\n"


def test_every_requested_sample_gets_a_file_even_with_no_variants(tmp_path):
    """The load-bearing case: absence must still produce a file."""
    ids = ["sample_01", "sample_02", "sample_03"]
    proc, out = _run(tmp_path, ids, [_row("sample_01"), _row("sample_01")])
    assert proc.returncode == 0, proc.stderr
    produced = sorted(p.name for p in out.iterdir())
    assert produced == [f"{s}.filtered.maf" for s in ids]

    # The two with no variants are header-only, not absent and not empty.
    for sid in ("sample_02", "sample_03"):
        lines = (out / f"{sid}.filtered.maf").read_text().splitlines()
        assert lines == [_HEADER]


def test_rows_are_routed_to_the_right_sample(tmp_path):
    ids = ["sample_01", "sample_02"]
    rows = [_row("sample_01"), _row("sample_02"), _row("sample_02")]
    proc, out = _run(tmp_path, ids, rows)
    assert proc.returncode == 0, proc.stderr
    assert len((out / "sample_01.filtered.maf").read_text().splitlines()) == 2
    assert len((out / "sample_02.filtered.maf").read_text().splitlines()) == 3


def test_barcode_matching_is_case_insensitive(tmp_path):
    """FILTER_MAF matched case-insensitively; the split must not diverge."""
    proc, out = _run(tmp_path, ["sample_01"], [_row("SAMPLE_01")])
    assert proc.returncode == 0, proc.stderr
    assert len((out / "sample_01.filtered.maf").read_text().splitlines()) == 2


def test_a_sample_absent_from_the_maf_does_not_steal_rows(tmp_path):
    proc, out = _run(tmp_path, ["sample_01", "sample_09"], [_row("sample_02")])
    assert proc.returncode == 0, proc.stderr
    for sid in ("sample_01", "sample_09"):
        assert (out / f"{sid}.filtered.maf").read_text().splitlines() == [_HEADER]


def test_single_sample_mode_passes_every_row_to_every_sample(tmp_path):
    """Mode 2 semantics: no filtering, all rows pass through.

    Grouped on `single_sample` by the caller, so a whole group agrees -- which
    is why batching cannot mix pass-through and filtered samples.
    """
    ids = ["sample_01", "sample_02"]
    proc, out = _run(tmp_path, ids, [_row("whatever"), _row("other")], single=True)
    assert proc.returncode == 0, proc.stderr
    for sid in ids:
        assert len((out / f"{sid}.filtered.maf").read_text().splitlines()) == 3


def test_comment_lines_are_stripped(tmp_path):
    maf = tmp_path / "cohort.maf"
    maf.write_text("#version 2.4\n" + _HEADER + "\n" + _row("sample_01"))
    (tmp_path / "ids.txt").write_text("sample_01\n")
    script = tmp_path / "split.py"
    script.write_text(_script(["sample_01"], False, "cohort.maf"))
    proc = subprocess.run(
        [sys.executable, "split.py"], cwd=tmp_path, capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stderr
    text = (tmp_path / "split" / "sample_01.filtered.maf").read_text()
    assert not text.startswith("#")


def test_the_module_refuses_to_emit_a_short_set(tmp_path):
    """The guard that makes the caller's contract enforceable."""
    body = _MODULE.read_text()
    assert "produced != len(sample_ids)" in body, (
        "SPLIT_MAF no longer verifies it wrote one file per requested sample. "
        "Without that, a missing file silently drops the sample downstream."
    )


def test_the_caller_groups_on_single_sample_not_just_the_maf(tmp_path):
    """Batching pass-through and filtered samples together would corrupt both."""
    wf = Path(__file__).resolve().parents[2] / "nextflow/workflows/krewlyzer/main.nf"
    text = wf.read_text()
    assert re.search(r"\[\s*maf\s*,\s*meta\.single_sample", text), (
        "SPLIT_MAF grouping must key on (maf, single_sample). Keyed on the MAF "
        "alone, a pass-through sample and a filtered one sharing a file would "
        "be handled under one mode and one of them would be wrong."
    )
