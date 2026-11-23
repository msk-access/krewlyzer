"""Exercise FILTER_MAF's embedded Python without Nextflow.

Tool-level mode sends every MAF through FILTER_MAF before mFSD. The samplesheet guide says a
MAF may be gzipped, and mFSD's reader recognises gzip from the file's first bytes, but both
of FILTER_MAF's scripts opened the MAF as plain text, so a `.maf.gz` failed the task before
mFSD saw it. Each script is lifted out of the module and run on a plain and a gzipped copy of
the same MAF.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_MODULE = (
    Path(__file__).resolve().parents[2]
    / "nextflow/modules/local/krewlyzer/filter_maf/main.nf"
)
_MAF = (
    "#version 2.4\n"
    "Hugo_Symbol\tChromosome\tStart_Position\tTumor_Sample_Barcode\n"
    "TP53\tchr17\t7577120\tsample_01\n"
    "KRAS\tchr12\t25398284\tsample_02\n"
)


def _script(mode: str) -> str:
    """The Python block for `single` (MODE 2) or `multi` (MODE 1), Groovy templates bound."""
    text = _MODULE.read_text()
    marker = {"single": "// --- MODE 2", "multi": "// --- MODE 1"}[mode]
    body = text.split(marker, 1)[1].split('"""', 1)[1].split('"""', 1)[0]
    body = textwrap.dedent(body)
    body = body.replace("${prefix}", "sample_01").replace("${maf}", "in.maf")
    body = body.replace("${task.process}", "FILTER_MAF")
    return body.replace("\\\\n", "\\n").replace("\\\\t", "\\t")


def _run(folder: Path, mode: str, gz: bool) -> str:
    folder.mkdir()
    maf = folder / "in.maf"  # the name says plain; only the bytes say gzip
    if gz:
        maf.write_bytes(gzip.compress(_MAF.encode()))
    else:
        maf.write_text(_MAF)
    (folder / "filter.py").write_text(_script(mode))
    proc = subprocess.run(
        [sys.executable, "filter.py"], cwd=folder, capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stderr
    return (folder / "sample_01.filtered.maf").read_text()


@pytest.mark.parametrize("mode", ["single", "multi"])
def test_a_gzipped_maf_filters_as_the_plain_one(tmp_path, mode):
    plain = _run(tmp_path / "plain", mode, gz=False)
    assert (
        "TP53" in plain
    ), "the plain run kept no variant: the comparison holds nothing"
    assert _run(tmp_path / "gz", mode, gz=True) == plain
