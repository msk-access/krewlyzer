"""What each tool writes is what its tool-level Nextflow module declares, in every format.

A module's `output:` block decides what Nextflow keeps and publishes. A table the tool writes
and the module does not declare is dropped from `--outdir` without a word. A required pattern
the tool does not write fails the task: eight modules required `*.X.tsv` under the Parquet
default, so tool-level mode could not run at all.

Each tool runs on the synthetic sample, staged the way `tool_level/main.nf` stages it
(`tool_level.py`), in three configurations that between them write every extension and both
modes. The patterns are read from the modules (`tests/nextflow_modules.py`), so a new table or a
renamed pattern fails here rather than on the cluster.

Blind to a table only some data produces: the synthetic sample has no on-target fragment in an
OCR, so `OCF.ontarget` is never written here. `test_tool_level_modules.py` holds the modules'
declarations to run-all's, which does not depend on the data.
"""

from __future__ import annotations

import fnmatch
from pathlib import Path
from typing import Dict, List

import pytest
import nextflow_modules as nf
from conftest import requires_data

from .runall import write_vcf
from .synth import BASELINE, make_sample
from .tool_level import run_tools

pytestmark = [pytest.mark.invariant, pytest.mark.slow, requires_data]

#: Files a tool leaves that are deliberately not declared, so not published, with the reason.
UNDECLARED = {
    # Rust intermediates `region-entropy` leaves behind (run-all removes them). Not tables
    # downstream reads; the single tool's parity gap with run-all is tracked separately.
    "region_entropy": ("*.raw.tsv",),
    # The counts FSC and FSR are both computed from. `fsc` publishes them; two tasks
    # publishing one file name into one folder would race.
    "fsr": ("*.fsc_counts.*",),
}

#: Inputs staged into the task folder, and what every task writes beside its tables.
_STAGED = ("S.bed.gz", "S.bed.gz.tbi", "*.correction_factors*", "*.log", "*.json")

CONFIGS = {
    "parquet-panel": (["--output-format", "parquet"], True),
    "tsv-gz-wgs": (["--output-format", "tsv", "--compress"], False),
    "both-panel": (["--output-format", "both"], True),
}


def _left(folder: Path, staged=_STAGED) -> List[str]:
    return sorted(
        p.name
        for p in folder.iterdir()
        if p.is_file() and not any(fnmatch.fnmatchcase(p.name, g) for g in staged)
    )


@pytest.fixture(scope="module", params=sorted(CONFIGS))
def written(request, tmp_path_factory) -> Dict[str, List[str]]:
    """``{module: files its tool left}`` for one configuration."""
    fmt, panel = CONFIGS[request.param]
    root = tmp_path_factory.mktemp(request.param)
    s = make_sample(root / "in", BASELINE)
    folders = run_tools(s, root, fmt, panel, write_vcf(root / "v.vcf"))
    return {
        # Extraction's own outputs are everything it wrote but its logs.
        module: _left(folder, ("*.log", "*.json") if module == "extract" else _STAGED)
        for module, folder in folders.items()
    }


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES])
def test_every_file_a_tool_writes_is_declared(written, module):
    if module not in written:
        pytest.skip(f"{module} does not run in this configuration")
    declared = nf.outputs(module)
    allowed = UNDECLARED.get(module, ())
    undeclared = [
        name
        for name in written[module]
        if not any(o.matches(name) for o in declared)
        and not any(fnmatch.fnmatchcase(name, g) for g in allowed)
    ]
    assert not undeclared, (
        f"`{module}` writes {undeclared}, which its module does not declare: "
        "Nextflow would drop them from --outdir"
    )


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES])
def test_every_required_output_is_written(written, module):
    if module not in written:
        pytest.skip(f"{module} does not run in this configuration")
    missing = [
        o.glob
        for o in nf.outputs(module)
        if not o.optional and not any(o.matches(n) for n in written[module])
    ]
    assert not missing, f"`{module}` requires {missing}, which the tool did not write"
