"""What each tool writes is what its tool-level Nextflow module declares, in every format.

A module's `output:` block decides what Nextflow keeps and publishes. A table the tool writes
and the module does not declare is dropped from `--outdir` without a word. A required pattern
the tool does not write fails the task: eight modules required `*.X.tsv` under the Parquet
default, so tool-level mode could not run at all.

Each tool runs on the synthetic sample, staged the way `tool_level/main.nf` stages it (the BED
and extract's GC factor tables beside it), in three configurations that between them write every
extension and both modes. The patterns are read from the modules (`tests/nextflow_modules.py`),
so a new table or a renamed pattern fails here rather than on the cluster.
"""

from __future__ import annotations

import fnmatch
from pathlib import Path
from typing import Dict, List

import pytest
from typer.testing import CliRunner

import nextflow_modules as nf
from conftest import requires_data
from krewlyzer.cli import app

from .runall import write_arms, write_vcf
from .synth import BASELINE, make_sample

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


def _run(*args) -> None:
    result = CliRunner().invoke(app, [str(a) for a in args])
    assert (
        result.exit_code == 0
    ), f"krewlyzer {args[0]} failed:\n{result.output[-2000:]}"


def _stage(extracted: Path, folder: Path) -> None:
    folder.mkdir(parents=True)
    for f in extracted.iterdir():
        if f.name.startswith("S.bed.gz") or "correction_factors" in f.name:
            (folder / f.name).symlink_to(f)


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
    arms = write_arms(s.directory / "arms.validated.bed", s)
    mode = ["-T", s.assets.targets, "-A", "xs1"] if panel else ["--skip-target-regions"]

    extracted = root / "extract"
    _run(
        "extract",
        "-i",
        s.bam,
        "-r",
        s.reference,
        "-o",
        extracted,
        "-s",
        "S",
        *mode,
        *fmt,
    )
    # Extraction's own outputs are everything it wrote but its logs.
    found = {"extract": _left(extracted, ("*.log", "*.json"))}

    bed_tools = {
        "fsc": ["fsc", "-b", s.assets.bins],
        "fsr": ["fsr", "-b", s.assets.bins],
        "fsd": ["fsd", "-a", arms],
        "ocf": ["ocf", "-r", s.assets.regions],
        "wps": ["wps", "--wps-anchors", s.assets.anchors, "-B", s.assets.alu],
        "region_entropy": [
            "region-entropy",
            "--tfbs-regions",
            s.assets.regions,
            "--atac-regions",
            s.assets.regions,
        ],
    }
    for module, (command, *args) in bed_tools.items():
        folder = root / module
        _stage(extracted, folder)
        _run(
            command,
            "-i",
            folder / "S.bed.gz",
            "-o",
            folder,
            "-s",
            "S",
            *args,
            "--skip-pon",
            *mode,
            *fmt,
        )
        found[module] = _left(folder)

    vcf = write_vcf(root / "v.vcf")
    bam_tools = {
        "motif": ["motif", "-i", s.bam, "-r", s.reference, "--skip-pon", *mode],
        "mfsd": ["mfsd", "-i", s.bam, "-V", vcf, "-g", s.reference],
        "uxm": ["uxm", "-i", s.bam],
    }
    for module, args in bam_tools.items():
        folder = root / module
        _run(*args, "-o", folder, "--sample-name", "S", *fmt)
        found[module] = _left(folder)
    if panel:  # region-MDS reads an assay's gene BED
        folder = root / "region_mds"
        _run(
            "region-mds",
            s.bam,
            s.reference,
            folder,
            "--sample-name",
            "S",
            "-a",
            "xs1",
            "--skip-pon",
            *fmt,
        )
        found["region_mds"] = _left(folder)
    return found


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
