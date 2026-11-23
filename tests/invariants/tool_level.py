"""Run each tool on the synthetic sample the way tool-level Nextflow mode runs it.

`tool_level/main.nf` runs `extract`, then each tool in its own task. A BED-reading tool's task
holds the BED and extract's GC factor tables, staged as links. The tools get the sample's mode
and the format options, and `extract` and `motif`, which read the BAM, the synthetic blacklist
that run-all is given. Shared by the tests that hold tool-level mode to its declarations
(`test_tool_level_outputs.py`) and to run-all (`test_tool_level_matches_runall.py`).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from typer.testing import CliRunner

from krewlyzer.cli import app

from .runall import SAMPLE, write_arms
from .synth import SynthSample


def run(*args) -> None:
    result = CliRunner().invoke(app, [str(a) for a in args])
    assert (
        result.exit_code == 0
    ), f"krewlyzer {args[0]} failed:\n{result.output[-2000:]}"


def mode_options(s: SynthSample, panel: bool) -> List[str]:
    return ["-T", s.assets.targets, "-A", "xs1"] if panel else ["--skip-target-regions"]


def _stage(extracted: Path, folder: Path) -> None:
    folder.mkdir(parents=True)
    for f in extracted.iterdir():
        if f.name.startswith(f"{SAMPLE}.bed.gz") or "correction_factors" in f.name:
            (folder / f.name).symlink_to(f)


def run_tools(
    s: SynthSample, root: Path, fmt: List[str], panel: bool, vcf: Path
) -> Dict[str, Path]:
    """``{module: the folder its task wrote}``, extract's included."""
    mode = mode_options(s, panel)
    arms = write_arms(s.directory / "arms.validated.bed", s)
    folders = {"extract": root / "extract"}
    run(
        "extract",
        "-i",
        s.bam,
        "-r",
        s.reference,
        "-o",
        folders["extract"],
        "-s",
        SAMPLE,
        "-x",
        s.assets.blacklist,
        *mode,
        *fmt,
    )

    bed_tools = {
        "fsc": ["fsc", "-b", s.assets.bins],
        "fsr": ["fsr", "-b", s.assets.bins],
        "fsd": ["fsd", "-a", arms],
        "ocf": ["ocf", "-r", s.assets.regions],
        "wps": ["wps", "--wps-anchors", s.assets.anchors, "-B", s.assets.alu],
        "region_entropy": ["region-entropy"],
    }
    for module, (command, *args) in bed_tools.items():
        folder = folders[module] = root / module
        _stage(folders["extract"], folder)
        run(
            command,
            "-i",
            folder / f"{SAMPLE}.bed.gz",
            "-o",
            folder,
            "-s",
            SAMPLE,
            *args,
            "--skip-pon",
            *mode,
            *fmt,
        )

    bam_tools = {
        "motif": [
            "motif",
            "-i",
            s.bam,
            "-r",
            s.reference,
            "-x",
            s.assets.blacklist,
            "--skip-pon",
            *mode,
        ],
        "mfsd": ["mfsd", "-i", s.bam, "-V", vcf, "-g", s.reference],
        "uxm": ["uxm", "-i", s.bam],
    }
    for module, args in bam_tools.items():
        folders[module] = root / module
        run(*args, "-o", folders[module], "--sample-name", SAMPLE, *fmt)
    if panel:  # run-all and tool-level run region-MDS only with an assay
        folders["region_mds"] = root / "region_mds"
        run(
            "region-mds",
            s.bam,
            s.reference,
            folders["region_mds"],
            "--sample-name",
            SAMPLE,
            "-a",
            "xs1",
            "--skip-pon",
            *fmt,
        )
    return folders
