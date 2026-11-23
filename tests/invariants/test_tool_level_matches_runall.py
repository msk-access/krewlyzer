"""Tool-level mode publishes run-all's tables: the same files, with the same contents.

Invariant #6. Run-all derives some settings inside the CLI rather than in its Nextflow module, so
a test that compares the modules' flags cannot see them. Two did not match: panel FSC/FSR used
5 Mb windows where run-all uses 100 kb bins, and `motif` read blacklisted fragments that run-all
excludes. Only running both and comparing what they wrote catches that.

On the synthetic sample, panel and WGS: run-all once, then each tool as `tool_level/main.nf`
stages it (`tool_level.py`), keeping only the files its module declares (what Nextflow
publishes).

Blind to mFSD's GC correction: run-all corrects mFSD with extract's factors and tool-level does
not (`docs/nextflow/outputs.md`), but on the synthetic sample the correction changes nothing.
"""

from __future__ import annotations

import fnmatch
import gzip
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
import pytest

import nextflow_modules as nf
from conftest import requires_data

from .runall import synth_run_all, write_vcf
from .synth import BASELINE, make_sample
from .tool_level import run_tools

pytestmark = [pytest.mark.invariant, pytest.mark.slow, requires_data]

#: What run-all writes and tool-level does not, deliberately, with the reason.
RUN_ALL_ONLY = {
    "*.metadata.*": "the record: written by the finalize step (the next PR)",
    "*.features.json": "run-all's summary",
    "*.validation.json": "run-all's own validation",
    "*.fingerprint.json": "run-all's cohort fingerprint",
    "*.TFBS.sync.*": "msk-access/krewlyzer#100",
    "*.TFBS.ontarget.sync.*": "msk-access/krewlyzer#100",
    "*.ATAC.sync.*": "msk-access/krewlyzer#100",
    "*.ATAC.ontarget.sync.*": "msk-access/krewlyzer#100",
}


def _ignored(name: str, table: Dict[str, str]) -> bool:
    return any(fnmatch.fnmatchcase(name, g) for g in table)


@pytest.fixture(scope="module", params=["panel", "wgs"])
def both(request, tmp_path_factory) -> Tuple[Path, Dict[str, Path]]:
    """Run-all's folder, and ``{file name: the tool folder that published it}``."""
    panel = request.param == "panel"
    root = tmp_path_factory.mktemp(request.param)
    s = make_sample(root / "in", BASELINE)
    vcf = write_vcf(root / "v.vcf")
    runall = root / "runall"
    mode = (
        dict(skip_target_regions=False, target_regions=s.assets.targets, assay="xs1")
        if panel
        else {}
    )
    code = synth_run_all(
        s,
        runall,
        no_tfbs=False,
        no_atac=False,
        variants=vcf,
        bisulfite_bam=s.bam,
        **mode,
    )
    assert code in (None, 0), f"run-all exited {code}"

    published: Dict[str, Path] = {}
    for module, folder in run_tools(
        s, root, ["--output-format", "parquet"], panel, vcf
    ).items():
        declared = nf.outputs(module)
        for f in folder.iterdir():
            if f.is_symlink() or not any(o.matches(f.name) for o in declared):
                continue
            if module == "extract" and ".metadata." in f.name:
                continue  # excluded from publishing
            assert f.name not in published, f"{f.name} published by two tasks"
            published[f.name] = folder
    return runall, published


def test_tool_level_publishes_the_files_run_all_writes(both):
    runall, published = both
    wrote = {
        f.name
        for f in runall.iterdir()
        if f.is_file()
        and not _ignored(f.name, RUN_ALL_ONLY)
        and not f.name.endswith(".log")
    }
    assert not sorted(
        wrote - set(published)
    ), "run-all writes, tool-level does not publish"
    assert not sorted(
        set(published) - wrote
    ), "tool-level publishes, run-all does not write"


def _read(path: Path):
    if path.name.endswith(".parquet"):
        return pd.read_parquet(path).reset_index(drop=True)
    if path.name.endswith(".bed.gz"):
        return gzip.open(path).read()
    return path.read_bytes()


def test_every_published_table_equals_run_all(both):
    runall, published = both
    differ = []
    for name, folder in sorted(published.items()):
        if not (runall / name).exists():
            continue
        ours, theirs = _read(folder / name), _read(runall / name)
        if isinstance(ours, pd.DataFrame):
            try:
                pd.testing.assert_frame_equal(
                    ours, theirs, check_dtype=False, rtol=1e-9
                )
            except AssertionError:
                differ.append(f"{name} ({len(ours)} rows, run-all {len(theirs)})")
        elif ours != theirs:
            differ.append(name)
    assert not differ, f"tool-level tables differ from run-all's: {differ}"
