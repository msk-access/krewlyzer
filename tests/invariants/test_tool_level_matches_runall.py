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


#: Compared tables that carry no data on the synthetic sample, by mode, with the reason. Two
#: empty tables are equal, so for these the comparison holds nothing; each entry says why it
#: cannot. A table off this list must carry data, and an entry whose table does is stale.
NO_DATA_ON_THE_SYNTHETIC_SAMPLE = {
    ("wgs", "*.FSC.parquet"): (
        "6 bins of 100 kb, and a WGS row is WGS_CONTINUE_N (50) bins: no row. "
        "Panel mode, one bin per row, compares FSC with rows"
    ),
    ("wgs", "*.FSR.parquet"): "as FSC",
    ("panel", "*.MDS.exon.parquet"): (
        "region-MDS reads the bundled xs1 gene BED, on hg19 coordinates, none on the "
        "synthetic contigs; run-all takes no other"
    ),
    ("panel", "*.MDS.gene.parquet"): "as MDS.exon",
    ("panel", "*.WPS.panel.parquet"): (
        "the bundled xs1 WPS anchors: hg19, none on the synthetic contigs"
    ),
    ("panel", "*.TFBS.ontarget.parquet"): (
        "the bundled xs1 TFBS regions: hg19, none on the synthetic contigs"
    ),
    ("panel", "*.ATAC.ontarget.parquet"): (
        "the bundled xs1 ATAC regions: hg19, none on the synthetic contigs"
    ),
    ("panel", "*.UXM.parquet"): (
        "the bundled hg19 markers, and a BAM without bisulfite conversion: U, X and M are 0"
    ),
    ("wgs", "*.UXM.parquet"): "as in panel mode",
}

#: Where a row is rather than what was measured there: a table of 1725 exons with
#: `n_fragments` and `mds` all 0 has nonzero `start` and `end`, and carries nothing.
POSITION_COLUMNS = {"start", "end", "center"}


def _ignored(name: str, table: Dict[str, str]) -> bool:
    return any(fnmatch.fnmatchcase(name, g) for g in table)


def _carries_data(frame: pd.DataFrame) -> bool:
    measured = [
        c for c in frame.select_dtypes("number").columns if c not in POSITION_COLUMNS
    ]
    return any((frame[c].fillna(0) != 0).any() for c in measured)


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


def test_every_compared_table_carries_data(both, request):
    """`test_every_published_table_equals_run_all` holds nothing for a table empty on both
    sides. Each one is named, with its reason, in NO_DATA_ON_THE_SYNTHETIC_SAMPLE."""
    mode = request.node.callspec.params["both"]
    runall, published = both
    expected = {g for m, g in NO_DATA_ON_THE_SYNTHETIC_SAMPLE if m == mode}
    empty, stale = [], []
    for name in sorted(published):
        if not name.endswith(".parquet") or not (runall / name).exists():
            continue
        has_data = _carries_data(pd.read_parquet(runall / name))
        allowed = any(fnmatch.fnmatchcase(name, g) for g in expected)
        if not has_data and not allowed:
            empty.append(name)
        elif has_data and allowed:
            stale.append(name)
    assert not empty, f"{mode}: compared tables with no data: {empty}"
    assert (
        not stale
    ), f"{mode}: carry data now, off NO_DATA_ON_THE_SYNTHETIC_SAMPLE: {stale}"
