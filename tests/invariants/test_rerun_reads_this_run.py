"""A rerun into a folder an earlier run used ships what a run into a fresh folder ships.

The Rust pass writes FSC counts, OCF and the TFBS/ATAC sync tables as plain TSV, and Python reads
each back to normalise it and write the run's formats. Read back Parquet-first (`read_table`), a
rerun of the same format found the *previous* run's ``.parquet`` beside this run's TSV and shipped
that instead: FSC and FSR computed from the first run's counts, OCF and the sync tables copied
from the first run. The format guard refuses only a change of format, so it stops none of this,
and the stale file carries this run's name and format, so `validate-output` cannot see it.

Each table the first run converted is planted with values that differ from a fresh run's
(invariant 1); the second run must then equal a run into a fresh folder.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from krewlyzer.assets import AssetManager

from .runall import SAMPLE, synth_run_all
from .synth import BASELINE, make_sample

pytestmark = pytest.mark.invariant

#: Tables Python converts from a Rust TSV in run-all, and what the second run ships from them.
PLANTED = (
    "fsc_counts",
    "fsc_counts.ontarget",
    "OCF",
    "OCF.sync",
    "TFBS.sync",
    "ATAC.sync",
)
SHIPPED = PLANTED + ("FSC", "FSR", "FSC.ontarget", "FSR.ontarget")


def _run(sample, out):
    regions = sample.assets.regions
    with pytest.MonkeyPatch.context() as mp:
        # The bundled TFBS/ATAC regions are on hg19 coordinates; the synthetic contigs are not.
        mp.setattr(AssetManager, "tfbs_regions", property(lambda self: regions))
        mp.setattr(AssetManager, "atac_regions", property(lambda self: regions))
        code = synth_run_all(
            sample,
            out,
            output_format="parquet",
            skip_target_regions=False,
            target_regions=sample.assets.targets,
            no_tfbs=False,
            no_atac=False,
        )
    assert code is None, f"run-all exited {code}"


def _table(out, stem):
    return pd.read_parquet(out / f"{SAMPLE}.{stem}.parquet")


def _plant(path):
    """Change every float column, row-dependently, so no per-sample normalisation cancels it."""
    df = pd.read_parquet(path)
    weight = np.arange(len(df)) + 2.0
    for col in df.columns:
        if pd.api.types.is_float_dtype(df[col]):
            df[col] = df[col].fillna(1.0) * weight + 12345.0
        elif pd.api.types.is_integer_dtype(df[col]) and col not in ("start", "end"):
            df[col] = df[col] * weight.astype(df[col].dtype) + 7
    df.to_parquet(path, index=False)


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("rerun")
    sample = make_sample(root / "in", BASELINE)
    fresh, reused = root / "fresh", root / "reused"
    _run(sample, fresh)
    _run(sample, reused)
    for stem in PLANTED:
        _plant(reused / f"{SAMPLE}.{stem}.parquet")
    planted = {stem: _table(reused, stem) for stem in PLANTED}
    _run(sample, reused)
    return fresh, reused, planted


@pytest.mark.parametrize("stem", PLANTED)
def test_the_planted_table_differs_from_a_fresh_one(runs, stem):
    """Without this the comparison below could pass for a test that sees nothing."""
    fresh, _, planted = runs
    assert not planted[stem].equals(
        _table(fresh, stem)
    ), f"{stem}: planting changed nothing"


@pytest.mark.parametrize("stem", SHIPPED)
def test_a_rerun_ships_what_a_fresh_folder_ships(runs, stem):
    fresh, reused, _ = runs
    pd.testing.assert_frame_equal(
        _table(reused, stem),
        _table(fresh, stem),
        obj=f"{stem} after a rerun into a used folder",
    )
