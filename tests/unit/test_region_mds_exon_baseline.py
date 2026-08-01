"""`MDS.exon` had no baseline at all, so its score had nothing to mean.

It is the finest localisation krewlyzer produces — which exon, not which
gene — 1,006 rows on xs2 and 1,725 on xs1, and it shipped a raw MDS value
with nothing to compare it against. Every other feature table carries a
PON-derived column; this one carried none because no block existed.

Measured on the 0.8.3 corpus before building it, and it refuted the
assumption that exon-level data would be sparse: every exon appears in every
sample of its assay (7/7 xs1, 19/19 xs2) with a measurable spread, and under
0.25% carry fewer than 10 fragments. So this needs no sparsity machinery —
only the same NaN-not-floor rule as everywhere else.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from krewlyzer.pon.build import (
    MIN_SAMPLES_PER_KEY,
    _compute_region_mds_exon_baseline,
    _save_pon_model,
)
from krewlyzer.pon.model import PonModel, RegionMdsExonBaseline

pytestmark = pytest.mark.unit


def _write(directory: Path, index: int, rows: list[tuple[str, str, float]]) -> str:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"s{index}.MDS.exon.tsv"
    pd.DataFrame(rows, columns=["gene", "name", "mds"]).to_csv(
        path, sep="\t", index=False
    )
    return str(path)


def _cohort(directory: Path, per_sample: list[list[tuple[str, str, float]]]):
    return [_write(directory, i, rows) for i, rows in enumerate(per_sample)]


# ---------------------------------------------------------------------------
# building
# ---------------------------------------------------------------------------


def test_the_baseline_is_fitted_from_the_data(tmp_path):
    paths = _cohort(
        tmp_path,
        [
            [("TP53", "ex1", 0.90)],
            [("TP53", "ex1", 0.94)],
            [("TP53", "ex1", 0.98)],
        ],
    )
    entry = _compute_region_mds_exon_baseline(paths).exon_baseline[("TP53", "ex1")]
    assert entry["mds_mean"] == pytest.approx(0.94)
    assert entry["mds_std"] == pytest.approx(np.std([0.90, 0.94, 0.98], ddof=1))
    assert entry["n_samples"] == 3


def test_two_cohorts_give_two_baselines(tmp_path):
    """Anti-degeneracy, applied to the block itself (invariant #1)."""
    a = _compute_region_mds_exon_baseline(
        _cohort(tmp_path / "a", [[("TP53", "ex1", v)] for v in (0.90, 0.92, 0.94)])
    )
    b = _compute_region_mds_exon_baseline(
        _cohort(tmp_path / "b", [[("TP53", "ex1", v)] for v in (0.70, 0.72, 0.74)])
    )
    assert a.exon_baseline[("TP53", "ex1")]["mds_mean"] != pytest.approx(
        b.exon_baseline[("TP53", "ex1")]["mds_mean"]
    )


def test_the_key_is_gene_and_name_together(tmp_path):
    """`name` alone is not unique across genes.

    Panel BEDs reuse exon names -- two genes both having an `ex1` is the
    normal case, not an edge one. Keying on `name` alone would silently pool
    two unrelated distributions into one baseline.
    """
    paths = _cohort(
        tmp_path,
        [
            [("TP53", "ex1", 0.90), ("EGFR", "ex1", 0.50)],
            [("TP53", "ex1", 0.92), ("EGFR", "ex1", 0.52)],
            [("TP53", "ex1", 0.94), ("EGFR", "ex1", 0.54)],
        ],
    )
    baseline = _compute_region_mds_exon_baseline(paths).exon_baseline
    assert baseline[("TP53", "ex1")]["mds_mean"] == pytest.approx(0.92)
    assert baseline[("EGFR", "ex1")]["mds_mean"] == pytest.approx(0.52)


def test_an_exon_below_the_sample_floor_is_dropped(tmp_path):
    """The same floor FSC gene, FSC region and WPS use, named once."""
    paths = _cohort(
        tmp_path,
        [
            [("TP53", "ex1", 0.90), ("RARE", "ex9", 0.10)],
            [("TP53", "ex1", 0.92)],
            [("TP53", "ex1", 0.94)],
        ],
    )
    baseline = _compute_region_mds_exon_baseline(paths).exon_baseline
    assert ("TP53", "ex1") in baseline
    assert (
        "RARE",
        "ex9",
    ) not in baseline, (
        f"an exon seen once was kept despite MIN_SAMPLES_PER_KEY={MIN_SAMPLES_PER_KEY}"
    )


def test_a_renamed_source_column_is_fatal(tmp_path):
    """Refuse, do not build a partial baseline.

    This is the `wps_background` defect: the builder looked for a column that
    did not exist, fell back to a default, and shipped it four times. A
    missing column here stops the build instead.
    """
    path = tmp_path / "s.MDS.exon.tsv"
    pd.DataFrame({"gene": ["TP53"], "exon": ["ex1"], "mds": [0.9]}).to_csv(
        path, sep="\t", index=False
    )
    with pytest.raises(ValueError, match="name"):
        _compute_region_mds_exon_baseline([str(path)])


def test_an_unmeasurable_spread_is_nan_not_floored(tmp_path):
    """Three identical observations measure no spread. NaN, never a floor."""
    paths = _cohort(tmp_path, [[("TP53", "ex1", 0.9)] for _ in range(3)])
    entry = _compute_region_mds_exon_baseline(paths).exon_baseline[("TP53", "ex1")]
    assert entry["mds_std"] == 0.0 or math.isnan(entry["mds_std"])
    assert math.isnan(
        RegionMdsExonBaseline(exon_baseline={("TP53", "ex1"): entry}).compute_zscore(
            "TP53", "ex1", 0.95
        )
    ), "a zero or NaN spread must not yield a finite z"


# ---------------------------------------------------------------------------
# scoring and round-trip
# ---------------------------------------------------------------------------


def test_scoring_separates_the_two_absences():
    """`None` = not in the baseline; NaN = there, but nothing measured.

    Both write NaN to the column. Only the log can tell "rebuild the PON for
    this panel" from "this exon has no variance in the cohort".
    """
    baseline = RegionMdsExonBaseline(
        exon_baseline={
            ("TP53", "ex1"): {"mds_mean": 0.9, "mds_std": 0.01, "n_samples": 5},
            ("TP53", "ex2"): {"mds_mean": 0.9, "mds_std": float("nan"), "n_samples": 5},
        }
    )
    assert baseline.compute_zscore("TP53", "ex1", 0.92) == pytest.approx(2.0)
    assert math.isnan(baseline.compute_zscore("TP53", "ex2", 0.92))
    assert baseline.compute_zscore("EGFR", "ex1", 0.92) is None


def test_the_block_survives_a_parquet_round_trip(tmp_path):
    """The block must load back as it was written, or the rebuild is wasted."""
    model = PonModel(
        schema_version="1.0",
        assay="xs2",
        build_date="2026-01-01",
        n_samples=19,
        reference="ref",
        panel_mode=True,
        target_regions_file="xs2.targets.bed.gz",
        region_mds_exon=RegionMdsExonBaseline(
            exon_baseline={
                ("TP53", "ex1"): {"mds_mean": 0.90, "mds_std": 0.010, "n_samples": 19},
                ("EGFR", "ex7"): {"mds_mean": 0.55, "mds_std": 0.004, "n_samples": 19},
            }
        ),
    )
    out = tmp_path / "t.pon.parquet"
    # `_save_pon_model` is the only writer; `PonModel.save` was a second,
    # metadata-only one and has been removed.
    _save_pon_model(model, out)

    loaded = PonModel.load(out)
    assert loaded.region_mds_exon is not None, "the block did not survive the write"
    got = loaded.region_mds_exon.exon_baseline
    assert set(got) == {("TP53", "ex1"), ("EGFR", "ex7")}
    assert got[("TP53", "ex1")]["mds_mean"] == pytest.approx(0.90)
    assert got[("EGFR", "ex7")]["mds_std"] == pytest.approx(0.004)
    assert loaded.region_mds_exon.compute_zscore("TP53", "ex1", 0.92) == pytest.approx(
        2.0
    )
