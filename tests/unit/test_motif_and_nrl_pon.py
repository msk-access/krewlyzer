"""Two baselines that were built into every PON and produced no column.

**k-mer (B3).** `mds_baseline` carries 625 k-mer means and sigmas — every
4-mer over ACGTN — and nothing read them. `EndMotif.parquet` shipped `Motif,
Frequency` alone, so a shift in one motif was invisible unless it moved the
whole-sample MDS enough to notice.

**NRL / periodicity (B4).** This one was fully plumbed and still produced
nothing, for a reason worth recording: `compute_nrl_zscore` defaulted to
`group_id="all"`, and no PON has ever contained a group by that name — they are
`Global_All`, `Chr1_All` … `Family_AluY`. Every lookup missed, the function
returned `None`, and `None` was indistinguishable from "this PON has no
baseline". Two defects had to stack for `nrl_z` to be absent: a fabricated
baseline (`4cd634b`) *and* a default naming nothing.
"""

from __future__ import annotations

import math
import pathlib

import numpy as np
import pandas as pd
import pytest

from krewlyzer.core.motif_pon import apply_motif_pon
from krewlyzer.pon.model import (
    GENOME_WIDE_GROUP,
    MdsBaseline,
    PonModel,
    WpsBackgroundBaseline,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# B4: the group id that matched nothing
# ---------------------------------------------------------------------------


def _background(**overrides):
    groups = pd.DataFrame(
        {
            "group_id": [GENOME_WIDE_GROUP, "Chr1_All", "Family_AluY"],
            "nrl_mean": [190.0, 188.0, 192.0],
            "nrl_std": [5.0, 4.0, 6.0],
            "periodicity_mean": [0.40, 0.38, 0.42],
            "periodicity_std": [0.05, 0.04, 0.06],
        }
    )
    for key, value in overrides.items():
        groups[key] = value
    return WpsBackgroundBaseline(groups=groups)


def test_the_default_group_resolves():
    """The whole defect: the old default was `"all"`, which never existed."""
    baseline = _background()
    assert baseline.compute_nrl_zscore(200.0) == pytest.approx(2.0)
    assert GENOME_WIDE_GROUP in set(baseline.groups["group_id"])
    assert GENOME_WIDE_GROUP != "all"


def test_a_group_that_is_absent_is_still_none():
    """Distinguishable from a group present with an unmeasurable spread."""
    assert _background().compute_nrl_zscore(200.0, "NoSuchGroup") is None


def test_a_group_with_no_spread_is_nan_not_none():
    baseline = _background(nrl_std=0.0)
    assert math.isnan(baseline.compute_nrl_zscore(200.0))


def test_every_group_can_be_scored_not_only_the_genome_wide_one():
    """The baseline holds 28 groups; only `Global_All` was ever scored.

    Per-chromosome and per-Alu-family NRL drift is the reason those 27 other
    baselines are built at all.
    """
    baseline = _background()
    scores = [
        baseline.compute_nrl_zscore(200.0, group)
        for group in baseline.groups["group_id"]
    ]
    assert all(s is not None for s in scores)
    assert len(set(np.round(scores, 6))) > 1, "each group has its own baseline"


# ---------------------------------------------------------------------------
# B3: per-motif z, and the normalisation the join must correct
# ---------------------------------------------------------------------------


def _motif_pon(expected, stds):
    return PonModel(
        schema_version="1.0",
        assay="xs1",
        build_date="2026-01-01",
        n_samples=47,
        reference="r",
        panel_mode=True,
        target_regions_file="t",
        mds_baseline=MdsBaseline(
            kmer_expected=expected, kmer_std=stds, mds_mean=0.95, mds_std=0.01
        ),
    )


def _write(tmp_path, motifs, freqs):
    tmp_path = pathlib.Path(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "s.EndMotif.parquet"
    pd.DataFrame({"Motif": motifs, "Frequency": freqs}).to_parquet(path)
    return path


def test_the_join_is_on_the_motif_string_not_position(tmp_path):
    """The baseline has 625 keys and the output 256; they are not aligned.

    Pairing by row order would score every motif against a different motif's
    baseline. Proven directly: reversing the baseline's insertion order must
    change nothing.
    """
    motifs = ["ACGT", "TTTT", "GGGG"]
    freqs = [0.40, 0.20, 0.40]
    expected = {"ACGT": 0.40, "TTTT": 0.20, "GGGG": 0.30, "NACG": 0.10}
    stds = {k: 0.01 for k in expected}

    def score(baseline_expected):
        path = _write(tmp_path / str(id(baseline_expected)), motifs, freqs)
        apply_motif_pon(
            path,
            _motif_pon(baseline_expected, stds),
            output_base=path.with_suffix(""),
            output_format="parquet",
        )
        return pd.read_parquet(path).set_index("Motif")["frequency_z"]

    forward = score(expected)
    reversed_order = score({k: expected[k] for k in reversed(list(expected))})
    pd.testing.assert_series_equal(forward, reversed_order)

    # And the values follow the baseline, not the row: after renormalising to
    # the 0.90 the three shared motifs carry, ACGT's 0.40 sits below its
    # 0.40/0.90 share while GGGG's 0.40 sits above its 0.30/0.90.
    assert forward["ACGT"] < 0
    assert forward["GGGG"] > 0


def test_the_baselines_partial_mass_is_corrected_for(tmp_path):
    """Sample frequencies sum to 1; the baseline's ACGT subset sums to 0.972.

    The missing 2.79% sits in N-containing k-mers the output never reports.
    Comparing them directly biases every z upward — measured on a real sample,
    naive z has median +0.37 against −0.21 after correction.
    """
    motifs = [f"{a}{b}" for a in "ACGT" for b in "ACGT"]
    # baseline covers only these, summing to 0.8; sample sums to 1.0
    expected = {m: 0.8 / len(motifs) for m in motifs}
    expected["NN"] = 0.2  # the unreported mass
    stds = {m: 0.001 for m in expected}
    path = _write(tmp_path, motifs, [1.0 / len(motifs)] * len(motifs))

    apply_motif_pon(
        path,
        _motif_pon(expected, stds),
        output_base=path.with_suffix(""),
        output_format="parquet",
    )
    z = pd.to_numeric(pd.read_parquet(path)["frequency_z"], errors="coerce")
    assert abs(z.median()) < 0.5, (
        f"median z is {z.median():+.2f}; the sample was not put on the "
        "baseline's scale, so every motif carries the normalisation offset"
    )


def test_a_motif_absent_from_the_baseline_is_nan(tmp_path):
    path = _write(tmp_path, ["ACGT", "ZZZZ"], [0.5, 0.5])
    apply_motif_pon(
        path,
        _motif_pon({"ACGT": 1.0}, {"ACGT": 0.01}),
        output_base=path.with_suffix(""),
        output_format="parquet",
    )
    out = pd.read_parquet(path).set_index("Motif")
    assert not math.isnan(out.loc["ACGT", "frequency_z"])
    assert math.isnan(out.loc["ZZZZ", "frequency_z"])


def test_no_baseline_means_no_column_rather_than_a_fabricated_one(tmp_path):
    path = _write(tmp_path, ["ACGT"], [1.0])
    scored = apply_motif_pon(
        path,
        _motif_pon({}, {}),
        output_base=path.with_suffix(""),
        output_format="parquet",
    )
    assert scored == 0
    assert "frequency_z" not in pd.read_parquet(path).columns
