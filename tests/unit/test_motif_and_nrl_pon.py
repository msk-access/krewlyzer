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


# ---------------------------------------------------------------------------
# Breakpoint motifs are not end motifs
# ---------------------------------------------------------------------------


def test_breakpoint_motifs_are_scored_against_their_own_baseline(tmp_path):
    """Until 0.9.0 all four motif tables used `mds_baseline`.

    That is an **end-motif, genome-wide** block. An end motif is the 4-mer at
    the fragment's 5′ terminus — what the nuclease left. A breakpoint motif
    spans the cut site and includes reference bases *not present in the
    fragment*. On a real sample the two frequency vectors correlate 0.696, and
    `BreakPointMotif` came out at median |z| 5.85 on XS1 and 11.25 on XS2, with
    70 and 136 of 256 motifs beyond |z| = 10.

    A correctly fitted baseline gives a median near 0.67 — the half-normal
    median — because that is what scoring a sample against its own cohort
    means.

    Full 256-mer baselines, not a two-motif toy: the scorer renormalises over
    the shared motifs, so a partial baseline makes the renormalisation dominate
    and the test measure the wrong thing. The code says so itself, warning when
    the shared mass is small.
    """
    rng = np.random.default_rng(0)
    kmers = [
        a + b + c + d for a in "ACGT" for b in "ACGT" for c in "ACGT" for d in "ACGT"
    ]

    # Two genuinely different distributions, as end and breakpoint motifs are.
    end_p = rng.dirichlet(np.ones(256) * 40)
    bp_p = rng.dirichlet(np.ones(256) * 40)
    sample = bp_p * (1 + rng.normal(0, 0.02, 256))
    sample /= sample.sum()

    path = pathlib.Path(tmp_path) / "s.BreakPointMotif.parquet"
    pd.DataFrame({"Motif": kmers, "Frequency": sample}).to_parquet(path)

    pon = _motif_pon(
        dict(zip(kmers, end_p)),
        dict(zip(kmers, end_p * 0.02)),
    )
    pon.breakpoint_motif_baseline = MdsBaseline(
        kmer_expected=dict(zip(kmers, bp_p)),
        kmer_std=dict(zip(kmers, bp_p * 0.02)),
        mds_mean=float("nan"),
        mds_std=float("nan"),
    )

    def _score(attr):
        apply_motif_pon(
            path,
            pon,
            output_base=path.with_suffix(""),
            output_format="parquet",
            baseline_attr=attr,
        )
        z = pd.read_parquet(path)["frequency_z"].abs()
        return float(z.median()), int((z > 10).sum())

    right_median, right_big = _score("breakpoint_motif_baseline")
    wrong_median, wrong_big = _score("mds_baseline")

    assert right_median < 2.0, (
        f"against its own baseline the median |z| should sit near 0.67, "
        f"got {right_median:.2f}"
    )
    assert (
        right_big == 0
    ), f"{right_big} motifs beyond |z|=10 against the right baseline"
    assert wrong_median > right_median * 3, (
        f"the end-motif baseline should be visibly worse "
        f"({wrong_median:.2f} vs {right_median:.2f}); if not, this fixture no "
        "longer reproduces the defect"
    )
    assert wrong_big > 20, f"only {wrong_big} motifs misplaced by the wrong baseline"


def test_every_motif_table_gets_its_own_baseline():
    """Three of the four tables were scored against the wrong block.

    `mds_baseline_ontarget` already existed and was already used for the
    whole-sample MDS z — just not for the per-motif frequencies.
    """
    import inspect

    from krewlyzer.core import sample_processor

    source = inspect.getsource(sample_processor)
    for key, attr in (
        ('("edm", "mds_baseline")', "edm"),
        ('("bpm", "breakpoint_motif_baseline")', "bpm"),
        ('("edm_ontarget", "mds_baseline_ontarget")', "edm_ontarget"),
        ('("bpm_ontarget", "breakpoint_motif_baseline_ontarget")', "bpm_ontarget"),
    ):
        assert key in source, f"{attr} is no longer paired with its own baseline"


def test_a_breakpoint_record_carries_no_mds_score(tmp_path):
    """MDS is defined on end motifs.

    Computing one over breakpoint 4-mers would be a different statistic under
    the same name, so the collector returns None rather than a number.
    """
    from krewlyzer.pon.from_outputs import _motif_record

    d = pathlib.Path(tmp_path)
    pd.DataFrame({"Motif": ["ACGT"], "Frequency": [1.0]}).to_parquet(
        d / "S.BreakPointMotif.parquet"
    )
    pd.DataFrame({"Motif": ["ACGT"], "Frequency": [1.0]}).to_parquet(
        d / "S.EndMotif.parquet"
    )

    bp = _motif_record(d, "S", table="BreakPointMotif")
    end = _motif_record(d, "S", table="EndMotif")
    assert bp is not None and bp["mds"] is None
    assert end is not None and end["mds"] is not None


def test_a_pon_without_the_breakpoint_block_gives_no_z(tmp_path):
    """Every PON built before 0.9.0 lacks it.

    The honest result is no `frequency_z` column at all — which is visibly
    different from EndMotif, which does get one. That asymmetry says "no
    baseline for this table", which is exactly the situation.
    """
    path = tmp_path / "s.BreakPointMotif.parquet"
    pd.DataFrame({"Motif": ["ACGT"], "Frequency": [0.007]}).to_parquet(path)

    scored = apply_motif_pon(
        path,
        _motif_pon({"ACGT": 0.0008}, {"ACGT": 4e-5}),
        output_base=path.with_suffix(""),
        output_format="parquet",
        baseline_attr="breakpoint_motif_baseline",
    )
    assert scored == 0
    assert "frequency_z" not in pd.read_parquet(path).columns
