"""`krewlyzer report` — the verdict, the charts, and what neither may claim.

Two properties matter more than the rest, and both are about restraint:

**Absence is never disagreement.** An axis without a PON z-score is *not
assessable*. Counting it as "does not agree" would make a thinner run read as a
healthier one, which is the most dangerous direction for this error.

**A missing measurement is never drawn as a value.** Several columns use zero to
mean "nothing was observed" — `ALT_MeanSize` when no ALT fragment exists,
`mds_e1` before it became NaN. On axes where lower is more abnormal, plotting
that zero puts a fabricated point at the most tumour-like end.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.validate import plots
from krewlyzer.validate.describe import describe_sample
from krewlyzer.validate.verdict import (
    DEFAULT_Z_THRESHOLD,
    Support,
    compute_verdict,
)

SAMPLE = "P-0000000-T01-XS1"


def _write(directory: Path, suffix: str, frame: pd.DataFrame) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(directory / f"{SAMPLE}{suffix}")


@pytest.fixture
def sample_dir(tmp_path: Path) -> Path:
    d = tmp_path / SAMPLE
    _write(
        d,
        ".MDS.parquet",
        pd.DataFrame({"Sample": [SAMPLE], "MDS": [0.81], "mds_z": [-3.2]}),
    )
    _write(
        d,
        ".OCF.offtarget.parquet",
        pd.DataFrame(
            {"tissue": ["Liver", "Tcell"], "OCF": [1.1, 0.4], "ocf_z": [2.9, -1.0]}
        ),
    )
    _write(
        d,
        ".ATAC.parquet",
        pd.DataFrame(
            {
                "label": ["LUAD", "COAD"],
                "count": [10, 20],
                "mean_size": [166.0, 167.0],
                "entropy": [6.9, 6.8],
                "z_score": [0.4, 0.2],
            }
        ),
    )
    _write(
        d, ".FSR.parquet", pd.DataFrame({"region": ["1:0-1"], "short_long_log2": [0.1]})
    )
    return d


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------


def test_all_four_axes_are_reported(sample_dir):
    names = [a.name for a in compute_verdict(sample_dir).axes]
    assert names == [
        "Fragment size",
        "Nuclease cutting",
        "Tissue shedding",
        "Chromatin accessibility",
    ]


def test_direction_is_respected_per_axis(sample_dir):
    """MDS agrees when it goes *down*; OCF when it goes *up*.

    A single shared comparison would have marked `mds_z = -3.2` as opposing —
    the exact inversion that had MDS documented backwards for a year.
    """
    by = {a.name: a for a in compute_verdict(sample_dir).axes}
    assert (
        by["Nuclease cutting"].support is Support.AGREES
    ), "lower MDS is the tumour direction"
    assert (
        by["Tissue shedding"].support is Support.AGREES
    ), "higher OCF is the tumour direction"
    assert by["Chromatin accessibility"].support is Support.QUIET


def test_an_axis_moving_the_wrong_way_is_reported_not_hidden(tmp_path):
    d = tmp_path / SAMPLE
    _write(d, ".MDS.parquet", pd.DataFrame({"mds_z": [4.0]}))
    axis = {a.name: a for a in compute_verdict(d).axes}["Nuclease cutting"]
    assert axis.support is Support.OPPOSES


def test_a_missing_axis_is_not_assessable_rather_than_disagreeing(tmp_path):
    """The most dangerous error to get wrong, so it gets its own test."""
    d = tmp_path / SAMPLE
    d.mkdir()
    verdict = compute_verdict(d)
    assert all(a.support is Support.NOT_ASSESSABLE for a in verdict.axes)
    assert verdict.assessable == []
    assert "no axis had the inputs it needs" in verdict.summary
    assert verdict.headline == "No axis could be assessed"


def test_not_assessable_axes_are_excluded_from_the_denominator(sample_dir):
    """ "2 of 4" when two cannot be judged would understate agreement."""
    verdict = compute_verdict(sample_dir)
    assert len(verdict.assessable) == 4
    assert verdict.headline == "2 of 4 assessable axes agree"


def test_one_agreeing_axis_is_called_a_hypothesis(tmp_path):
    d = tmp_path / SAMPLE
    _write(d, ".MDS.parquet", pd.DataFrame({"mds_z": [-5.0]}))
    assert "hypothesis" in compute_verdict(d).summary


def test_the_threshold_is_an_argument_not_a_constant(sample_dir):
    """It is conventional and tunable, and the report says so."""
    loose = compute_verdict(sample_dir, z_threshold=0.1)
    strict = compute_verdict(sample_dir, z_threshold=10.0)
    assert len(loose.agreeing) > len(strict.agreeing)
    assert strict.agreeing == []
    assert DEFAULT_Z_THRESHOLD == 2.0


def test_the_verdict_never_claims_a_clinical_finding(sample_dir):
    """It reports convergence between measurements, nothing more."""
    text = compute_verdict(sample_dir).summary.lower()
    for word in ("diagnos", "positive for", "patient has", "malignan"):
        assert word not in text


# ---------------------------------------------------------------------------
# Charts refuse rather than approximate
# ---------------------------------------------------------------------------


def test_a_constant_column_is_reported_not_plotted():
    """A flat line looks like a measurement; "constant" does not."""
    fsr = pd.DataFrame({"short_long_log2": [0.5] * 40})
    chart = plots.short_long_spread(fsr)
    assert not chart.drawn
    assert "constant" in chart.reason


def test_an_absent_table_gives_a_reason_not_an_empty_axis():
    chart = plots.fragment_size_density(None)
    assert not chart.drawn and chart.reason
    chart = plots.short_long_spread(pd.DataFrame({"region": ["1:0-1"]}))
    assert "pon-model" in chart.reason.lower()


def test_charts_are_tagged_with_the_table_they_explain():
    """They render inside that table's section, not in a gallery."""
    for chart in plots.build_charts({}):
        assert chart.suffix.startswith(".") and chart.suffix.endswith(".parquet")


def test_every_act_of_the_report_structure_is_covered():
    """The plan is a five-act structure; earlier revisions covered seven of it."""
    suffixes = {c.suffix for c in plots.build_charts({})}
    for expected in (
        ".FSD.parquet",
        ".FSC.parquet",
        ".FSR.parquet",
        ".mFSD.parquet",
        ".EndMotif.parquet",
        ".EndMotif1mer.parquet",
        ".MDS.gene.parquet",
        ".WPS_background.parquet",
        ".ATAC.parquet",
        ".TFBS.parquet",
        ".UXM.parquet",
    ):
        assert expected in suffixes, f"no chart for {expected}"


def test_a_drawn_chart_renders_its_title():
    """Titles were only emitted on the not-drawn branch.

    Every figure that worked arrived without a heading, leaving the reader to
    infer from the caption what they were looking at.
    """
    from krewlyzer.validate.htmlreport import _figure_html

    drawn = plots.fragment_size_density(
        pd.DataFrame({"region": ["chr1"], "65-69": [10.0], "70-74": [20.0]})
    )
    assert drawn.drawn
    assert "Fragment size distribution" in _figure_html(drawn, 0)

    absent = plots.fragment_size_density(None)
    assert "Fragment size distribution" in _figure_html(absent, 1)


# ---------------------------------------------------------------------------
# mFSD: the zero that is not a measurement
# ---------------------------------------------------------------------------


def _mfsd(rows) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_a_variant_with_no_alt_fragment_is_never_placed_on_the_size_axis():
    """`ALT_MeanSize` is 0 when nothing was observed.

    `0 - REF` is about -190 bp, which on this axis is the most tumour-like
    value there is. It has to be shown off-scale or not at all.
    """
    df = _mfsd(
        [
            dict(Chrom="17", Pos=1, REF_MeanSize=190.0, ALT_MeanSize=0.0, ALT_Count=0),
            dict(
                Chrom="5", Pos=2, REF_MeanSize=190.0, ALT_MeanSize=196.0, ALT_Count=58
            ),
        ]
    )
    chart = plots.mutant_fragment_sizes(df)
    assert chart.drawn
    measured, unobserved = chart.figure.data[0], chart.figure.data[1]
    assert list(measured.x) == [6.0], "only the observed variant is on the size axis"
    assert list(unobserved.x) == [0] and list(unobserved.y) == [0]
    assert "not a measurement" in unobserved.hovertemplate


def test_unobserved_variants_are_shown_rather_than_dropped():
    """Their absence is a finding; silently omitting them hides it."""
    df = _mfsd(
        [
            dict(Chrom="1", Pos=1, REF_MeanSize=190.0, ALT_MeanSize=0.0, ALT_Count=0),
            dict(
                Chrom="2", Pos=2, REF_MeanSize=190.0, ALT_MeanSize=170.0, ALT_Count=30
            ),
        ]
    )
    chart = plots.mutant_fragment_sizes(df)
    assert "no ALT-supporting fragment" in chart.caption
    assert len(chart.figure.data) == 2


def test_mfsd_scales_to_many_variants_without_growing_unboundedly():
    """One row per variant made a 12,000-pixel figure at n=300."""
    many = _mfsd(
        [
            dict(
                Chrom="1",
                Pos=i,
                REF_MeanSize=190.0,
                ALT_MeanSize=180.0 + i % 20,
                ALT_Count=10 + i,
            )
            for i in range(300)
        ]
    )
    few = _mfsd(
        [dict(Chrom="1", Pos=1, REF_MeanSize=190.0, ALT_MeanSize=180.0, ALT_Count=10)]
    )
    tall, short = (
        plots.mutant_fragment_sizes(many).figure.layout.height,
        plots.mutant_fragment_sizes(few).figure.layout.height,
    )
    assert tall == short, "height must not scale with variant count"
    assert tall < 600


def test_a_log_axis_is_not_used_when_it_would_hide_the_zero_points():
    """Zero cannot be drawn on a log axis, so those markers would vanish."""
    with_zero = _mfsd(
        [
            dict(Chrom="1", Pos=1, REF_MeanSize=190.0, ALT_MeanSize=0.0, ALT_Count=0),
            dict(Chrom="2", Pos=2, REF_MeanSize=190.0, ALT_MeanSize=170.0, ALT_Count=5),
            dict(
                Chrom="3", Pos=3, REF_MeanSize=190.0, ALT_MeanSize=175.0, ALT_Count=5000
            ),
        ]
    )
    assert plots.mutant_fragment_sizes(with_zero).figure.layout.yaxis.type != "log"


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def test_the_report_renders_end_to_end(sample_dir):
    from krewlyzer.core.output_utils import read_table
    from krewlyzer.validate.htmlreport import render_html

    described = describe_sample(sample_dir)
    tables = {
        t.suffix: read_table(sample_dir / f"{SAMPLE}{t.suffix}")
        for t in described.tables
    }
    html = render_html(
        described, compute_verdict(sample_dir), plots.build_charts(tables), tables
    )

    assert html.startswith("<!doctype html>") and html.rstrip().endswith("</html>")
    assert 'data-mode="dark"' in html, "theme toggle missing"
    assert "no panel of normals" in html.lower() or "Panel of normals" in html


def test_the_report_does_not_claim_what_kreview_reads(sample_dir):
    """That is a fact about another repository, kept in sync by nothing.

    The report says what *this* repo gates on instead, which is true by
    construction.
    """
    from krewlyzer.validate.htmlreport import render_html

    html = render_html(describe_sample(sample_dir), compute_verdict(sample_dir), [], {})
    assert "read downstream" not in html
    assert ">gated<" in html or ">inventory<" in html


def test_a_run_without_a_pon_is_flagged_prominently(tmp_path):
    """Without one, most of the interpretable surface is absent."""
    from krewlyzer.validate.htmlreport import detect_pon, render_html

    d = tmp_path / SAMPLE
    _write(d, ".MDS.parquet", pd.DataFrame({"Sample": [SAMPLE], "MDS": [0.8]}))
    described = describe_sample(d)
    assert detect_pon(described) is False
    html = render_html(described, compute_verdict(d), [], {})
    assert "No panel of normals in this run" in html


def test_a_pon_run_is_detected_from_its_output(sample_dir):
    from krewlyzer.validate.htmlreport import detect_pon

    assert detect_pon(describe_sample(sample_dir)) is True


# ---------------------------------------------------------------------------
# Whole-genome runs
# ---------------------------------------------------------------------------
#
# A panel run splits OCF into on- and off-target; a whole-genome run writes
# neither and emits a plain `.OCF.parquet`. Every reader here originally looked
# only for the panel pair, so a WGS sample reported OCF absent with the file
# sitting in the directory -- and losing an axis makes the verdict *look*
# stronger, because the denominator shrinks with it.


def _wgs_ocf(tmp_path: Path) -> Path:
    d = tmp_path / SAMPLE
    _write(
        d,
        ".OCF.parquet",
        pd.DataFrame(
            {
                "tissue": ["Liver", "Tcell", "Lung"],
                "OCF": [12.0, 4.0, 6.0],
                "ocf_z": [3.4, -0.2, 0.8],
            }
        ),
    )
    return d


def test_a_whole_genome_ocf_table_is_read(tmp_path):
    directory = _wgs_ocf(tmp_path)
    axis = next(
        a for a in compute_verdict(directory).axes if a.name == "Tissue shedding"
    )
    assert axis.support is Support.AGREES, axis.reason
    assert axis.value == pytest.approx(3.4)
    assert "Liver" in (axis.detail or "")


def test_a_whole_genome_run_keeps_the_tissue_axis_in_the_denominator(tmp_path):
    """The failure mode is silent: fewer axes, so a higher fraction agreeing."""
    verdict = compute_verdict(_wgs_ocf(tmp_path))
    assert "Tissue shedding" in {a.name for a in verdict.assessable}


def test_the_panel_split_still_wins_when_both_exist(tmp_path):
    """Off-target is the unbiased view where capture applies, so it ranks first."""
    directory = _wgs_ocf(tmp_path)
    _write(
        directory,
        ".OCF.offtarget.parquet",
        pd.DataFrame({"tissue": ["Colon"], "OCF": [9.0], "ocf_z": [5.5]}),
    )
    axis = next(
        a for a in compute_verdict(directory).axes if a.name == "Tissue shedding"
    )
    assert axis.value == pytest.approx(5.5)
    assert "Colon" in (axis.detail or "")


def test_the_ocf_chart_reads_the_whole_genome_table(tmp_path):
    """verdict and plots must not disagree about which file exists."""
    directory = _wgs_ocf(tmp_path)
    described = describe_sample(directory)
    from krewlyzer.core.output_utils import read_table

    tables = {
        t.suffix: read_table(directory / f"{described.sample_id}{t.suffix}")
        for t in described.tables
    }
    ocf_charts = [c for c in plots.build_charts(tables) if "OCF" in c.suffix]
    assert ocf_charts, "no chart is tagged with an OCF table"
    assert all(c.drawn for c in ocf_charts), [c.reason for c in ocf_charts]


def test_every_optional_output_can_be_interpreted():
    """An ungated table still gets a section, so it still needs a meaning.

    All six shipped without one: `report` drew the mFSD chart and put no
    explanation beside it, because the registry test keyed on CONTRACT alone.
    """
    from krewlyzer.validate.contract import NOT_CONSUMED
    from krewlyzer.validate.meaning import MEANINGS

    missing = sorted(s for s in NOT_CONSUMED if s not in MEANINGS)
    assert not missing, f"optional outputs with no interpretation: {missing}"
