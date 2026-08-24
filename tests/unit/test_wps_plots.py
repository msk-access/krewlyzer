"""The WPS panels must plot base pairs, and say what they rest on.

Two defects, both found by reading a real report rather than the code:

* every WPS profile is 200 bins of 10 bp over a 2000 bp window, and the plots
  used the *bin index* as the x value while labelling the axis "bp" — so a TSS
  dip at bin +15 read as +15 bp when it is +150, and a ±1000 bp window read as
  ±100. A reader has no way to catch that from the figure;
* `capture_mask` marks positions inside a bait and away from its edges, and
  nothing in the Python package read it. On a targeted panel the genome-wide
  anchors are almost entirely off-bait — one XS2 sample had 506 captured bins
  of 11,704,400 — so the curve is an off-target readout presented as a clean
  chromatin profile with nothing saying so.

Neither made a number wrong. Both made a correct number easy to misread, which
is the failure mode this repository treats as the expensive one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from krewlyzer.validate.plots import (
    WPS_BIN_BP,
    nucleosome_profile,
    wps_anchor_profile,
    wps_panel_profile,
)

pytest.importorskip("plotly")

N_BINS = 200


def _profile(depth: float = 1.0) -> list:
    """A dip at the centre, so the degeneracy guard lets the chart draw."""
    x = np.arange(N_BINS) - N_BINS // 2
    return list(depth * (-10.0 * np.exp(-(x**2) / 200.0) - 1.0))


def _wps_frame(captured: int = 0, n: int = 4) -> pd.DataFrame:
    mask = [1] * captured + [0] * (N_BINS - captured)
    return pd.DataFrame(
        {
            "region_type": ["TSS"] * n + ["CTCF"] * n,
            "wps_nuc": [_profile(1.0)] * n + [_profile(0.4)] * n,
            "capture_mask": [mask] * (2 * n),
        }
    )


@pytest.mark.parametrize(
    "builder, frame",
    [
        (wps_anchor_profile, _wps_frame()),
        (wps_panel_profile, _wps_frame()),
        (
            nucleosome_profile,
            pd.DataFrame({"group_id": ["Alu"], "stacked_wps_nuc": [_profile()]}),
        ),
    ],
    ids=["genome-wide", "panel", "alu-background"],
)
def test_the_x_axis_is_base_pairs_not_bin_index(builder, frame):
    """±1000 bp across 200 bins, so the step is WPS_BIN_BP and the span is ±1000."""
    chart = builder(frame)
    assert chart.drawn, chart.reason
    for trace in chart.figure.data:
        x = list(trace.x)
        assert len(x) == N_BINS
        assert x[1] - x[0] == WPS_BIN_BP, "x is bin index, not bp"
        assert x[0] == -(N_BINS // 2) * WPS_BIN_BP == -1000
        assert x[-1] == 990
    assert "(bp)" in chart.figure.layout.xaxis.title.text


def test_bin_width_matches_the_window_the_rust_side_writes():
    """WPS_BIN_BP is WINDOW_SIZE / NUM_BINS; claims.py pins both against wps.rs.

    Restated here so a reader of this file sees where 10 comes from, and so the
    arithmetic itself is checked rather than only the two endpoints.
    """
    assert WPS_BIN_BP * N_BINS == 2000


def test_an_off_bait_profile_says_so():
    chart = wps_anchor_profile(_wps_frame(captured=0))
    assert "Bait coverage:" in chart.caption
    assert "off-target coverage" in chart.caption


def test_a_fully_captured_profile_says_nothing():
    """A WGS run has no baits, so the notion does not apply and a note is noise."""
    chart = wps_anchor_profile(_wps_frame(captured=N_BINS))
    assert "Bait coverage" not in chart.caption


def test_the_note_counts_anchors_not_only_bins():
    """The bin fraction alone misleads, in both directions.

    A WPS window is 2000 bp against a bait of order 100 bp, so even a perfectly
    bait-centred anchor captures a small minority of its bins. Anchors-touching
    -bait is the number that separates a panel anchor set from a genome-wide
    one: 5 of 166 against 33 of 58,522 on the sample this was written from.
    """
    # 3 of 8 anchors carry a captured bin.
    frame = _wps_frame(captured=0, n=4)
    mask_hit = [1] * 20 + [0] * (N_BINS - 20)
    frame.loc[[0, 1, 4], "capture_mask"] = pd.Series([mask_hit] * 3, index=[0, 1, 4])
    caption = wps_anchor_profile(frame).caption
    assert "3 of 8 anchors" in caption
    assert "60 of 1,600 bins" in caption


def test_the_capture_fraction_never_rounds_away_to_zero():
    """506 of 11,704,400 is 0.00432%, which every fixed-decimal format hides.

    Reporting it as "0.00%" reads as a rounding artefact rather than a fact, so
    the note carries counts and a significant-digit percentage.
    """
    frame = _wps_frame(captured=1, n=250)  # 250 captured bins of 100,000
    chart = wps_anchor_profile(frame)
    assert "0.00%" not in chart.caption
    assert "500 of 100,000 bins" in chart.caption


def test_the_note_does_not_call_panel_anchors_unintended():
    """An earlier draft said "not the assay's intended capture".

    That was wrong for `.WPS.panel.parquet`: those anchors are precisely what
    the assay targets, and the low bin fraction is window-versus-bait geometry.
    """
    caption = wps_panel_profile(_wps_frame(captured=1)).caption
    assert "intended capture" not in caption


def test_no_capture_column_means_no_claim_about_capture():
    """The Alu background has no mask; inventing a figure would be worse."""
    frame = _wps_frame().drop(columns=["capture_mask"])
    assert "Bait coverage" not in wps_anchor_profile(frame).caption


def test_the_legend_says_how_many_anchors_each_curve_averages():
    """n is the difference between a mean over 105 anchors and over 34,851."""
    chart = wps_anchor_profile(_wps_frame(n=7))
    assert all("n=7" in (t.name or "") for t in chart.figure.data)


def test_the_two_anchor_sets_are_separate_charts():
    """They are different measurements; only one used to be plotted."""
    a = wps_anchor_profile(_wps_frame())
    b = wps_panel_profile(_wps_frame())
    assert a.suffix == ".WPS.parquet"
    assert b.suffix == ".WPS.panel.parquet"
    assert a.title != b.title
