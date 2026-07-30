"""Two different inputs must produce two different outputs.

The single assertion that would have caught the largest defect in this
codebase. `nrl_bp` was 150.0 for every sample ever produced, `periodicity_score`
0.3333, `adjusted_score` 0.0 -- present, plausible, and constant, which passes
any test that only checks a column exists. Four of the five columns one output
table contributes carried no information for a year.

Cheap to run and needs no real data: synthesize two samples that differ in a
way the metric is supposed to detect, and require the metric to notice.
"""

from __future__ import annotations

import pandas as pd
import pytest

from .pipeline import run_pipeline
from .synth import BASELINE, SHORT_SHIFTED, ULTRA_LONG_HEAVY

pytestmark = [pytest.mark.invariant, pytest.mark.rust, pytest.mark.slow]


def _counts(tmp_path, profile) -> pd.DataFrame:
    run = run_pipeline(tmp_path / profile.name, profile)
    return pd.read_csv(run.output_dir / f"{profile.name}.fsc_counts.tsv", sep="\t")


@pytest.fixture(scope="module")
def baseline_vs_short(tmp_path_factory):
    root = tmp_path_factory.mktemp("antidegeneracy")
    return _counts(root, BASELINE), _counts(root, SHORT_SHIFTED)


def test_size_channels_respond_to_the_size_distribution(baseline_vs_short):
    """Shifting the fragment-length distribution must move the channels.

    A counter wired to a constant, or reading the wrong field, would produce
    identical channel proportions for both profiles.
    """
    baseline, shifted = baseline_vs_short

    def share(df, channel):
        return df[channel].sum() / max(df["total"].sum(), 1)

    assert share(shifted, "core_short") > share(baseline, "core_short"), (
        "a distribution centred at 120bp instead of 167bp must put a larger "
        "share in core_short (101-149bp)"
    )
    assert share(shifted, "mono_nucl") < share(
        baseline, "mono_nucl"
    ), "and a smaller share in mono_nucl (150-220bp)"


def test_channel_profile_is_not_identical_across_inputs(baseline_vs_short):
    """The blunt form: some channel, somewhere, has to differ."""
    baseline, shifted = baseline_vs_short
    channels = [
        "ultra_short",
        "core_short",
        "mono_nucl",
        "di_nucl",
        "long",
        "ultra_long",
    ]

    identical = [c for c in channels if baseline[c].sum() == shifted[c].sum()]
    assert len(identical) < len(channels), (
        "every channel total is identical across two deliberately different "
        "size distributions -- the counter is not reading the data"
    )


def test_long_tail_moves_the_ultra_long_channel(tmp_path):
    """A necrosis-like long tail must be visible where it is supposed to be."""
    baseline = _counts(tmp_path, BASELINE)
    heavy = _counts(tmp_path, ULTRA_LONG_HEAVY)

    assert heavy["ultra_long"].sum() > baseline["ultra_long"].sum(), (
        "a profile with 40% of fragments at ~600bp did not raise ultra_long "
        "relative to a 167bp baseline"
    )


def test_totals_track_the_fragment_count(tmp_path):
    """Guards the opposite failure: a metric that varies for the wrong reason.

    Two profiles with the same fragment count must agree on the total, or the
    differences above could be explained by one sample simply having more
    fragments.
    """
    baseline = _counts(tmp_path, BASELINE)
    shifted = _counts(tmp_path, SHORT_SHIFTED)

    assert BASELINE.n_fragments == SHORT_SHIFTED.n_fragments, (
        "the profiles must be matched on count for the comparison to isolate "
        "the size distribution"
    )
    # Extraction filters 65-1000bp, so allow the tails to differ slightly.
    assert abs(baseline["total"].sum() - shifted["total"].sum()) <= 0.05 * max(
        baseline["total"].sum(), 1
    ), "matched profiles produced very different fragment totals"
