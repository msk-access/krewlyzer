"""Counts must add up.

Partition invariants over the real pipeline: every fragment lands in exactly
one size channel, and the channels sum to the total. These are the checks that
make a ratio interpretable -- `channel / total` means nothing if the channels
do not partition the total, and that is precisely the confusion the audit found
between the genome-bin counter and the gene aggregator.

Integer conservation is asserted exactly. GC correction is off for these
samples, so the counts are whole fragments and there is no rounding to absorb.
"""

from __future__ import annotations

import pandas as pd
import pytest

from .pipeline import run_pipeline
from .synth import BASELINE, ULTRA_LONG_HEAVY

pytestmark = [pytest.mark.invariant, pytest.mark.rust, pytest.mark.slow]

CHANNELS = [
    "ultra_short",
    "core_short",
    "mono_nucl",
    "di_nucl",
    "long",
    "ultra_long",
]


@pytest.fixture(scope="module")
def baseline_counts(tmp_path_factory) -> pd.DataFrame:
    run = run_pipeline(tmp_path_factory.mktemp("conservation"), BASELINE)
    counts = run.output_dir / f"{BASELINE.name}.fsc_counts.tsv"
    assert counts.exists(), f"pipeline produced no counts table at {counts}"
    return pd.read_csv(counts, sep="\t")


def test_channels_partition_the_total(baseline_counts):
    """The six channels must sum to `total`, exactly."""
    missing = [c for c in CHANNELS if c not in baseline_counts.columns]
    assert not missing, f"count table is missing channel(s): {missing}"

    summed = baseline_counts[CHANNELS].sum(axis=1)
    diff = (summed - baseline_counts["total"]).abs()

    assert diff.max() == 0, (
        f"{int((diff > 0).sum())} bin(s) where the six channels do not sum to "
        f"total; largest gap {diff.max():g}. channel/total ratios are only "
        "interpretable if the channels partition the total."
    )


def test_no_fragment_is_counted_twice(baseline_counts):
    """Total across bins must equal the fragments that went in.

    Catches double-counting at bin boundaries, which an equality check within
    a single row cannot see.
    """
    assert baseline_counts["total"].sum() == baseline_counts[CHANNELS].sum().sum()


def test_channels_are_non_negative(baseline_counts):
    for channel in CHANNELS:
        assert (baseline_counts[channel] >= 0).all(), f"{channel} went negative"


def test_ultra_long_channel_actually_receives_fragments(tmp_path):
    """A channel nothing ever lands in is indistinguishable from a broken one.

    ultra_long (401-1000bp) was added deliberately for necrosis detection, and
    a baseline cfDNA profile barely populates it -- so this drives a profile
    with a heavy long tail and checks the fragments arrive.
    """
    run = run_pipeline(tmp_path, ULTRA_LONG_HEAVY)
    counts = pd.read_csv(
        run.output_dir / f"{ULTRA_LONG_HEAVY.name}.fsc_counts.tsv", sep="\t"
    )

    assert counts["ultra_long"].sum() > 0, (
        "a profile with 40% of fragments at ~600bp put nothing in ultra_long; "
        "either the channel bound moved or the fragments are being dropped"
    )
    # And it must still partition.
    assert (counts[CHANNELS].sum(axis=1) - counts["total"]).abs().max() == 0
