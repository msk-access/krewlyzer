"""Gene-level and genome-level FSC must use the same size bands.

They did not. The genome-bin counter split at <=100/<=149/<=220/<=260/<=400
with a sixth ultra_long channel; the gene/region aggregator split at
<100/<150/<260/<400 with no sixth channel. A column named ``di_nucl`` therefore
meant 221-260bp in ``FSC.parquet`` and 260-399bp in ``FSC.gene.parquet``, and
gene ``long`` held what the genome table called ``ultra_long``.

Three separate places already documented the *correct* bands -- the feature
table in ``fsc.md``, the column table in ``output-files.md``, and a line in the
same file claiming "146 genes x 6 channels". Only the implementation disagreed,
so this is a correction rather than a redefinition.

``rust/src/fsc.rs`` pins the classifier in unit tests. These tests pin what
comes out the other end of the real writer, which is what a consumer actually
reads.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd
import pytest

from krewlyzer import _core

from .pipeline import extract_to_bedgz
from .synth import BASELINE, Contig, SynthProfile, make_sample

pytestmark = [pytest.mark.invariant, pytest.mark.rust, pytest.mark.slow]

CONTIGS = (Contig("chr1", 200_000),)

CHANNELS = [
    "ultra_short",
    "core_short",
    "mono_nucl",
    "di_nucl",
    "long",
    "ultra_long",
]
RATIOS = [f"{c}_ratio" for c in CHANNELS]

#: A profile spanning every band, so no channel is empty by construction.
#:
#: One mode per band, centred well inside it with a tight spread, so ordinary
#: jitter cannot starve a channel and make the partition checks vacuous.
ALL_BANDS = SynthProfile(
    name="all_bands",
    n_fragments=1800,
    length_mixture=(
        (82, 5, 1.0),  # ultra_short  65-100
        (125, 6, 1.0),  # core_short  101-149
        (185, 8, 1.0),  # mono_nucl   150-220
        (240, 5, 1.0),  # di_nucl     221-260
        (330, 15, 1.0),  # long        261-400
        (600, 40, 1.0),  # ultra_long  401-1000
    ),
    seed=20260728,
)


def _aggregate(directory: Path, profile: SynthProfile, mode: str) -> pd.DataFrame:
    """Run the real gene/region aggregator and read back what it wrote."""
    sample = make_sample(directory / "input", profile, CONTIGS)
    bed_gz = extract_to_bedgz(sample, directory / "extract")

    out = directory / f"{mode}.tsv"
    written = _core.aggregate_by_gene(
        str(bed_gz),
        str(sample.assets.genes),
        str(out),
        None,  # no GC model for a synthetic reference
        mode,
    )
    assert written > 0, f"{mode} aggregation wrote nothing; the test is vacuous"
    assert out.exists(), f"{mode} aggregation reported {written} rows but wrote no file"
    return pd.read_csv(out, sep="\t")


@pytest.fixture(scope="module")
def gene_table(tmp_path_factory) -> pd.DataFrame:
    return _aggregate(tmp_path_factory.mktemp("gene"), ALL_BANDS, "gene")


@pytest.fixture(scope="module")
def region_table(tmp_path_factory) -> pd.DataFrame:
    return _aggregate(tmp_path_factory.mktemp("region"), ALL_BANDS, "region")


@pytest.mark.parametrize("table", ["gene_table", "region_table"])
def test_all_six_channels_are_emitted(table, request):
    """Five channels is the old shape; ultra_long must be present."""
    df = request.getfixturevalue(table)
    missing = [c for c in CHANNELS + RATIOS if c not in df.columns]
    assert not missing, (
        f"{table} is missing {missing}; the gene path emitted five channels "
        "before the bands were aligned, folding everything over 400bp into 'long'"
    )


@pytest.mark.parametrize("table", ["gene_table", "region_table"])
def test_channels_partition_the_total(table, request):
    """`channel / total` is meaningless unless the channels partition total."""
    df = request.getfixturevalue(table)
    # Counts are written at 2dp; six of them can drift by at most 6 * 0.005.
    tol = 6 * 0.005
    diff = (df[CHANNELS].sum(axis=1) - df["total"]).abs()
    assert diff.max() <= tol, (
        f"{int((diff > tol).sum())} row(s) where the six channels do not sum "
        f"to total (largest gap {diff.max():g}, tolerance {tol:g})"
    )


@pytest.mark.parametrize("table", ["gene_table", "region_table"])
def test_the_six_ratios_sum_to_one(table, request):
    """The relation the contract gate asserts, checked against real output."""
    df = request.getfixturevalue(table)
    summed = df[RATIOS].sum(axis=1)
    # Only rows with fragments; an empty region divides by a 1e-9 floor.
    populated = summed[df["total"] > 0]
    assert not populated.empty, "no populated rows; the fixture is vacuous"
    assert (
        populated - 1.0
    ).abs().max() <= 6 * 5e-5, (
        f"six ratios sum to {populated.min():.6f}..{populated.max():.6f}, not 1"
    )


@pytest.mark.parametrize("table", ["gene_table", "region_table"])
def test_the_five_legacy_ratios_sum_to_one_minus_ultra_long(table, request):
    """What a consumer reading only the pre-0.9.0 columns will see.

    kreview selects the five ratios that existed before ``ultra_long_ratio``.
    Those now sum to ``1 - ultra_long_ratio``, and the gate says so explicitly
    rather than letting a downstream renormalisation quietly use the wrong base.
    """
    df = request.getfixturevalue(table)
    legacy = df[[f"{c}_ratio" for c in CHANNELS[:-1]]].sum(axis=1)
    expected = 1.0 - df["ultra_long_ratio"]
    populated = (legacy - expected).abs()[df["total"] > 0]
    assert populated.max() <= 6 * 5e-5, (
        f"the five pre-ultra_long ratios deviate from 1 - ultra_long_ratio by "
        f"up to {populated.max():.6f}"
    )


#: (length, channel it must land in, channel the pre-0.9.0 gene path used).
#:
#: Every one of these is a length where the two band sets disagreed. Checking
#: the partition alone cannot see this: a fragment in the wrong channel still
#: sums correctly, so the totals and the ratios all look fine.
MISROUTED_BEFORE = [
    (100, "ultra_short", "core_short"),
    (230, "di_nucl", "mono_nucl"),
    (255, "di_nucl", "mono_nucl"),
    (300, "long", "di_nucl"),
    (390, "long", "di_nucl"),
    (500, "ultra_long", "long"),
]


@pytest.mark.parametrize("length,expected,previously", MISROUTED_BEFORE)
def test_a_known_length_lands_in_the_documented_channel(
    tmp_path, length, expected, previously
):
    """Route one length at a time and check where it came out.

    This is the assertion that actually distinguishes the old bands from the
    new ones end to end.
    """
    profile = SynthProfile(
        name=f"len{length}",
        n_fragments=300,
        # Zero spread: every fragment is exactly `length`, so a single
        # non-zero channel is unambiguous.
        length_mixture=((length, 0, 1.0),),
        seed=length,
    )
    df = _aggregate(tmp_path, profile, "gene")

    populated = [c for c in CHANNELS if df[c].sum() > 0]
    assert populated == [expected], (
        f"{length}bp landed in {populated or 'no channel'}, expected "
        f"['{expected}']"
        + (
            f"; '{previously}' is where the pre-0.9.0 gene bands put it"
            if previously in populated
            else ""
        )
    )


def test_ultra_long_is_not_folded_into_long(gene_table):
    """The specific regression: >400bp fragments used to land in `long`.

    ALL_BANDS puts a sixth of its fragments at 600bp, so a gene path that
    still folded them into `long` would leave `ultra_long` empty.
    """
    assert gene_table["ultra_long"].sum() > 0, (
        "no fragments in ultra_long despite a profile with a 600bp mode; "
        "the gene path is still folding >400bp into 'long'"
    )


def test_each_band_receives_its_own_mode(gene_table):
    """Every channel is populated, so none of the assertions above is vacuous.

    A band that never receives a fragment cannot fail a partition check, so
    without this the suite could pass on a classifier that routes everything
    into one channel.
    """
    empty = [c for c in CHANNELS if gene_table[c].sum() == 0]
    assert not empty, (
        f"channel(s) {empty} received nothing though ALL_BANDS has a mode in "
        "each; either a band boundary moved or the profile no longer spans them"
    )


def test_gene_totals_are_the_sum_of_their_regions(tmp_path):
    """The gene table is a rollup of the region table, channel by channel.

    Guards `GeneResult::merge`, where forgetting to accumulate the new sixth
    channel would lose ultra-long counts only for genes spanning more than one
    region -- invisible in any single-region fixture.
    """
    gene = _aggregate(tmp_path / "g", ALL_BANDS, "gene").set_index("gene")
    region = _aggregate(tmp_path / "r", ALL_BANDS, "region")

    assert (
        region.groupby("gene").size() > 1
    ).any(), "no gene spans multiple regions; this test cannot see a merge bug"

    rolled = region.groupby("gene")[CHANNELS + ["total"]].sum()
    shared: List[str] = sorted(set(gene.index) & set(rolled.index))
    assert shared, "gene and region tables share no genes"

    for column in CHANNELS + ["total"]:
        diff = (gene.loc[shared, column] - rolled.loc[shared, column]).abs()
        # Each region contributes its own 2dp rounding to the gene sum.
        tol = 0.005 * region.groupby("gene").size().loc[shared] + 0.005
        assert (diff <= tol).all(), (
            f"'{column}' in the gene table does not match the sum over its "
            f"regions (largest gap {diff.max():g})"
        )


def test_a_baseline_profile_still_partitions(tmp_path):
    """ALL_BANDS is deliberately unusual; the ordinary case must hold too."""
    df = _aggregate(tmp_path, BASELINE, "gene")
    diff = (df[CHANNELS].sum(axis=1) - df["total"]).abs()
    assert diff.max() <= 6 * 0.005
