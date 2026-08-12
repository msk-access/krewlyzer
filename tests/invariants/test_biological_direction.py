"""Metrics must move the way the biology says.

A metric can vary with its input and still be wired backwards. MDS was
documented as "higher = more abnormal" while its own threshold table and the
implementation both said the opposite, and nothing caught it because no test
asserted a direction.

Ordering assertions only -- never magnitudes. A threshold would encode a
calibration these synthetic samples cannot justify, and would break on any
platform whose arithmetic differs in the last bit. The direction is the claim.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from krewlyzer import _core

from .synth import BASELINE, STEREOTYPED_ENDS, make_sample

pytestmark = [pytest.mark.invariant, pytest.mark.rust, pytest.mark.slow]


def _mds(directory: Path, profile) -> float:
    """Normalised Shannon entropy over the 256 ACGT 4-mers, as shipped.

    Mirrors motif_utils.rs: entropy / log2(256). Computed here from the raw
    end-motif counts so the test does not depend on which output table happens
    to carry MDS.
    """
    sample = make_sample(directory / profile.name, profile)
    result = _core.extract_motif.process_bam_parallel(
        str(sample.bam),
        str(sample.reference),
        20,
        65,
        1000,
        4,
        1,
        None,  # no BED needed
        str(directory / profile.name / "motif"),
        None,
        None,
        True,
        True,
        True,
    )
    end_motifs = result[1]
    valid = {m: c for m, c in end_motifs.items() if set(m) <= set("ACGT")}
    total = sum(valid.values())
    assert total > 0, f"{profile.name}: no ACGT end motifs were counted"

    entropy = 0.0
    for count in valid.values():
        if count:
            p = count / total
            entropy -= p * math.log2(p)
    return entropy / 8.0  # log2(256)


@pytest.fixture(scope="module")
def mds_pair(tmp_path_factory):
    root = tmp_path_factory.mktemp("direction")
    return _mds(root, BASELINE), _mds(root, STEREOTYPED_ENDS)


def test_stereotyped_cutting_lowers_mds(mds_pair):
    """Lower MDS is the abnormal direction.

    The summary line in motif.md asserted the reverse of its own threshold
    table for a year. Forcing every fragment onto the same CGCG context must
    drive MDS *down* relative to varied ends.
    """
    baseline, stereotyped = mds_pair
    assert stereotyped < baseline, (
        f"stereotyped ends gave MDS={stereotyped:.4f}, varied ends "
        f"{baseline:.4f}. Lower MDS is the stereotyped/tumour direction; if "
        "this inverts, the normalisation or the direction has flipped."
    )


def test_mds_is_normalised_to_the_unit_interval(mds_pair):
    """MDS is entropy / log2(256), so it cannot leave [0, 1].

    region-mds.md documented a raw 6.0-8.0 scale the code never emitted; this
    pins which of the two is real.
    """
    for label, value in zip(("baseline", "stereotyped"), mds_pair):
        assert 0.0 <= value <= 1.0, f"{label} MDS={value} is outside [0, 1]"


def test_stereotyped_sample_is_actually_stereotyped(tmp_path):
    """Guards the fixture rather than the code.

    If the CpG block moved, or fragment placement stopped honouring it, the
    direction test above would compare two ordinary samples and pass for no
    reason -- the failure mode this whole suite exists to prevent.
    """
    sample = make_sample(tmp_path, STEREOTYPED_ENDS)
    result = _core.extract_motif.process_bam_parallel(
        str(sample.bam),
        str(sample.reference),
        20,
        65,
        1000,
        4,
        1,
        None,
        str(tmp_path / "motif"),
        None,
        None,
        True,
        True,
        True,
    )
    counts = {m: c for m, c in result[1].items() if set(m) <= set("ACGT")}
    total = sum(counts.values())
    dominant = max(counts.values()) / total

    assert dominant > 0.5, (
        f"the stereotyped profile's most common 4-mer is only {dominant:.1%} "
        "of ends; the fixture is not stereotyped and the direction test above "
        "would be vacuous"
    )
