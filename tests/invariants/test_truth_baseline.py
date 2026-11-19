"""Today's extractor against the simulated truth: the Step 0 baseline.

Spec §2f asks every value-changing commit for its error against simulated
truth *before and after*. This file is the "before", pinned exactly: the set of
cases where extraction deviates from the truth is asserted, not just counted,
so a change that fixes a case and one that breaks a case both fail here until
the record below is updated -- in the same commit, with the evidence the spec
asks for. Run with production's settings (``require_proper_pair=False``).

What the 0.9.x extractor emits is a BED of (chrom, start, end, gc), so that is
what can be scored today; the §2a per-fragment fields are scored once Step 1
emits them.
"""

from __future__ import annotations

import pytest

from .pipeline import extract_bam_to_bedgz
from .truth import CASES, make_truth_sample, score_bed

#: Cases where 0.9.x extraction deviates from the truth, and how.
KNOWN_DEVIATIONS = {
    # gc is 0.0 where no A/C/G/T base exists to resolve it; the spec requires
    # null. A 0.0 is a plausible value and enters the GC model as one.
    "all_n_insert": "gc_reported_unresolvable",
    "contig_absent_from_fasta": "gc_reported_unresolvable",
}

#: The BED stores gc to 4 decimals; anything above this is a real error.
GC_ROUNDING = 5e-5


@pytest.fixture(scope="module")
def scored(tmp_path_factory):
    d = tmp_path_factory.mktemp("truth_baseline")
    sample = make_truth_sample(d / "sample")
    bed = extract_bam_to_bedgz(
        sample.bam, sample.reference, d / "extract", require_proper_pair=False
    )
    return score_bed(sample.truth, bed)


def _deviations(scored):
    out = {}
    for case, r in scored.items():
        if case == "_extra_rows":
            continue
        if r["emitted_at_true_interval"] != r["expected_at_interval"]:
            out[case] = "coordinates"
        elif r["gc_reported_unresolvable"]:
            out[case] = "gc_reported_unresolvable"
        elif r["gc_error"] is not None and r["gc_error"] > GC_ROUNDING:
            out[case] = "gc_error"
    return out


def test_every_case_was_scored(scored):
    assert set(scored) - {"_extra_rows"} == {c.name for c in CASES}


def test_no_fragment_appears_at_an_interval_the_truth_does_not_have(scored):
    assert scored["_extra_rows"]["count"] == 0


def test_the_deviations_are_exactly_the_recorded_baseline(scored):
    assert _deviations(scored) == KNOWN_DEVIATIONS, (
        "extraction against simulated truth changed. If this is intended, update "
        "KNOWN_DEVIATIONS in the same commit and report the before/after (spec §2f)."
    )
