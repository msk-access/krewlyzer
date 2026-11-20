"""Contig-name mismatches must be loud (redesign item 0).

Extract matched contigs by exact name. A BAM naming `chr1` read with a FASTA naming `1` gave
every fragment GC 0.0 -- no length bin could be fitted, nothing was corrected, nothing said so.
A target BED naming `1` against a `chr1` BAM put every fragment off-target, silently, while the
feature code (which strips `chr`) still split on/off-target: the two disagreed on what "on
target" meant.
"""

import pytest

from krewlyzer import _core
from krewlyzer.core.gc_assets import (
    CONTIG_MISMATCH_FAIL_SHARE,
    CONTIG_MISMATCH_WARN_SHARE,
    check_reference_contigs,
)

from .synth import make_reference, make_sample

pytestmark = [pytest.mark.invariant, pytest.mark.rust]


def _extract(sample, reference, out, targets=None):
    return _core.extract_motif.process_bam_parallel(
        str(sample.bam),
        str(reference),
        20,
        65,
        1000,
        4,
        1,
        str(out),
        None,
        None,
        str(targets) if targets else None,
        True,
        True,
        True,
    )


def _renamed_reference(sample, directory):
    """The same sequences under `1`, `2` instead of `chr1`, `chr2`."""
    renamed = [type(c)(c.name.removeprefix("chr"), c.length) for c in sample.contigs]
    return make_reference(directory / "nochr.fa", renamed)


def test_fragments_on_contigs_missing_from_the_reference_are_counted(tmp_path):
    sample = make_sample(tmp_path / "s")
    matched = _extract(sample, sample.reference, tmp_path / "a.bed")
    mismatched = _extract(
        sample, _renamed_reference(sample, tmp_path), tmp_path / "b.bed"
    )
    assert matched[0] > 0 and matched[7]["fragments_on_missing_contig"] == 0
    assert mismatched[7]["fragments_on_missing_contig"] == mismatched[0]


def test_a_gross_mismatch_fails_and_a_small_one_warns(caplog):
    assert check_reference_contigs(0, 1000) == 0.0
    share = check_reference_contigs(20, 1000)  # 2%: e.g. decoys absent from the FASTA
    assert CONTIG_MISMATCH_WARN_SHARE < share < CONTIG_MISMATCH_FAIL_SHARE
    assert "missing from the reference" in caplog.text
    with pytest.raises(RuntimeError, match="chr"):
        check_reference_contigs(990, 1000)


def test_target_names_match_with_or_without_chr(tmp_path):
    sample = make_sample(tmp_path / "s")
    targets = tmp_path / "targets.bed"
    first = sample.contigs[0]
    targets.write_text(f"{first.name.removeprefix('chr')}\t0\t{first.length}\tT1\n")
    result = _extract(sample, sample.reference, tmp_path / "t.bed", targets=targets)
    on_target_observations = sum(result[6].values())
    assert on_target_observations > 0, "a `1`-named target BED must match a `chr1` BAM"
