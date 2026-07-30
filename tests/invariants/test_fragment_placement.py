"""Fragment intervals must not depend on which mate is R1.

The BED writer computes ``pos() + |tlen|``, which is correct only when R1 is
the leftmost mate. For a reverse-strand R1 -- R1 rightmost, negative insert
size -- the interval is shifted right by roughly ``tlen - read_length``.

Nothing caught it because every paired fixture in this repository uses flags
99/147 (R1 forward). Flags 83/163 appear nowhere else, and the one coordinate
assertion in the suite is on a flag-99 read.

The interval feeds OCF's ±1bp end phasing, WPS fragment centres, FSC/FSD/
region-entropy overlap queries, and the GC value in BED column 4 -- so this is
upstream of most positional features. Fragment *lengths* stay correct either
way, which is why size-only features are unaffected.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import List, Tuple

import pytest

from krewlyzer import _core

from .synth import ALL_FORWARD, ALL_REVERSE, Contig, FragmentSpec, make_sample
from .synth import make_reference, write_bam

pytestmark = [pytest.mark.invariant, pytest.mark.rust]

CONTIGS = (Contig("chr1", 200_000),)


def _read_bed(path: Path) -> List[Tuple[str, int, int]]:
    with open(path, "rb") as fh:
        is_gz = fh.read(2) == b"\x1f\x8b"
    opener = gzip.open if is_gz else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, start, end = line.split("\t")[:3]
            rows.append((chrom, int(start), int(end)))
    return sorted(rows)


def _extract_bed(tmp_path: Path, sample) -> List[Tuple[str, int, int]]:
    """Run the real extractor and return the intervals it wrote."""
    bed_out = tmp_path / "out.bed"
    _core.extract_motif.process_bam_parallel(
        str(sample.bam),
        str(sample.reference),
        20,  # mapq
        65,  # min_len
        1000,  # max_len
        4,  # kmer
        1,  # threads
        str(bed_out),  # <- the writer under test; previously passed None
        None,  # motif prefix
        None,  # exclude
        None,  # target regions
        True,  # skip duplicates
        True,  # require proper pair
        True,  # silent
    )
    for candidate in (Path(str(bed_out) + ".gz"), bed_out):
        if candidate.exists():
            return _read_bed(candidate)
    raise AssertionError(f"extraction wrote no BED (looked for {bed_out}[.gz])")


def test_forward_r1_intervals_match_truth(tmp_path):
    """The case that has always worked, pinned so the fix cannot break it."""
    sample = make_sample(tmp_path / "fwd", ALL_FORWARD, CONTIGS)
    assert _extract_bed(tmp_path / "fwd", sample) == _read_bed(sample.truth_bed)


@pytest.mark.xfail(
    strict=True,
    reason="pos() + |tlen| is wrong for reverse-strand R1; the fragment is "
    "shifted right by tlen - read_length. Strict, so the fix must remove "
    "this marker rather than quietly satisfying it.",
)
def test_reverse_r1_intervals_match_truth(tmp_path):
    """The same fragments, differing only in which mate carries the R1 flag."""
    sample = make_sample(tmp_path / "rev", ALL_REVERSE, CONTIGS)
    assert _extract_bed(tmp_path / "rev", sample) == _read_bed(sample.truth_bed)


@pytest.mark.xfail(
    strict=True,
    reason="same defect, stated as the invariant it violates: a fragment's "
    "coordinates cannot depend on which mate was sequenced first.",
)
def test_interval_is_independent_of_mate_orientation(tmp_path):
    """The sharpest form of the claim.

    Two BAMs describing the *identical* physical fragment, differing only in
    the R1 flag, must yield the identical interval. Reproduced by hand on real
    code: 1:1000-1150 emitted as both 1:1000-1150 and 1:1100-1250.
    """
    contigs = (Contig("chr1", 50_000),)
    reference = make_reference(tmp_path / "ref.fa", contigs)
    frag = FragmentSpec("chr1", 1_000, 150)

    intervals = {}
    for label, reverse in (("forward", False), ("reverse", True)):
        directory = tmp_path / label
        directory.mkdir()
        bam = write_bam(
            directory / "s.bam",
            [FragmentSpec(frag.chrom, frag.start, frag.length, r1_reverse=reverse)],
            contigs,
            reference,
        )
        sample = type("S", (), {"bam": bam, "reference": reference})()
        intervals[label] = _extract_bed(directory, sample)

    assert intervals["forward"] == intervals["reverse"], (
        f"the same fragment produced {intervals['forward']} with a forward R1 "
        f"and {intervals['reverse']} with a reverse R1"
    )


def test_fragment_lengths_survive_either_orientation(tmp_path):
    """Lengths come from |tlen| and are correct regardless of orientation.

    Worth pinning explicitly: it bounds the blast radius of the coordinate bug
    to positional features, and stops a future fix from disturbing the size
    distributions that FSD/FSR/FSC depend on.
    """
    lengths = {}
    for label, profile in (("fwd", ALL_FORWARD), ("rev", ALL_REVERSE)):
        sample = make_sample(tmp_path / label, profile, CONTIGS)
        rows = _extract_bed(tmp_path / label, sample)
        lengths[label] = sorted(end - start for _, start, end in rows)

    assert (
        lengths["fwd"] == lengths["rev"]
    ), "fragment lengths must not depend on mate orientation"
