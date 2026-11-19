"""The simulated truth must agree with the BAM it describes.

``truth.py`` computes each case's expected record from the case specification.
Here the BAM is read back independently -- pysam alignments, the FASTA, the
targets BED -- and every truth column that the BAM determines is checked
against it. A truth table that disagrees with its own BAM would turn every
"error before and after" figure (spec §2f) into a measurement of the
generator, so this has to pass before any extractor is scored against it.
"""

from __future__ import annotations

from collections import defaultdict

import pysam
import pytest

from .truth import CASES, CHUNK, N_BLOCK, make_truth_sample, revcomp


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_truth_sample(tmp_path_factory.mktemp("truth"))


@pytest.fixture(scope="module")
def reads(sample):
    primaries, supplementary = defaultdict(list), defaultdict(list)
    with pysam.AlignmentFile(str(sample.bam)) as bam:
        for r in bam.fetch(until_eof=True):
            (supplementary if r.is_supplementary else primaries)[r.query_name].append(r)
    return primaries, supplementary


def _rows(sample):
    return {t["case"]: t for t in sample.truth}


def _aligned_bases(r):
    seq = r.query_sequence
    return "".join(seq[q] for q, _ in r.get_aligned_pairs(matches_only=True))


def _soft(r, side):
    cig = [c for c in r.cigartuples if c[0] != 5]  # hard clips sit outside soft clips
    op, n = cig[0] if side == "left" else cig[-1]
    return n if op == 4 else 0


def test_every_spec_case_is_present_once(sample):
    names = [c.name for c in CASES]
    assert len(names) == len(set(names))
    assert {t["case"] for t in sample.truth} == set(names)


def test_coordinates_spans_and_indels_match_the_alignments(sample, reads):
    primaries, _ = reads
    for t in sample.truth:
        recs = sorted(primaries[t["qname"]], key=lambda r: r.reference_start)
        assert len(recs) == (2 if t["mate_seen"] else 1), t["case"]
        left = recs[0]
        assert left.reference_start == t["start"], t["case"]
        assert left.reference_length == t["left_span"], t["case"]
        assert (
            abs(left.template_length) == t["tlen_reported"] == t["end"] - t["start"]
        ), t["case"]
        ins = sum(n for op, n in left.cigartuples if op == 1)
        dele = sum(n for op, n in left.cigartuples if op == 2)
        assert (ins, dele) == (t["ins_bp_left"], t["del_bp_left"]), t["case"]
        assert left.is_read1 == t["r1_is_left"], t["case"]
        assert left.get_tag("NM") == t["nm_left"], t["case"]
        if t["mate_seen"]:
            right = recs[1]
            assert right.reference_end == t["end"], t["case"]
            assert right.reference_length == t["right_span"], t["case"]
            assert right.is_reverse and not left.is_reverse, t["case"]
            ins = sum(n for op, n in right.cigartuples if op == 1)
            dele = sum(n for op, n in right.cigartuples if op == 2)
            assert (ins, dele) == (t["ins_bp_right"], t["del_bp_right"]), t["case"]
            assert right.get_tag("NM") == t["nm_right"], t["case"]
            assert (
                left.get_tag("MC") == right.cigarstring
                and right.get_tag("MC") == left.cigarstring
            )
        else:
            # Coordinates come from TLEN when the mate is absent.
            assert left.reference_start + left.template_length == t["end"], t["case"]


def test_clips_are_on_the_side_the_truth_says(sample, reads):
    primaries, _ = reads
    for t in sample.truth:
        recs = sorted(primaries[t["qname"]], key=lambda r: r.reference_start)
        left = recs[0]
        # Forward left mate: 5' (fragment end) is its leading side.
        assert _soft(left, "left") == t["clip_end_left"], t["case"]
        assert _soft(left, "right") == t["clip_inner_left"], t["case"]
        hard = any(op == 5 for r in recs for op, _ in r.cigartuples)
        assert hard == t["hard_clipped"], t["case"]
        if t["mate_seen"]:
            right = recs[1]
            # Reverse right mate: 5' is its trailing side in reference order.
            assert _soft(right, "right") == t["clip_end_right"], t["case"]
            assert _soft(right, "left") == t["clip_inner_right"], t["case"]


def test_end_kmers_and_mismatches_match_reads_and_reference(sample, reads):
    primaries, _ = reads
    fasta = pysam.FastaFile(str(sample.reference))
    for t in sample.truth:
        if not t["reference_available"]:
            continue
        recs = sorted(primaries[t["qname"]], key=lambda r: r.reference_start)
        lk = _aligned_bases(recs[0])[:8]
        assert lk == t["kmer_left_ref"], t["case"]
        ref_l = fasta.fetch(t_chrom(t), t["start"], t["start"] + 8)
        if "N" not in lk:
            assert sum(a != b for a, b in zip(lk, ref_l)) == t["end_mismatch_left"], t[
                "case"
            ]
            assert t["read_kmer_left"] == lk
        else:
            assert t["read_kmer_left"] is None and t["end_mismatch_left"] is None
        if t["mate_seen"]:
            rk = _aligned_bases(recs[1])[-8:]
            assert rk == t["kmer_right_ref"], t["case"]
            if "N" not in rk:
                ref_r = fasta.fetch(t_chrom(t), t["end"] - 8, t["end"])
                assert (
                    sum(a != b for a, b in zip(rk, ref_r)) == t["end_mismatch_right"]
                ), t["case"]
                assert t["read_kmer_right"] == revcomp(rk)
    assert _rows(sample)["n_in_end_kmer_left"]["read_kmer_left"] is None
    assert _rows(sample)["mismatch_in_end_kmer_right"]["end_mismatch_right"] == 1


def t_chrom(t):
    return "chr" + t["chrom"]


def test_quality_tags_and_families_match(sample, reads):
    primaries, supplementary = reads
    for t in sample.truth:
        recs = sorted(primaries[t["qname"]], key=lambda r: r.reference_start)
        assert recs[0].mapping_quality == t["mapq_left"], t["case"]
        has_sa = any(r.has_tag("SA") for r in recs) or bool(
            supplementary.get(t["qname"])
        )
        assert has_sa == t["has_supplementary"], t["case"]
        if t["family_source"] == "fgbio":
            r = recs[0]
            got = (r.get_tag("aD"), r.get_tag("bD"), r.get_tag("cD"), r.get_tag("cM"))
            assert got == (
                t["family_pos"],
                t["family_neg"],
                t["consensus_depth_max"],
                t["consensus_depth_min"],
            )
        elif t["family_source"] == "marianas":
            w = t["qname"].split(":")
            assert (int(w[4]), int(w[5])) == (t["family_pos"], t["family_neg"]), t[
                "case"
            ]
            assert int(w[3]) - 1 == t["name_start_left"] == recs[0].reference_start, t[
                "case"
            ]
            assert int(w[7]) - 1 == t["name_start_right"] == recs[1].reference_start, t[
                "case"
            ]
        if t["mate_seen"]:
            assert recs[1].mapping_quality == t["mapq_right"], t["case"]
            assert (
                recs[0].get_tag("MQ") == t["mapq_right"]
                and recs[1].get_tag("MQ") == t["mapq_left"]
            )


def test_gc_and_target_relation_recomputed_independently(sample):
    fasta = pysam.FastaFile(str(sample.reference))
    targets = [line.split("\t") for line in sample.targets.read_text().splitlines()]
    for t in sample.truth:
        if t["reference_available"]:
            s = fasta.fetch(t_chrom(t), t["start"], t["end"])
            acgt = sum(s.count(b) for b in "ACGT")
            want = None if acgt == 0 else (s.count("G") + s.count("C")) / acgt
            assert t["gc"] == want, t["case"]
        else:
            assert t["gc"] is None and t_chrom(t) not in fasta.references
        same = [(int(a), int(b)) for c, a, b in targets if c == t_chrom(t)]
        overlap = sum(max(0, min(t["end"], b) - max(t["start"], a)) for a, b in same)
        assert overlap == t["target_overlap_bp"], t["case"]
    rows = _rows(sample)
    assert (
        rows["target_abut"]["target_overlap_bp"] == 0
        and rows["target_abut"]["target_distance"] == 0
    )
    assert rows["target_overlap"]["target_overlap_bp"] > 0
    assert rows["target_far"]["target_distance"] > 10_000


def test_the_reference_cases_are_what_they_claim(sample, reads):
    primaries, _ = reads
    fasta = pysam.FastaFile(str(sample.reference))
    rows = _rows(sample)
    n = rows["all_n_insert"]
    assert set(fasta.fetch("chr1", n["start"], n["end"])) == {"N"}
    assert N_BLOCK[0] <= n["start"] and n["end"] <= N_BLOCK[1]
    assert n["gc"] is None
    chunk = sorted(
        primaries[rows["mate_in_next_chunk"]["qname"]], key=lambda r: r.reference_start
    )
    assert chunk[0].reference_start < CHUNK <= chunk[1].reference_start
    a, b = rows["identical_coords_a"], rows["identical_coords_b"]
    assert (a["start"], a["end"]) == (b["start"], b["end"]) and a["qname"] != b["qname"]
