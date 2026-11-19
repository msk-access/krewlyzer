"""Simulated truth for the fragment model (spec §2f, test 1).

Writes a BAM in which every fragment is one named case from the spec's §4 list
-- clips, indels, a missing mate, a mate across a chunk boundary, low mate MAPQ,
N and mismatches in the end 8-mers, supplementary alignments, fgbio tags and
Marianas read names, target relations, identical coordinates -- together with
the record extraction *should* produce for it, in the §2a column names.

The truth is computed from the case specification, never read back from the
BAM; ``test_truth_generator.py`` reads the BAM independently with pysam and
checks the two agree, because a truth table that disagrees with its own BAM
makes every error figure meaningless. ``score_bed`` then compares an
extraction against the truth case by case: run on today's extractor it is the
Step 0 baseline every value-changing commit must improve on.

Conventions the spec leaves open, stated here so each is decided once and
visibly (flagged for the step that has to pin it):

* ``read_kmer_right`` is the reverse complement of the right mate's last 8
  aligned bases -- read 5'->3' on the reverse strand, as end-motif analyses
  conventionally report it. ``kmer_*_ref`` keep the reference orientation.
  (Step 1.)
* Marianas names encode each mate's start **1-based**; ``name_start_*`` are the
  0-based starts the names imply. Check against real BAMs in Step 1.
* ``target_distance`` is the gap in bases to the nearest target, so an
  abutting fragment has distance 0 and overlap 0. (Step 5.)
* ``gc`` is G+C over A/C/G/T bases of the reference insert, None when the
  insert has no A/C/G/T base -- "null when unresolved, never 0.0". (Step 3.)

``make_reference`` in ``synth.py`` is not reused: these cases need an N block,
a contig past the extractor's 10 Mb chunk, and a contig absent from the FASTA,
and changing that function's random sequence would silently change the
fixtures of every test built on it.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pysam

from .synth import R1_FWD, R1_REV, R2_FWD, R2_REV, Contig, _header

READ_LEN = 50
QUAL = 40  # 'I'

#: The extractor fetches in 10 Mb chunks (extract_motif.rs); CHUNK_CONTIG is
#: long enough to put a fragment's mates on either side of that boundary.
CHUNK = 10_000_000
MAIN = Contig("chr1", 200_000)
CHUNK_CONTIG = Contig("chr2", CHUNK + 1_000)
ABSENT = Contig("chrZ", 10_000)  # in the BAM header, not in the FASTA
CONTIGS = (MAIN, CHUNK_CONTIG, ABSENT)

N_BLOCK = (150_000, 150_600)  # reference N run on chr1, for the all-N insert
TARGETS = (("chr1", 20_000, 20_200), ("chr1", 40_000, 40_200), ("chr1", 60_000, 60_200))

_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _other_base(b: str) -> str:
    return "ACGT"[("ACGT".index(b) + 1) % 4] if b in "ACGT" else "A"


@dataclass(frozen=True)
class Mate:
    """One mate, in terms of the fragment rather than the CIGAR.

    ``clip5`` is the soft clip at the mate's 5' end -- the **fragment end** --
    and ``clip3`` at its 3' end, the inner side. ``ins``/``dele`` are
    ``(offset, length)`` from the mate's leftmost aligned base, reference
    orientation. ``kmer_edit`` is ``(kind, k)``: ``"N"`` or ``"mm"`` at the
    k-th aligned base counted inward from the fragment end.
    """

    aligned: int = READ_LEN  # aligned query bases (M)
    clip5: int = 0
    clip3: int = 0
    hard5: int = 0
    ins: Optional[Tuple[int, int]] = None
    dele: Optional[Tuple[int, int]] = None
    mapq: int = 60
    kmer_edit: Optional[Tuple[str, int]] = None
    present: bool = True
    sa: bool = False  # carries an SA tag, and a supplementary record is written

    @property
    def ref_span(self) -> int:
        return self.aligned + (self.dele[1] if self.dele else 0)


@dataclass(frozen=True)
class Case:
    name: str
    chrom: str
    start: int  # outermost aligned base of the left mate, 0-based
    length: int  # end - start
    r1_left: bool = True
    left: Mate = field(default_factory=Mate)
    right: Mate = field(default_factory=Mate)
    family: Optional[Tuple] = (
        None  # ("marianas", umi, ps1, ns1, ps2, ns2) | ("fgbio", umi, aD, bD, cD, cM, cE)
    )
    qname: Optional[str] = None

    @property
    def end(self) -> int:
        return self.start + self.length


def _cases() -> List[Case]:
    """The spec's §4 cases. Positions are spread so no two interact unless
    the case is about interaction."""
    c: List[Case] = []
    pos = iter(range(2_000, 140_000, 2_500))

    def at(name: str, length: int = 167, **kw) -> Case:
        return Case(name, MAIN.name, next(pos), length, **kw)

    c += [
        at("clean_r1_left"),
        at("clean_r1_right", r1_left=False),
        at("mate_missing", right=Mate(present=False)),
        at("softclip_end_left", left=Mate(aligned=45, clip5=5)),
        at("softclip_end_right", right=Mate(aligned=45, clip5=5)),
        at("softclip_inner_left", left=Mate(aligned=44, clip3=6)),
        at("hardclip_left", left=Mate(aligned=46, hard5=4)),
        at("insertion_left", left=Mate(aligned=47, ins=(25, 3))),
        at("deletion_right", right=Mate(aligned=50, dele=(20, 4))),
        at("mates_overlap", length=80),
        at("mate_low_mapq", right=Mate(mapq=5)),
        at("n_in_end_kmer_left", left=Mate(kmer_edit=("N", 2))),
        at("mismatch_in_end_kmer_right", right=Mate(kmer_edit=("mm", 1))),
        at("long_fragment", length=450),
        at("supplementary_left", left=Mate(sa=True)),
        at("fgbio_tagged", family=("fgbio", "ACGT-TTGA", 3, 2, 5, 2, 0.01)),
        at("marianas_duplex", family=("marianas", "ACG+TTA", 2, 1, 1, 2)),
        at("marianas_simplex", family=("marianas", "CCA+GTT", 3, 0, 3, 0)),
        at("marianas_singleton", family=("marianas", "GAT+CAG", 1, 0, 1, 0)),
    ]
    # Target relation: overlapping, abutting (end == target start), far.
    c += [
        Case("target_overlap", MAIN.name, 20_100, 167),
        Case("target_abut", MAIN.name, 40_000 - 167, 167),
        Case("target_far", MAIN.name, 100_500, 167),
    ]
    # Two fragments overlapping each other; two distinct ones at identical coordinates.
    c += [
        Case("overlapping_pair_a", MAIN.name, 142_000, 167),
        Case("overlapping_pair_b", MAIN.name, 142_080, 167),
        Case(
            "identical_coords_a",
            MAIN.name,
            145_000,
            167,
            family=("fgbio", "AAAA-CCCC", 2, 2, 4, 2, 0.0),
        ),
        Case(
            "identical_coords_b",
            MAIN.name,
            145_000,
            167,
            family=("fgbio", "GGGG-TTTT", 2, 2, 4, 2, 0.0),
        ),
    ]
    # Reference trouble: an insert of N only; a contig the FASTA lacks.
    c += [
        Case("all_n_insert", MAIN.name, N_BLOCK[0] + 100, 167),
        Case("contig_absent_from_fasta", ABSENT.name, 2_000, 167),
    ]
    # A fragment whose mates straddle the extractor's chunk boundary.
    c.append(Case("mate_in_next_chunk", CHUNK_CONTIG.name, CHUNK - 200, 400))
    return c


CASES: Tuple[Case, ...] = tuple(_cases())


# --------------------------------------------------------------------------
# Reference
# --------------------------------------------------------------------------


def make_truth_reference(path: Path, seed: int = 11) -> Dict[str, str]:
    """FASTA for chr1 and chr2 (not chrZ). Returns the sequences."""
    rng = random.Random(seed)
    seqs: Dict[str, str] = {}
    for contig in (MAIN, CHUNK_CONTIG):
        s = rng.choices("ACGT", k=contig.length)
        if contig is MAIN:
            s[N_BLOCK[0] : N_BLOCK[1]] = "N" * (N_BLOCK[1] - N_BLOCK[0])
        seqs[contig.name] = "".join(s)
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return seqs


def _ref(seqs: Dict[str, str], chrom: str, a: int, b: int) -> str:
    # The absent contig has no sequence; a deterministic stand-in keeps the
    # reads well-formed, and the truth records that no reference exists.
    if chrom not in seqs:
        return ("ACGT" * ((b - a) // 4 + 1))[: b - a]
    return seqs[chrom][a:b]


# --------------------------------------------------------------------------
# One mate: CIGAR, query sequence, and what the fragment end looks like
# --------------------------------------------------------------------------


@dataclass
class _Built:
    ref_start: int
    cigar: List[Tuple[int, int]]
    seq: str
    aligned_ref_orient: (
        str  # aligned query bases, reference orientation, no clips/insertions
    )
    nm: int


def _build(
    m: Mate, ref_start: int, reverse: bool, seqs: Dict[str, str], chrom: str
) -> _Built:
    """Lay one mate on the reference. Reference order: for a forward mate the
    5' side is on the left, for a reverse mate on the right."""
    ref = _ref(seqs, chrom, ref_start, ref_start + m.ref_span)
    q, cig, aligned, nm, r = [], [], [], 0, 0
    q_aligned_done = 0

    def match(n: int) -> None:
        nonlocal r, q_aligned_done
        part = list(ref[r : r + n])
        q.extend(part)
        aligned.extend(part)
        cig.append((0, n))
        r += n
        q_aligned_done += n

    events = []
    if m.ins:
        events.append((m.ins[0], "I", m.ins[1]))
    if m.dele:
        events.append((m.dele[0], "D", m.dele[1]))
    events.sort()
    for off, kind, n in events:
        match(off - q_aligned_done)
        if kind == "I":
            q.extend("A" * n)
            cig.append((1, n))
        else:
            cig.append((2, n))
            r += n
        nm += n
    match(m.aligned - q_aligned_done)

    # The end k-mer edit, k bases inward from the fragment end.
    if m.kmer_edit:
        kind, k = m.kmer_edit
        idx_aligned = (len(aligned) - 1 - k) if reverse else k
        # Map aligned index to query index: no insertion precedes the first
        # or follows the last 8 aligned bases in any case defined above.
        idx_q = idx_aligned if not reverse else len(q) - 1 - k
        new = "N" if kind == "N" else _other_base(aligned[idx_aligned])
        aligned[idx_aligned] = new
        q[idx_q] = new
        nm += 1

    left_clip, right_clip = (m.clip3, m.clip5) if reverse else (m.clip5, m.clip3)
    left_hard, right_hard = (0, m.hard5) if reverse else (m.hard5, 0)
    cigar = (
        [(5, left_hard)] * bool(left_hard)
        + [(4, left_clip)] * bool(left_clip)
        + cig
        + [(4, right_clip)] * bool(right_clip)
        + [(5, right_hard)] * bool(right_hard)
    )
    seq = "T" * left_clip + "".join(q) + "T" * right_clip
    return _Built(ref_start, cigar, seq, "".join(aligned), nm)


# --------------------------------------------------------------------------
# Truth
# --------------------------------------------------------------------------


def _target_relation(chrom: str, start: int, end: int) -> Tuple[int, int]:
    overlap, distance = 0, None
    for c, a, b in TARGETS:
        if c != chrom:
            continue
        overlap += max(0, min(end, b) - max(start, a))
        gap = max(0, a - end, start - b)
        distance = gap if distance is None else min(distance, gap)
    return overlap, (0 if overlap else distance)


def _gc(seqs: Dict[str, str], chrom: str, start: int, end: int) -> Optional[float]:
    if chrom not in seqs:
        return None
    s = seqs[chrom][start:end]
    acgt = sum(s.count(b) for b in "ACGT")
    return None if acgt == 0 else (s.count("G") + s.count("C")) / acgt


def _cigar_str(cigar: Sequence[Tuple[int, int]]) -> str:
    return "".join(f"{n}{'MIDNSHP=X'[op]}" for op, n in cigar)


def _kmer(bases: str) -> Optional[str]:
    return None if "N" in bases else bases


def _truth(case: Case, left: _Built, right: _Built, seqs: Dict[str, str]) -> Dict:
    lk_ref = left.aligned_ref_orient[:8]
    rk_ref = right.aligned_ref_orient[-8:]
    ref_l = _ref(seqs, case.chrom, case.start, case.start + 8)
    ref_r = _ref(seqs, case.chrom, case.end - 8, case.end)
    overlap, distance = _target_relation(case.chrom, case.start, case.end)
    fam = case.family or ()
    seen = case.right.present  # every case that drops a mate drops the right one
    row = {
        "case": case.name,
        "chrom": case.chrom.removeprefix("chr"),
        "start": case.start,
        "end": case.end,
        "left_span": case.left.ref_span,
        "right_span": case.right.ref_span if seen else None,
        "ins_bp_left": case.left.ins[1] if case.left.ins else 0,
        "del_bp_left": case.left.dele[1] if case.left.dele else 0,
        "ins_bp_right": (case.right.ins[1] if case.right.ins else 0) if seen else None,
        "del_bp_right": (
            (case.right.dele[1] if case.right.dele else 0) if seen else None
        ),
        "tlen_reported": case.length,
        "mate_seen": case.left.present and seen,
        "r1_is_left": case.r1_left,
        "left_reverse": False,
        "right_reverse": True if seen else None,
        "proper_pair": True,
        "has_supplementary": case.left.sa or case.right.sa,
        "mapq_left": case.left.mapq,
        "mapq_right": case.right.mapq if seen else None,
        "nm_left": left.nm,
        "nm_right": right.nm if seen else None,
        "mean_bq_left": QUAL,
        "mean_bq_right": QUAL if seen else None,
        "clip_end_left": case.left.clip5,
        "clip_end_right": case.right.clip5 if seen else None,
        "clip_inner_left": case.left.clip3,
        "clip_inner_right": case.right.clip3 if seen else None,
        "hard_clipped": bool(case.left.hard5 or case.right.hard5),
        "kmer_left_ref": lk_ref,
        "kmer_right_ref": rk_ref if seen else None,
        "read_kmer_left": _kmer(lk_ref),
        "read_kmer_right": _kmer(revcomp(rk_ref)) if seen else None,
        "end_mismatch_left": (
            None if "N" in lk_ref else sum(a != b for a, b in zip(lk_ref, ref_l))
        ),
        "end_mismatch_right": (
            None
            if "N" in rk_ref or not seen
            else sum(a != b for a, b in zip(rk_ref, ref_r))
        ),
        "family_source": fam[0] if fam else None,
        "family_pos": fam[2] if fam else None,
        "family_neg": fam[3] if fam else None,
        "consensus_depth_max": fam[4] if fam and fam[0] == "fgbio" else None,
        "consensus_depth_min": fam[5] if fam and fam[0] == "fgbio" else None,
        "consensus_error": fam[6] if fam and fam[0] == "fgbio" else None,
        "name_start_left": case.start if fam and fam[0] == "marianas" else None,
        "name_start_right": (
            (case.end - case.right.ref_span) if fam and fam[0] == "marianas" else None
        ),
        "gc": _gc(seqs, case.chrom, case.start, case.end),
        "target_overlap_bp": overlap,
        "target_distance": distance,
        "reference_available": case.chrom in seqs,
    }
    return row


# --------------------------------------------------------------------------
# BAM
# --------------------------------------------------------------------------


def _qname(case: Case, i: int) -> str:
    if case.qname:
        return case.qname
    if case.family and case.family[0] == "marianas":
        _, umi, ps1, ns1, ps2, ns2 = case.family
        right_start = case.end - case.right.ref_span
        # Marianas DuplexUMIBamToCollapsedFastqSecondPass.makeReadName layout;
        # positions 1-based (convention, see module docstring).
        return f"Marianas:{umi}:{case.chrom}:{case.start + 1}:{ps1}:{ns1}:{case.chrom}:{right_start + 1}:{ps2}:{ns2}"
    return f"truth{i}_{case.name}"


@dataclass
class TruthSample:
    directory: Path
    reference: Path
    bam: Path
    targets: Path
    truth: List[Dict]


def make_truth_sample(directory: Path, cases: Sequence[Case] = CASES) -> TruthSample:
    directory.mkdir(parents=True, exist_ok=True)
    reference = directory / "truth_ref.fa"
    seqs = make_truth_reference(reference)
    targets = directory / "truth_targets.bed"
    targets.write_text("".join(f"{c}\t{a}\t{b}\n" for c, a, b in TARGETS))

    tid = {c.name: i for i, c in enumerate(CONTIGS)}
    unsorted = directory / "truth.unsorted.bam"
    rows: List[Dict] = []
    with pysam.AlignmentFile(str(unsorted), "wb", header=_header(CONTIGS)) as bam:
        for i, case in enumerate(cases):
            left = _build(case.left, case.start, False, seqs, case.chrom)
            right = _build(
                case.right, case.end - case.right.ref_span, True, seqs, case.chrom
            )
            name = _qname(case, i)
            rows.append({**_truth(case, left, right, seqs), "qname": name})
            mates = [(left, case.left, False), (right, case.right, True)]
            for built, spec, is_right in mates:
                if not spec.present:
                    continue
                other = right if not is_right else left
                if is_right:
                    flag = R2_REV if case.r1_left else R1_REV
                else:
                    flag = R1_FWD if case.r1_left else R2_FWD
                # A missing mate keeps its partner's flags and TLEN: the spec's case
                # is a mate absent from the BAM or filtered out, not one unmapped.
                rec = pysam.AlignedSegment()
                rec.query_name = name
                rec.flag = flag
                rec.reference_id = tid[case.chrom]
                rec.reference_start = built.ref_start
                rec.mapping_quality = spec.mapq
                rec.cigartuples = built.cigar
                rec.query_sequence = built.seq
                rec.query_qualities = pysam.qualitystring_to_array(
                    chr(QUAL + 33) * len(built.seq)
                )
                rec.next_reference_id = tid[case.chrom]
                rec.next_reference_start = other.ref_start
                rec.template_length = case.length if not is_right else -case.length
                rec.set_tag("NM", built.nm)
                rec.set_tag("AS", len(built.aligned_ref_orient) - 4 * built.nm)
                other_spec = case.right if not is_right else case.left
                if other_spec.present:
                    rec.set_tag("MQ", other_spec.mapq)
                    rec.set_tag("MC", _cigar_str(other.cigar))
                if case.family and case.family[0] == "fgbio":
                    _, umi, aD, bD, cD, cM, cE = case.family
                    rec.set_tag("RX", umi)
                    rec.set_tag("aD", aD)
                    rec.set_tag("bD", bD)
                    rec.set_tag("cD", cD)
                    rec.set_tag("cM", cM)
                    rec.set_tag("cE", float(cE))
                if spec.sa:
                    rec.set_tag(
                        "SA", f"{case.chrom},{built.ref_start + 50_001},+,30M20S,60,0;"
                    )
                bam.write(rec)
                if spec.sa:
                    sup = pysam.AlignedSegment()
                    sup.query_name = name
                    sup.flag = (flag | 0x800) & ~0x2
                    sup.reference_id = tid[case.chrom]
                    sup.reference_start = built.ref_start + 50_000
                    sup.mapping_quality = 60
                    sup.cigartuples = [(0, 30), (5, 20)]
                    sup.query_sequence = built.seq[:30]
                    sup.query_qualities = pysam.qualitystring_to_array(
                        chr(QUAL + 33) * 30
                    )
                    sup.next_reference_id = tid[case.chrom]
                    sup.next_reference_start = other.ref_start
                    sup.template_length = 0
                    bam.write(sup)

    bam_path = directory / "truth.bam"
    pysam.sort("-o", str(bam_path), str(unsorted))
    pysam.index(str(bam_path))
    unsorted.unlink(missing_ok=True)
    return TruthSample(directory, reference, bam_path, targets, rows)


# --------------------------------------------------------------------------
# Scoring an extraction against the truth
# --------------------------------------------------------------------------


def read_bed(path: Path) -> List[Tuple[str, int, int, Optional[float]]]:
    opener = gzip.open if path.read_bytes()[:2] == b"\x1f\x8b" else open
    out = []
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gc = float(f[3]) if len(f) > 3 and f[3] not in ("", "NA", "nan") else None
            out.append((f[0].removeprefix("chr"), int(f[1]), int(f[2]), gc))
    return out


def score_bed(truth: Sequence[Dict], bed: Path) -> Dict[str, Dict]:
    """Per case: how many BED rows match the true interval exactly, and the
    gc error when one does. Identical-coordinate cases share an interval, so
    their expected count is the number of cases at it."""
    rows = read_bed(bed)
    by_interval: Dict[Tuple[str, int, int], List[Optional[float]]] = {}
    for c, s, e, gc in rows:
        by_interval.setdefault((c, s, e), []).append(gc)
    expected: Dict[Tuple[str, int, int], int] = {}
    for t in truth:
        k = (t["chrom"], t["start"], t["end"])
        expected[k] = expected.get(k, 0) + 1
    out: Dict[str, Dict] = {}
    for t in truth:
        k = (t["chrom"], t["start"], t["end"])
        got = by_interval.get(k, [])
        gc = got[0] if got else None
        out[t["case"]] = {
            "emitted_at_true_interval": len(got),
            "expected_at_interval": expected[k],
            "gc_error": (None if gc is None or t["gc"] is None else abs(gc - t["gc"])),
            "gc_emitted": gc,
            # A value where none is resolvable: the spec says null, never 0.0.
            "gc_reported_unresolvable": bool(got)
            and t["gc"] is None
            and gc is not None,
        }
    matched = {k for k in by_interval if k in expected}
    out["_extra_rows"] = {
        "count": sum(len(v) for k, v in by_interval.items() if k not in matched)
    }
    return out


__all__ = [
    "CASES",
    "Case",
    "Mate",
    "TruthSample",
    "make_truth_sample",
    "score_bed",
    "read_bed",
    "revcomp",
]
