"""Synthetic cfDNA inputs with known ground truth.

Builds a reference FASTA, paired-end BAMs whose fragment lengths **and mate
orientation** are exactly controlled, and the small BED assets the pipeline
needs -- so the invariant suite runs offline, with no Git LFS assets.

Two properties the previous generator lacked, both of which hid real bugs:

**Orientation.** It emitted R1 forward every time (flags 99/147). Flags 83/163
-- R1 on the reverse strand -- appear nowhere in this repository's tests. The
BED writer computes ``pos() + |tlen|``, which is only correct when R1 is the
leftmost mate, so a whole branch of the coordinate logic was unreachable by any
test. ``reverse_r1_fraction`` defaults to 0.5 for that reason.

**The BED writer itself.** It passed ``None`` for the BED path, so the fragment
writer never ran at all. :func:`write_truth_bed` gives the expected intervals so
a test can diff against what extraction actually emits.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pysam

# Flag combinations for a proper FR pair. Named because the numbers are
# opaque and getting them wrong silently produces the orientation you were
# trying to avoid testing.
R1_FWD, R2_REV = 99, 147  # R1 leftmost, forward
R1_REV, R2_FWD = 83, 163  # R1 rightmost, reverse -- the untested case

DEFAULT_READ_LEN = 50

# The CG-repeat block make_reference lays down on the first contig.
CPG_BLOCK_START = 100_000
CPG_BLOCK_LEN = 20_000


@dataclass(frozen=True)
class FragmentSpec:
    """One fragment, with the orientation of its R1 pinned."""

    chrom: str
    start: int
    length: int
    r1_reverse: bool = False
    xm: Optional[str] = None  # bismark methylation tag, for UXM

    @property
    def end(self) -> int:
        return self.start + self.length


@dataclass(frozen=True)
class Contig:
    name: str
    length: int


@dataclass(frozen=True)
class SynthProfile:
    """A recipe for one synthetic sample.

    Profiles exist so tests can assert that two *different* inputs produce two
    different outputs -- the check that would have caught the NRL degeneracy.
    """

    name: str
    n_fragments: int = 400
    # (mean, sd, weight) mixture of fragment lengths
    length_mixture: Sequence[Tuple[int, int, float]] = ((167, 25, 1.0),)
    reverse_r1_fraction: float = 0.5
    stereotyped_ends: bool = False
    """Place every fragment on the CpG block so all 5' end motifs are CGCG.

    Gives the motif-diversity direction test a contrast: a sample whose
    cutting is maximally stereotyped must score lower on MDS than one with
    varied ends, and MDS scored the wrong way round in the docs for a year.
    """
    seed: int = 0
    read_len: int = DEFAULT_READ_LEN


BASELINE = SynthProfile("baseline", seed=1)
SHORT_SHIFTED = SynthProfile("short_shifted", length_mixture=((120, 20, 1.0),), seed=2)
ULTRA_LONG_HEAVY = SynthProfile(
    "ultra_long_heavy",
    length_mixture=((167, 25, 0.6), (600, 80, 0.4)),
    seed=3,
)
# A matched pair: identical seed, so identical fragments, differing *only* in
# which mate carries the R1 flag. Any difference in output is then attributable
# to orientation and nothing else -- which is the whole point of the pair, and
# is lost the moment the seeds diverge.
_ORIENTATION_SEED = 4
ALL_FORWARD = SynthProfile(
    "all_forward", reverse_r1_fraction=0.0, seed=_ORIENTATION_SEED
)
ALL_REVERSE = SynthProfile(
    "all_reverse", reverse_r1_fraction=1.0, seed=_ORIENTATION_SEED
)
STEREOTYPED_ENDS = SynthProfile("stereotyped_ends", stereotyped_ends=True, seed=6)


def make_reference(path: Path, contigs: Sequence[Contig], seed: int = 7) -> Path:
    """Reproducible reference. A CpG-rich block at 100k-120k on the first
    contig gives UXM something to classify."""
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for contig in contigs:
            fh.write(f">{contig.name}\n")
            bases = []
            for i in range(contig.length):
                if contig.name == contigs[0].name and 100_000 <= i < 120_000:
                    bases.append("CG"[i % 2])
                else:
                    bases.append(rng.choice("ACGT"))
            seq = "".join(bases)
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return path


def make_fragments(
    profile: SynthProfile, contigs: Sequence[Contig]
) -> List[FragmentSpec]:
    rng = random.Random(profile.seed)
    weights = [w for _, _, w in profile.length_mixture]
    frags: List[FragmentSpec] = []

    for i in range(profile.n_fragments):
        mean, sd, _ = rng.choices(profile.length_mixture, weights=weights)[0]
        length = max(65, min(1000, int(round(rng.gauss(mean, sd)))))
        contig = contigs[0] if profile.stereotyped_ends else contigs[i % len(contigs)]
        if profile.stereotyped_ends:
            # make_reference lays a CG dinucleotide repeat over 100k-120k, so
            # every even offset in it starts with the same CGCG 4-mer. Odd
            # offsets would give GCGC -- also uniform, but mixing the two would
            # halve the stereotypy for no reason.
            start = CPG_BLOCK_START + 2 * ((i * 7) % ((CPG_BLOCK_LEN - length) // 2))
        else:
            # Spread fragments out so interval-overlap tests are unambiguous,
            # and keep clear of the contig ends.
            span = contig.length - length - 2000
            start = 1000 + (i * 997) % max(span, 1)
        frags.append(
            FragmentSpec(
                chrom=contig.name,
                start=start,
                length=length,
                r1_reverse=rng.random() < profile.reverse_r1_fraction,
            )
        )
    return frags


def _header(contigs: Sequence[Contig]) -> Dict:
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c.name, "LN": c.length} for c in contigs],
    }


def write_bam(
    path: Path,
    fragments: Sequence[FragmentSpec],
    contigs: Sequence[Contig],
    reference: Path,
    read_len: int = DEFAULT_READ_LEN,
) -> Path:
    """Write a coordinate-sorted, indexed PE BAM.

    Both mates are emitted for every fragment. Which mate is R1 is decided by
    ``FragmentSpec.r1_reverse``; the fragment occupies the *same* interval
    either way, which is what makes the placement test meaningful.
    """
    contig_ids = {c.name: i for i, c in enumerate(contigs)}
    unsorted = path.with_suffix(".unsorted.bam")
    fasta = pysam.FastaFile(str(reference))

    with pysam.AlignmentFile(str(unsorted), "wb", header=_header(contigs)) as bam:
        for i, frag in enumerate(fragments):
            left_start = frag.start
            right_start = frag.end - read_len
            left_seq = fasta.fetch(
                frag.chrom, left_start, left_start + read_len
            ).upper()
            right_seq = fasta.fetch(
                frag.chrom, right_start, right_start + read_len
            ).upper()

            # The leftmost mate is R1 unless the spec says otherwise.
            if frag.r1_reverse:
                mates = [
                    (right_start, right_seq, R1_REV, True, -frag.length),
                    (left_start, left_seq, R2_FWD, False, frag.length),
                ]
            else:
                mates = [
                    (left_start, left_seq, R1_FWD, False, frag.length),
                    (right_start, right_seq, R2_REV, True, -frag.length),
                ]

            for pos, seq, flag, _is_rev, tlen in mates:
                mate_pos = right_start if pos == left_start else left_start
                rec = pysam.AlignedSegment()
                rec.query_name = f"frag{i}"
                rec.query_sequence = seq
                rec.flag = flag
                rec.reference_id = contig_ids[frag.chrom]
                rec.reference_start = pos
                rec.mapping_quality = 60
                rec.cigartuples = [(0, read_len)]
                rec.next_reference_id = contig_ids[frag.chrom]
                rec.next_reference_start = mate_pos
                rec.template_length = tlen
                rec.query_qualities = pysam.qualitystring_to_array("I" * read_len)
                if frag.xm is not None:
                    rec.set_tag("XM", frag.xm)
                bam.write(rec)

    fasta.close()
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    unsorted.unlink(missing_ok=True)
    return path


def write_truth_bed(path: Path, fragments: Sequence[FragmentSpec]) -> Path:
    """The intervals extraction *should* emit: leftmost to rightmost base of
    the pair, 0-based half-open, independent of which mate is R1."""
    rows = sorted((f.chrom, f.start, f.end) for f in fragments)
    with open(path, "w") as fh:
        for chrom, start, end in rows:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    return path


@dataclass
class AssetBundle:
    """Minimal BED assets on the synthetic contigs.

    Built rather than borrowed so the suite never reaches for
    ``src/krewlyzer/data`` -- those are Git LFS pointers in a fresh clone, and
    a test that silently needs them is a test that silently skips.
    """

    bins: Path
    arms: Path
    targets: Path
    genes: Path
    alu: Path
    anchors: Path
    regions: Path


def make_assets(directory: Path, contigs: Sequence[Contig]) -> AssetBundle:
    directory.mkdir(parents=True, exist_ok=True)
    bin_size = 100_000

    bins = directory / "bins.bed"
    with open(bins, "w") as fh:
        for c in contigs:
            for start in range(0, c.length, bin_size):
                fh.write(f"{c.name}\t{start}\t{min(start + bin_size, c.length)}\n")

    arms = directory / "arms.bed"
    with open(arms, "w") as fh:
        for c in contigs:
            mid = c.length // 2
            fh.write(f"{c.name}\t0\t{mid}\t{c.name}p\n")
            fh.write(f"{c.name}\t{mid}\t{c.length}\t{c.name}q\n")

    targets = directory / "targets.bed"
    genes = directory / "genes.bed"
    with open(targets, "w") as ft, open(genes, "w") as fg:
        for gi, c in enumerate(contigs):
            for ei in range(4):
                start = 5_000 + ei * 20_000
                end = start + 2_000
                ft.write(f"{c.name}\t{start}\t{end}\n")
                # 8-column WGS-style gene BED so strand is present.
                strand = "+" if gi % 2 == 0 else "-"
                fg.write(
                    f"{c.name}\t{start}\t{end}\tENST{gi:05d}\tNM_{gi:05d}\t"
                    f"GENE{gi}\t{ei}\t{strand}\n"
                )

    alu = directory / "alu.bed"
    with open(alu, "w") as fh:
        for c in contigs:
            for start in range(20_000, c.length - 5_000, 10_000):
                fh.write(f"{c.name}\t{start}\t{start + 300}\tAluY\t0\t+\n")

    anchors = directory / "anchors.bed"
    with open(anchors, "w") as fh:
        for c in contigs:
            for i, start in enumerate(range(30_000, c.length - 5_000, 25_000)):
                kind = "TSS" if i % 2 == 0 else "CTCF"
                fh.write(f"{c.name}\t{start}\t{start + 1}\t{kind}\t0\t+\n")

    regions = directory / "regions.bed"
    with open(regions, "w") as fh:
        for c in contigs:
            for start in range(40_000, c.length - 5_000, 30_000):
                fh.write(f"{c.name}\t{start}\t{start + 2_000}\tTISSUE_A\n")

    return AssetBundle(
        bins=bins,
        arms=arms,
        targets=targets,
        genes=genes,
        alu=alu,
        anchors=anchors,
        regions=regions,
    )


@dataclass
class SynthSample:
    profile: SynthProfile
    directory: Path
    reference: Path
    bam: Path
    truth_bed: Path
    fragments: List[FragmentSpec]
    contigs: List[Contig]
    assets: AssetBundle


DEFAULT_CONTIGS = (Contig("chr1", 300_000), Contig("chr2", 300_000))


def make_sample(
    directory: Path,
    profile: SynthProfile = BASELINE,
    contigs: Sequence[Contig] = DEFAULT_CONTIGS,
) -> SynthSample:
    directory.mkdir(parents=True, exist_ok=True)
    reference = make_reference(directory / "ref.fa", contigs)
    fragments = make_fragments(profile, contigs)
    bam = write_bam(
        directory / f"{profile.name}.bam",
        fragments,
        contigs,
        reference,
        read_len=profile.read_len,
    )
    truth = write_truth_bed(directory / f"{profile.name}.truth.bed", fragments)
    assets = make_assets(directory / "assets", contigs)
    return SynthSample(
        profile=profile,
        directory=directory,
        reference=reference,
        bam=bam,
        truth_bed=truth,
        fragments=list(fragments),
        contigs=list(contigs),
        assets=assets,
    )
