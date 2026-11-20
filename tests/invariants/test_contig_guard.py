"""One rule matches BAM contigs to the reference FASTA (rust/src/contigs.rs).

Extract matched contigs by exact name. A BAM naming `chr1` read with a FASTA naming `1` gave
every fragment GC 0.0 -- no length bin could be fitted, nothing was corrected, nothing said so.
Names now resolve under either `chr` spelling and lengths must agree, so the remaining failures
are a genuinely wrong reference (refused before any read is processed) and contigs the reference
lacks (decoys, alts: counted, and left out of the GC observations).

A target BED naming `1` against a `chr1` BAM put every fragment off-target, silently, while the
feature code (which strips `chr`) still split on/off-target; the bundled `chr`-named blacklist
was never applied to `1`-named BAMs. Both now match either spelling.
"""

import gzip

import pysam
import pytest

from krewlyzer import _core
from krewlyzer.core.gc_assets import CONTIG_MISMATCH_WARN_SHARE, check_reference_contigs

from .synth import DEFAULT_CONTIGS, Contig, make_reference, make_sample

pytestmark = [pytest.mark.invariant, pytest.mark.rust]

#: A contig a primary-assembly reference lacks, as decoys and alts are.
DECOY = Contig("chrUn_gl000220", 200_000)


def _extract(sample, reference, out=None, targets=None, exclude=None):
    return _core.extract_motif.process_bam_parallel(
        str(sample.bam),
        str(reference),
        20,
        65,
        1000,
        4,
        1,
        str(out) if out else None,
        "enable",  # end and breakpoint motifs: the breakpoint context is a FASTA fetch too
        str(exclude) if exclude else None,
        str(targets) if targets else None,
        True,
        True,
        True,
    )


def _bed_lines(path):
    with gzip.open(path, "rt") as fh:  # extract writes BGZF
        return fh.read().splitlines()


def _other_spelling(name):
    return name.removeprefix("chr") if name.startswith("chr") else f"chr{name}"


@pytest.mark.parametrize("bam_says_chr", [True, False], ids=["chr-bam", "nochr-bam"])
def test_a_reference_spelling_contigs_the_other_way_gives_the_same_result(
    tmp_path, bam_says_chr
):
    contigs = DEFAULT_CONTIGS
    if not bam_says_chr:
        contigs = tuple(
            Contig(_other_spelling(c.name), c.length) for c in DEFAULT_CONTIGS
        )
    sample = make_sample(tmp_path / "s", contigs=contigs)
    # Same seed and lengths, so the same sequences under the other names.
    other = make_reference(
        tmp_path / "other.fa",
        [Contig(_other_spelling(c.name), c.length) for c in contigs],
    )
    same = _extract(sample, sample.reference, tmp_path / "same.bed.gz")
    renamed = _extract(sample, other, tmp_path / "renamed.bed.gz")

    gc = [float(line.split("\t")[3]) for line in _bed_lines(tmp_path / "same.bed.gz")]
    assert (
        len(set(gc)) > 1 and same[2]
    ), "GC and breakpoint motifs must be real to compare"
    assert renamed[7] == {"fragments_on_missing_contig": 0, "gc_fetch_failed": 0}
    assert renamed[:7] == same[:7]  # counts, end and breakpoint motifs, GC observations
    assert _bed_lines(tmp_path / "renamed.bed.gz") == _bed_lines(
        tmp_path / "same.bed.gz"
    )


def test_a_contig_with_another_length_in_the_reference_is_refused(tmp_path):
    """Same names, other lengths: another assembly (hg38 against hg19), which no name rule sees."""
    sample = make_sample(tmp_path / "s")
    other = make_reference(
        tmp_path / "other.fa", [Contig(c.name, c.length - 1000) for c in sample.contigs]
    )
    with pytest.raises(ValueError, match="different assembly"):
        _extract(sample, other, tmp_path / "x.bed.gz")


def test_a_primary_contig_missing_from_the_reference_is_refused(tmp_path):
    sample = make_sample(tmp_path / "s")
    partial = make_reference(tmp_path / "partial.fa", sample.contigs[:1])
    with pytest.raises(ValueError, match="chr2 is in the BAM but not in the reference"):
        _extract(sample, partial, tmp_path / "x.bed.gz")


def test_fragments_on_a_contig_the_reference_lacks_are_counted_and_kept_out_of_gc(
    tmp_path,
):
    sample = make_sample(tmp_path / "s", contigs=(*DEFAULT_CONTIGS, DECOY))
    primary = make_reference(tmp_path / "primary.fa", DEFAULT_CONTIGS)
    result = _extract(sample, primary, tmp_path / "d.bed.gz")

    on_decoy = sum(
        1
        for line in _bed_lines(tmp_path / "d.bed.gz")
        if line.startswith(f"{DECOY.name}\t")
    )
    assert 0 < on_decoy < result[0]
    assert result[7] == {"fragments_on_missing_contig": on_decoy, "gc_fetch_failed": 0}
    observed = sum(result[3].values()) + sum(result[6].values())
    assert (
        observed == result[0] - on_decoy
    ), "an unresolved GC must not enter the model as 0%"


def test_contig_counts_fill_without_a_bed(tmp_path):
    """`krewlyzer motif` writes no BED: it reported 0 fragments, 0 of them on missing contigs."""
    sample = make_sample(tmp_path / "s", contigs=(*DEFAULT_CONTIGS, DECOY))
    primary = make_reference(tmp_path / "primary.fa", DEFAULT_CONTIGS)
    with_bed = _extract(sample, primary, tmp_path / "d.bed.gz")
    without = _extract(sample, primary)
    assert with_bed[7]["fragments_on_missing_contig"] > 0
    assert (without[0], without[7]) == (with_bed[0], with_bed[7])


def test_a_share_on_contigs_the_reference_lacks_is_warned_about_not_refused(caplog):
    assert check_reference_contigs(0, 0, 1000) == 0.0
    assert "reference" not in caplog.text
    share = check_reference_contigs(
        20, 0, 1000
    )  # 2%: decoys absent from a primary assembly
    assert share > CONTIG_MISMATCH_WARN_SHARE
    assert "missing from the reference" in caplog.text
    assert check_reference_contigs(990, 0, 1000) == pytest.approx(
        0.99
    )  # warned, not raised


def test_any_failed_reference_fetch_is_warned_about(caplog):
    check_reference_contigs(0, 3, 1000)
    assert "3 fragments" in caplog.text and "could not be fetched" in caplog.text


def _factor_table(path):
    """A strict-format factor table whose factor rises with GC, so GC weighting shows."""
    header = "length_bin_min\tlength_bin_max\tgc_percent\tobserved\texpected\tcorrection_factor"
    rows = [
        f"{lo}\t{lo + 4}\t{gc}\t100\t100\t{1 + gc / 100}"
        for lo in range(60, 1000, 5)
        for gc in range(101)
    ]
    path.write_text("\n".join([header, *rows]) + "\n")
    return path


def _mfsd(sample, reference, factors, out):
    """mFSD at a reference-matching SNV inside the first synthetic fragment."""
    frag = sample.fragments[0]
    pos = frag.start + 6  # 1-based, within the fragment's first read
    ref = pysam.FastaFile(str(sample.reference)).fetch(frag.chrom, pos - 1, pos)
    vcf = out.with_suffix(".vcf")
    vcf.write_text(f"{frag.chrom}\t{pos}\t.\t{ref}\t{'C' if ref != 'C' else 'G'}\n")
    _core.mfsd.calculate_mfsd(
        str(sample.bam),
        str(vcf),
        str(out),
        "vcf",
        20,
        65,
        1000,
        False,
        str(reference),
        str(factors) if factors else None,
        False,
        False,
        True,
        20,
    )
    return out.read_text()


def test_mfsd_gc_lookups_resolve_the_other_spelling(tmp_path):
    sample = make_sample(tmp_path / "s")
    other = make_reference(
        tmp_path / "other.fa",
        [Contig(_other_spelling(c.name), c.length) for c in sample.contigs],
    )
    factors = _factor_table(tmp_path / "f.tsv")
    same = _mfsd(sample, sample.reference, factors, tmp_path / "same.tsv")
    assert same != _mfsd(
        sample, sample.reference, None, tmp_path / "plain.tsv"
    ), "GC must weigh"
    assert _mfsd(sample, other, factors, tmp_path / "renamed.tsv") == same


def test_mfsd_refuses_a_reference_of_another_assembly(tmp_path):
    sample = make_sample(tmp_path / "s")
    other = make_reference(
        tmp_path / "other.fa", [Contig(c.name, c.length - 1000) for c in sample.contigs]
    )
    with pytest.raises(ValueError, match="different assembly"):
        _mfsd(sample, other, _factor_table(tmp_path / "f.tsv"), tmp_path / "m.tsv")


def test_the_expected_gc_reference_resolves_the_reference_spelling(tmp_path):
    """Valid regions name contigs without `chr`; an hg19-style FASTA names them `chr1`. The old
    fallback fetched `1` first, which rust-htslib answers with a null-pointer `Vec`, not `Err`.
    """
    reference = make_reference(tmp_path / "ref.fa", DEFAULT_CONTIGS)
    regions = tmp_path / "valid.bed"
    regions.write_text(
        "".join(f"{_other_spelling(c.name)}\t1000\t60000\n" for c in DEFAULT_CONTIGS)
    )
    total = _core.gc.generate_ref_genome_gc(
        str(reference), str(regions), str(tmp_path / "gc.parquet")
    )
    assert total > 0


def test_the_expected_gc_reference_refuses_a_region_the_reference_lacks(tmp_path):
    reference = make_reference(tmp_path / "ref.fa", DEFAULT_CONTIGS[:1])
    regions = tmp_path / "valid.bed"
    regions.write_text("1\t1000\t60000\n2\t1000\t60000\n")
    with pytest.raises(ValueError, match="contig 2 is not in the reference"):
        _core.gc.generate_ref_genome_gc(
            str(reference), str(regions), str(tmp_path / "gc.parquet")
        )


def test_target_names_match_with_or_without_chr(tmp_path):
    sample = make_sample(tmp_path / "s")
    targets = tmp_path / "targets.bed"
    first = sample.contigs[0]
    targets.write_text(f"{first.name.removeprefix('chr')}\t0\t{first.length}\tT1\n")
    result = _extract(sample, sample.reference, tmp_path / "t.bed", targets=targets)
    on_target_observations = sum(result[6].values())
    assert on_target_observations > 0, "a `1`-named target BED must match a `chr1` BAM"


def test_blacklist_names_match_with_or_without_chr(tmp_path):
    """The bundled hg19 blacklist names contigs `chr1`; MSK-ACCESS BAMs name them `1`. Extract
    compared names exactly, so in production the blacklist matched nothing and was never applied
    (on four evidence samples it removes 0.7-1.8% of fragments once it does)."""
    sample = make_sample(tmp_path / "s")
    first = sample.contigs[0]
    exclude = tmp_path / "blacklist.bed"
    exclude.write_text(f"{first.name.removeprefix('chr')}\t0\t{first.length}\n")
    kept = _extract(sample, sample.reference, tmp_path / "x.bed", exclude=exclude)
    everything = _extract(sample, sample.reference, tmp_path / "all.bed")
    on_first = sum(
        1
        for line in _bed_lines(tmp_path / "all.bed")
        if line.split("\t", 1)[0] == first.name
    )
    assert on_first > 0
    assert (
        kept[0] == everything[0] - on_first
    ), "a `1`-named blacklist must apply to a `chr1` BAM"
