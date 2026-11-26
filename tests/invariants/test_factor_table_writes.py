"""A GC factor table from extraction exists exactly when the run says it wrote one.

The tools reading a sample's folder find its factor table by existence (`resolve_gc_factors`), so a
table left behind is used, and a log line naming a table that was not written sends whoever reads
the log after one that does not exist. Three ways that went wrong:

- With `-F both`, Rust wrote the TSV table and then failed on the Parquet one. The record said
  `gc_correction_computed=False` while a `.tsv` table sat in the folder for every tool to use.
- The on-target factor writes in run-all and `process_sample` logged "Wrote 0 ... factors" (and the
  off-target one in `process_sample` "Computed 0") when nothing was fitted and nothing written.
- A population with no observations at all (a panel whose targets no fragment reaches) was
  skipped at every call site without a word.

The two strict xfails at the end pin two defects of the other writer of a factor table, the
unified pipeline's refit (`write_gc_factor_table`'s docstring says when it runs). They are fixed
where the refit is removed (#92), not here.
"""

from __future__ import annotations

import logging
import re

import pandas as pd
import pysam
import pytest

from conftest import requires_data

from .runall import SAMPLE, synth_run_all
from .synth import BASELINE, make_sample

pytestmark = [pytest.mark.invariant, requires_data]

#: Every (length bin, GC%) well populated: every bin the reference covers can be fitted.
DENSE = {(b, gc): 10_000 for b in range(188) for gc in range(101)}

#: A log line claiming a count of zero factors, written or computed.
ZERO_CLAIM = re.compile(r"\b(Wrote|Computed) 0\b")


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("factor_writes_in"), BASELINE)


def _extraction(directory, gc_observations):
    from krewlyzer.core.sample_processor import ExtractionResult

    raw = directory / f"{SAMPLE}.bed"
    raw.write_text("chr1\t1000\t1167\t0.5\t1.0\n")
    bed_gz = directory / f"{SAMPLE}.bed.gz"
    pysam.tabix_compress(str(raw), str(bed_gz), force=True)
    raw.unlink()
    return ExtractionResult(
        sample_id=SAMPLE,
        input_path=directory / "in.bam",
        fragment_count=1,
        bed_path=bed_gz,
        gc_observations=gc_observations,
        mapq=20,
        minlen=65,
        maxlen=1000,
        skip_duplicates=True,
        require_proper_pair=True,
        blacklist_sha256="",  # none applied (`blacklist_digest`)
    )


def test_a_format_that_fails_leaves_no_table_in_the_other(tmp_path):
    """`-F both`: the Parquet write fails (its path is a directory, so it cannot be created) after
    the TSV one succeeded. The record says no table was computed, and no table may be found.
    """
    from krewlyzer.core.gc_assets import factor_table_base, resolve_gc_factors
    from krewlyzer.core.sample_processor import write_extraction_outputs

    base = factor_table_base(tmp_path, SAMPLE, "off_target")
    (tmp_path / f"{base.name}.parquet").mkdir()

    write_extraction_outputs(
        _extraction(tmp_path, DENSE), tmp_path, genome="hg19", output_format="both"
    )

    row = pd.read_csv(tmp_path / f"{SAMPLE}.metadata.tsv", sep="\t").iloc[0]
    assert bool(row["gc_correction_computed"]) is False
    assert not (tmp_path / f"{base.name}.tsv").exists(), "half a table was left"
    (tmp_path / f"{base.name}.parquet").rmdir()  # the test's own obstacle, not a table
    bed = tmp_path / f"{SAMPLE}.bed.gz"
    assert not resolve_gc_factors(bed, tmp_path, SAMPLE, "off_target").found


def test_both_formats_are_written_when_both_succeed(tmp_path):
    """The rollback is for a failure: a `-F both` run that succeeds keeps both tables."""
    from krewlyzer.core.gc_assets import factor_table_base
    from krewlyzer.core.sample_processor import write_extraction_outputs

    write_extraction_outputs(
        _extraction(tmp_path, DENSE), tmp_path, genome="hg19", output_format="both"
    )

    base = factor_table_base(tmp_path, SAMPLE, "off_target")
    assert (tmp_path / f"{base.name}.tsv").exists()
    assert (tmp_path / f"{base.name}.parquet").is_file()


def _claims(caplog):
    return [r.getMessage() for r in caplog.records if ZERO_CLAIM.search(r.getMessage())]


def _warned_none(caplog, population):
    return [
        r.getMessage()
        for r in caplog.records
        if r.levelno >= logging.WARNING
        and population in r.getMessage()
        and "no factor table written" in r.getMessage()
    ]


def _assert_warned_no_fragments(caplog, population):
    """The missing table was warned about, and for this reason rather than an unfittable one."""
    warned = _warned_none(caplog, population)
    assert warned, f"a {population} population with no fragments went unsaid"
    assert all(f"no {population} fragments" in m.lower() for m in warned), warned


def _unreached_targets(directory):
    """One target no fragment of the synthetic sample reaches: no on-target observations."""
    targets = directory / "unreached.bed"
    targets.write_text("chr2\t299990\t299999\n")
    return targets


def _whole_genome_targets(sample, directory):
    """Targets covering every contig: no off-target observations."""
    targets = directory / "everything.bed"
    targets.write_text("".join(f"{c.name}\t0\t{c.length}\n" for c in sample.contigs))
    return targets


def _extract(sample, out, targets):
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    return CliRunner().invoke(
        app,
        [
            "extract",
            "-i",
            str(sample.bam),
            "-r",
            str(sample.reference),
            "-o",
            str(out),
            "-s",
            SAMPLE,
            "-x",
            str(sample.assets.blacklist),
            "-T",
            str(targets),
            "-F",
            "parquet",
        ],
    )


def _process_sample(sample, out, targets):
    """`build-pon`'s path with every feature off: extraction and its factor tables."""
    from krewlyzer.core.sample_processor import SampleParams, process_sample

    out.mkdir()
    process_sample(
        input_path=sample.bam,
        output_dir=out,
        sample_name=SAMPLE,
        reference=sample.reference,
        params=SampleParams(threads=1, exclude_regions=sample.assets.blacklist),
        target_regions=targets,
        enable_fsc=False,
        enable_fsr=False,
        enable_fsd=False,
        enable_wps=False,
        enable_ocf=False,
        enable_tfbs=False,
        enable_atac=False,
        pon_mode=True,
        output_format="tsv",
        compress=False,
    )


def test_run_all_names_no_on_target_table_it_did_not_write(sample, tmp_path, caplog):
    """A panel run-all on the synthetic sample: too few fragments to fit any factor."""
    out = tmp_path / "out"
    with caplog.at_level(logging.DEBUG):
        code = synth_run_all(
            sample,
            out,
            skip_target_regions=False,
            target_regions=sample.assets.targets,
        )
    assert code is None
    assert not any("correction_factors.ontarget" in p.name for p in out.iterdir())
    assert _claims(caplog) == []
    assert _warned_none(caplog, "on-target"), "an unfitted on-target table went unsaid"


def test_process_sample_names_no_table_it_did_not_write(sample, tmp_path, caplog):
    """`build-pon`'s path, both populations."""
    out = tmp_path / "out"
    with caplog.at_level(logging.DEBUG):
        _process_sample(sample, out, sample.assets.targets)
    assert not any("correction_factors" in p.name for p in out.iterdir())
    assert _claims(caplog) == []
    assert _warned_none(caplog, "off-target")
    assert _warned_none(caplog, "on-target")


def test_extract_names_no_on_target_table_it_did_not_write(sample, tmp_path, caplog):
    out = tmp_path / "out"
    with caplog.at_level(logging.DEBUG):
        result = _extract(sample, out, sample.assets.targets)
    assert result.exit_code == 0, result.output
    assert not any("correction_factors" in p.name for p in out.iterdir())
    assert _claims(caplog) == []
    assert _warned_none(caplog, "on-target")


# A population with no observations at all is still expected, so its missing table is warned
# about, with that reason. Each write site is pinned: a site that skipped the write when there
# was nothing to fit said nothing, and that is the defect these tests catch.


def test_extract_warns_when_no_fragment_reaches_a_target(sample, tmp_path, caplog):
    out = tmp_path / "out"
    with caplog.at_level(logging.DEBUG):
        result = _extract(sample, out, _unreached_targets(tmp_path))
    assert result.exit_code == 0, result.output
    assert not any("correction_factors.ontarget" in p.name for p in out.iterdir())
    _assert_warned_no_fragments(caplog, "on-target")


def test_run_all_warns_when_no_fragment_reaches_a_target(sample, tmp_path, caplog):
    out = tmp_path / "out"
    with caplog.at_level(logging.DEBUG):
        code = synth_run_all(
            sample,
            out,
            skip_target_regions=False,
            target_regions=_unreached_targets(tmp_path),
        )
    assert code is None
    assert not any("correction_factors.ontarget" in p.name for p in out.iterdir())
    _assert_warned_no_fragments(caplog, "on-target")


@pytest.mark.parametrize("population", ["off-target", "on-target"])
def test_process_sample_warns_when_a_population_has_no_fragments(
    sample, tmp_path, caplog, population
):
    targets = (
        _whole_genome_targets(sample, tmp_path)
        if population == "off-target"
        else _unreached_targets(tmp_path)
    )
    with caplog.at_level(logging.DEBUG):
        _process_sample(sample, tmp_path / "out", targets)
    _assert_warned_no_fragments(caplog, population)


def test_extraction_outputs_warn_when_nothing_is_off_target(tmp_path, caplog):
    """`write_extraction_outputs`, the off-target write of extract and run-all."""
    from krewlyzer.core.sample_processor import write_extraction_outputs

    with caplog.at_level(logging.DEBUG):
        write_extraction_outputs(
            _extraction(tmp_path, {}), tmp_path, genome="hg19", output_format="parquet"
        )
    assert not any("correction_factors" in p.name for p in tmp_path.iterdir())
    _assert_warned_no_fragments(caplog, "off-target")


# Known defects of the unified pipeline's refit, pinned until fragment-model Step 2 removes it
# (#92: the Phase-1 tally is deleted there). Each test asserts the correct behaviour and is strict:
# the day the defect goes, it XPASSes, fails the build, and its marker has to come out. Only an
# assertion counts as the expected failure (`raises`); a broken premise fails outright through
# `pytest.fail`, so neither pin can outlive its defect by failing for another reason.


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="run-all refits, from the whole BED (on-target included), the table its own "
    "extraction could not fit; removed with the Phase-1 tally (#92)",
)
def test_run_all_writes_no_table_its_extraction_could_not_fit(tmp_path, caplog):
    """A panel run-all whose targets hold almost every fragment: extraction, off-target only, fits
    nothing (the on-target table it fits is right). No off-target table may then appear in the
    folder, and the record says none, on a rerun into that folder too.
    """
    from krewlyzer.core.gc_assets import FACTOR_KINDS, factor_table_base

    from .synth import Contig, SynthProfile

    profile = SynthProfile(
        "dense", n_fragments=60_000, length_mixture=((167, 2, 1.0),), seed=7
    )
    contigs = (Contig("chr1", 300_000), Contig("chr2", 300_000))
    dense = make_sample(tmp_path / "in", profile, contigs=contigs)
    targets = tmp_path / "targets.bed"
    targets.write_text("chr1\t0\t300000\nchr2\t0\t298500\n")
    out = tmp_path / "out"
    base = factor_table_base(out, SAMPLE, "off_target").name
    off_target = {f"{base}{ext}" for _, ext in FACTOR_KINDS}
    for run in ("first", "rerun"):
        with caplog.at_level(logging.DEBUG):
            code = synth_run_all(
                dense, out, skip_target_regions=False, target_regions=targets
            )
        if code is not None:
            pytest.fail(f"{run} run-all exited {code}")
        if run == "first" and not _warned_none(caplog, "off-target"):
            pytest.fail("premise broken: extraction fitted an off-target table")
        tables = sorted(p.name for p in out.iterdir() if p.name in off_target)
        assert tables == [], f"{run}: {tables}"
        record = pd.read_parquet(out / f"{SAMPLE}.metadata.parquet").iloc[0]
        assert str(record["gc_factors_offtarget"]).startswith("none"), run


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="the refit numbers chromosomes in the order it first meets them, once for the valid "
    "regions and again for the BED (pipeline.rs Phase 1); removed with the Phase-1 tally (#92)",
)
def test_the_refit_does_not_depend_on_the_beds_contig_order(tmp_path):
    """The same fragments, written chr1-first and chr2-first, refit the same table. Each fragment
    is checked against its own chromosome's valid regions, whatever order the BED meets them in
    (and whichever contigs it lacks: one missing ahead of another shifts it the same way).
    """
    import random

    from krewlyzer import _core
    from krewlyzer.assets import AssetManager

    assets = AssetManager("hg19")
    rng = random.Random(3)
    rows = [("chr1", 5_000_000, 5_000_167, 0.5)]
    for _ in range(30_000):
        start = 100_000 + rng.randrange(190_000)
        gc = round(min(max(rng.gauss(0.5, 0.05), 0.0), 1.0), 2)
        rows.append(("chr2", start, start + 167, gc))

    tables = []
    for first, second in (("chr1", "chr2"), ("chr2", "chr1")):
        raw = tmp_path / f"{first}_first.bed"
        with open(raw, "w") as fh:
            for chrom in (first, second):
                for c, s, e, g in sorted(r for r in rows if r[0] == chrom):
                    fh.write(f"{c}\t{s}\t{e}\t{g}\t1.0\n")
        bed = tmp_path / f"{first}_first.bed.gz"
        pysam.tabix_compress(str(raw), str(bed), force=True)
        table = tmp_path / f"{first}_first.correction_factors.tsv"
        _core.run_unified_pipeline(
            str(bed),
            str(assets.resolve("gc_reference")),
            str(assets.resolve("valid_regions")),
            str(table),
            None,
            output_format="tsv",
            compress=False,
            silent=True,
        )
        if first == "chr1" and not table.exists():
            pytest.fail(
                "premise broken: the refit wrote no table from the chr1-first BED"
            )
        assert table.exists(), f"no table refitted from the {first}-first BED"
        tables.append(pd.read_csv(table, sep="\t"))
    pd.testing.assert_frame_equal(tables[0], tables[1])
