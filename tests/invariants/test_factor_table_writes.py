"""A GC factor table exists exactly when a run says it wrote one, and only extraction writes one.

The tools reading a sample's folder take the factor table beside the BED unless extraction's record
says it wrote none (`resolve_gc_factors`), so a table left behind can be used, and a log line naming
a table that was not written sends whoever reads the log after one that does not exist. The ways
that went wrong:

- With `-F both`, Rust wrote the TSV table and then failed on the Parquet one. The record said
  `gc_correction_computed=False` while a `.tsv` table sat in the folder for every tool to use.
- run-all's on-target factor write logged "Wrote 0 ... factors", and the PON build's two writes
  "Computed 0", when nothing was fitted and nothing written.
- A population with no observations at all (a panel whose targets no fragment reaches) was
  skipped at every call site without a word.
- A feature step that found no off-target table beside its BED refitted one from the whole BED,
  on-target fragments included, and corrected with it: after run-all's or the PON build's own
  extraction had fitted none, in every single tool, and in the next run-all over a folder a single
  tool had refitted into, which took that table for extraction's.
"""

from __future__ import annotations

import logging
import re

import pandas as pd
import pysam
import pytest

from conftest import requires_data

from .pon_donor import build_donor, dense_panel, extract, factor_tables
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
    assert not resolve_gc_factors(bed, "off_target").found


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
        build_donor(sample, out, sample.assets.targets)
    assert not any("correction_factors" in p.name for p in out.iterdir())
    assert _claims(caplog) == []
    assert _warned_none(caplog, "off-target")
    assert _warned_none(caplog, "on-target")


def test_extract_names_no_on_target_table_it_did_not_write(sample, tmp_path, caplog):
    out = tmp_path / "out"
    with caplog.at_level(logging.DEBUG):
        result = extract(sample, out, sample.assets.targets)
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
        result = extract(sample, out, _unreached_targets(tmp_path))
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
        build_donor(sample, tmp_path / "out", targets)
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


# No feature step fits a factor table: with none beside the BED, every path -- run-all, the PON
# build, a single tool -- runs without off-target correction, and only extraction's table is used.


def test_run_all_writes_no_table_its_extraction_could_not_fit(tmp_path, caplog):
    """Extraction fits no off-target table (the on-target one it fits is right), so every feature
    runs without off-target correction, a warning says how to get one, and the record says so; a
    rerun into the folder finds no table beside the BED and does the same.
    """
    dense, targets = dense_panel(tmp_path)
    out = tmp_path / "out"
    for run in ("first", "rerun"):
        caplog.clear()
        with caplog.at_level(
            logging.WARNING, logger="krewlyzer.core.unified_processor"
        ):
            code = synth_run_all(
                dense, out, skip_target_regions=False, target_regions=targets
            )
        assert code is None, f"{run} run-all exited {code}"
        warned = [
            r.getMessage() for r in caplog.records if r.levelno == logging.WARNING
        ]
        assert any(
            "GC factors (off-target)" in m and "re-extract with GC correction" in m
            # the commonest cause, which re-extracting does not cure
            and "nothing to fit" in m
            for m in warned
        ), f"{run}: no warning with the remedy for each cause"
        record = pd.read_parquet(out / f"{SAMPLE}.metadata.parquet").iloc[0]
        assert not record["gc_correction_computed"], f"{run}: extraction fitted one"
        assert factor_tables(out, "off_target") == [], run
        assert str(record["gc_factors_offtarget"]).startswith("none:"), run
        assert not record["gc_correction_applied_offtarget"], run


def test_run_all_corrects_with_the_table_its_extraction_fitted(tmp_path):
    """The contrast: targets on one contig leave the other off-target, extraction fits its table,
    and run-all records correcting with it. Without this, the test above passes on a run-all that
    never finds a table; that the table it finds changes the features is the next test's.
    """
    dense, _ = dense_panel(tmp_path)
    targets = tmp_path / "chr1.bed"
    targets.write_text("chr1\t0\t300000\n")
    out = tmp_path / "out"
    code = synth_run_all(dense, out, skip_target_regions=False, target_regions=targets)
    assert code is None
    record = pd.read_parquet(out / f"{SAMPLE}.metadata.parquet").iloc[0]
    assert record["gc_correction_computed"]
    assert factor_tables(out, "off_target") != []
    assert record["gc_factors_offtarget"] in ("parquet", "tsv.gz", "tsv")
    assert record["gc_correction_applied_offtarget"]


def test_process_sample_writes_no_table_its_extraction_could_not_fit(tmp_path):
    """The PON build's own extraction, as run-all's: no off-target table is fitted from the BED."""
    dense, targets = dense_panel(tmp_path)
    out = tmp_path / "out"
    build_donor(dense, out, targets)
    assert factor_tables(out, "off_target") == []


def _arms(sample):
    return sample.directory / "arms.validated.bed"  # written by synth_run_all


def _single_tool(bed, out, targets, arms, gc_correct=True):
    """`krewlyzer fsd`'s path: run_features as a single tool calls it."""
    from krewlyzer.core.unified_processor import run_features

    return run_features(
        bed,
        out,
        SAMPLE,
        enable_fsd=True,
        fsd_arms=arms,
        target_regions=targets,
        gc_correct=gc_correct,
        output_format="parquet",
        compress=False,
    )


def test_the_table_extraction_fitted_changes_the_features(tmp_path):
    """Invariant 1 for the table's way into the Rust pass: a record saying a table was applied
    proves nothing unless the features change with it. On the panel whose extraction fits an
    off-target table, FSD from that BED differs with and without GC correction, and run-all's
    is the corrected one.
    """
    from krewlyzer.core.output_utils import resolve_table_path

    dense, _ = dense_panel(tmp_path)
    chr1 = tmp_path / "chr1.bed"
    chr1.write_text("chr1\t0\t300000\n")
    out = tmp_path / "out"
    code = synth_run_all(dense, out, skip_target_regions=False, target_regions=chr1)
    assert code is None
    fsd = {}
    for gc_correct in (True, False):
        tool = tmp_path / f"gc_{gc_correct}"
        tool.mkdir()
        outputs = _single_tool(
            out / f"{SAMPLE}.bed.gz", tool, chr1, _arms(dense), gc_correct=gc_correct
        )
        assert outputs.gc_factors_choice["off_target"].found is gc_correct
        fsd[gc_correct] = pd.read_parquet(resolve_table_path(outputs.fsd))
    assert not fsd[True].equals(fsd[False]), "the table changed nothing"
    pd.testing.assert_frame_equal(
        fsd[True], pd.read_parquet(resolve_table_path(out / outputs.fsd.name))
    )


def test_a_single_tool_corrects_as_run_all_does(tmp_path):
    """Invariant 6 on the panel whose extraction fits no off-target table: a single tool reading
    that extraction's BED fits none either, so its table equals run-all's."""
    dense, targets = dense_panel(tmp_path)
    out = tmp_path / "out"
    code = synth_run_all(dense, out, skip_target_regions=False, target_regions=targets)
    assert code is None
    tool = tmp_path / "tool"
    tool.mkdir()
    from krewlyzer.core.output_utils import resolve_table_path

    outputs = _single_tool(out / f"{SAMPLE}.bed.gz", tool, targets, _arms(dense))
    assert factor_tables(tool, "off_target") == []
    assert not outputs.gc_factors_choice["off_target"].found
    # `outputs.fsd` names the .tsv; each run wrote the Parquet table it resolves to.
    pd.testing.assert_frame_equal(
        pd.read_parquet(resolve_table_path(outputs.fsd)),
        pd.read_parquet(resolve_table_path(out / outputs.fsd.name)),
    )


def test_a_single_tool_leaves_no_table_a_run_all_rerun_takes_for_extractions(tmp_path):
    """A single tool run into run-all's folder writes no table beside the BED, so the next run-all
    over that folder still finds none and records it."""
    dense, targets = dense_panel(tmp_path)
    out = tmp_path / "out"
    code = synth_run_all(dense, out, skip_target_regions=False, target_regions=targets)
    assert code is None
    _single_tool(out / f"{SAMPLE}.bed.gz", out, targets, _arms(dense))
    assert factor_tables(out, "off_target") == []
    code = synth_run_all(dense, out, skip_target_regions=False, target_regions=targets)
    assert code is None
    record = pd.read_parquet(out / f"{SAMPLE}.metadata.parquet").iloc[0]
    assert str(record["gc_factors_offtarget"]).startswith("none:")
    assert not record["gc_correction_applied_offtarget"]


def test_a_run_all_rerun_does_not_take_an_earlier_builds_refit_for_extractions(
    tmp_path, caplog
):
    """Leftover state from builds before the refit's removal: after extraction fitted no table,
    run-all's feature step refitted one into its own folder, which is the BED's, under
    extraction's name. A rerun reuses the BED its record vouches for, so it found that table
    and corrected with it, beside a record saying extraction computed none. Extraction's record
    is what says whether the table is its own: the rerun runs uncorrected, warns, and records it.
    """
    from krewlyzer.core.output_utils import resolve_table_path

    dense, targets = dense_panel(tmp_path)
    chr1 = tmp_path / "chr1.bed"
    chr1.write_text("chr1\t0\t300000\n")
    fitted = (
        tmp_path / "fitted"
    )  # a real table, standing in for the earlier build's refit
    assert (
        synth_run_all(dense, fitted, skip_target_regions=False, target_regions=chr1)
        is None
    )
    out = tmp_path / "out"
    assert (
        synth_run_all(dense, out, skip_target_regions=False, target_regions=targets)
        is None
    )
    fsd = resolve_table_path(out / f"{SAMPLE}.FSD.tsv")
    uncorrected = pd.read_parquet(fsd)
    for name in factor_tables(fitted, "off_target"):
        (out / name).write_bytes((fitted / name).read_bytes())
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="krewlyzer.core.unified_processor"):
        code = synth_run_all(
            dense, out, skip_target_regions=False, target_regions=targets
        )
    assert code is None
    record = pd.read_parquet(out / f"{SAMPLE}.metadata.parquet").iloc[0]
    assert not record["gc_correction_computed"]
    assert str(record["gc_factors_offtarget"]).startswith("none:")
    assert not record["gc_correction_applied_offtarget"]
    assert any(
        "not extraction's" in r.getMessage() and "re-extract" in r.getMessage()
        for r in caplog.records
        if r.levelno == logging.WARNING
    ), "the disowned table went unsaid"
    pd.testing.assert_frame_equal(pd.read_parquet(fsd), uncorrected)


def test_run_all_says_how_to_get_the_on_target_table_it_lacks(sample, tmp_path, caplog):
    """The on-target warning names the remedy, as the off-target one does."""
    out = tmp_path / "out"
    with caplog.at_level(logging.WARNING, logger="krewlyzer.core.unified_processor"):
        code = synth_run_all(
            sample,
            out,
            skip_target_regions=False,
            target_regions=sample.assets.targets,
        )
    assert code is None
    warned = [
        r.getMessage()
        for r in caplog.records
        if r.levelno == logging.WARNING and "GC factors (on-target)" in r.getMessage()
    ]
    assert warned, "a missing on-target table went unsaid"
    assert all(
        "re-extract with GC correction" in m and "nothing to fit" in m for m in warned
    ), warned


def test_a_table_in_the_output_folder_but_not_beside_the_bed_is_not_used(tmp_path):
    """A BED's table is the one beside it. Another extraction's table in the folder a tool writes
    to describes other fragments, so the BED runs uncorrected rather than corrected with it.
    """
    dense, targets = dense_panel(tmp_path)
    chr1 = tmp_path / "chr1.bed"
    chr1.write_text("chr1\t0\t300000\n")
    fitted = tmp_path / "fitted"
    code = synth_run_all(dense, fitted, skip_target_regions=False, target_regions=chr1)
    assert code is None
    assert factor_tables(fitted, "off_target") != []
    extracted = tmp_path / "extracted"
    result = extract(dense, extracted, targets)
    assert result.exit_code == 0, result.output
    tool = tmp_path / "tool"
    tool.mkdir()
    for name in factor_tables(fitted, "off_target"):
        (tool / name).write_bytes((fitted / name).read_bytes())
    outputs = _single_tool(extracted / f"{SAMPLE}.bed.gz", tool, targets, _arms(dense))
    assert not outputs.gc_factors_choice["off_target"].found


def test_process_sample_fits_no_table_from_a_pre_extracted_bed(tmp_path):
    """The PON build given extraction's BED instead of the BAM: as from the BAM, no table."""
    dense, targets = dense_panel(tmp_path)
    extracted = tmp_path / "extracted"
    result = extract(dense, extracted, targets)
    assert result.exit_code == 0, result.output
    out = tmp_path / "out"
    build_donor(dense, out, targets, input_path=extracted / f"{SAMPLE}.bed.gz")
    assert factor_tables(out, "off_target") == []


@pytest.mark.parametrize("given", ["bam", "bed"])
def test_the_pon_build_refuses_a_donor_with_no_off_target_table(tmp_path, given):
    """A donor whose extraction fitted no off-target table is refused by the PON's GC gate,
    whichever input the build was given (`build-pon` then stops, or skips it with
    `--allow-failures`)."""
    dense, targets = dense_panel(tmp_path)
    # synth_run_all writes the arms file FSD reads; its own outputs are not used here.
    code = synth_run_all(
        dense, tmp_path / "arms_run", skip_target_regions=False, target_regions=targets
    )
    assert code is None
    input_path = None
    if given == "bed":
        extracted = tmp_path / "extracted"
        result = extract(dense, extracted, targets)
        assert result.exit_code == 0, result.output
        input_path = extracted / f"{SAMPLE}.bed.gz"
    with pytest.raises(RuntimeError, match="no off-target GC factor table"):
        build_donor(
            dense, tmp_path / "out", targets, input_path=input_path, arms=_arms(dense)
        )
