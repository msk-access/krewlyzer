"""A sample's GC choice is extraction's, and every tool that reads its BED holds to it.

The record holds one GC choice per population (`gc_factors_*`, `gc_correction_applied_*`). The
first GC-correcting tool to reach it used to write it, and every later single tool kept it,
whatever that tool did:

- `fsd`, then `wps --no-gc-correct`: an uncorrected WPS under "applied", and `build-pon
  --from-outputs` accepted the donor;
- `fsd --no-gc-correct`, then `fsc`: a corrected FSC under "none";
- `motif` first wrote no GC fields and every later tool kept that, so the donor was refused as
  "written before the GC-lookup fix";
- a current tool filled GC fields into a 0.9.x extraction's record, vouching for tables 0.9.x made.

Extraction records the choice, by the resolver every tool uses. A single tool never writes it, and
one whose choice differs stops before any work.
"""

from __future__ import annotations

import pandas as pd
import pytest
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app

from .pon_donor import dense_panel, extract, factor_tables
from .runall import SAMPLE, write_arms

pytestmark = [pytest.mark.invariant, requires_data]

GC_FIELDS = (
    "gc_factors_offtarget",
    "gc_correction_applied_offtarget",
    "gc_factors_ontarget",
    "gc_correction_applied_ontarget",
)


@pytest.fixture(scope="module")
def dense(tmp_path_factory):
    """(sample, targets): targets holding all but a sliver of the fragments, so extraction fits an
    on-target table and no off-target one."""
    return dense_panel(tmp_path_factory.mktemp("gc_record_in"))


@pytest.fixture(scope="module")
def chr1(dense):
    """Targets on chr1 only: chr2's fragments are off-target, so extraction fits both tables."""
    sample, _ = dense
    targets = sample.directory.parent / "chr1.bed"
    targets.write_text("chr1\t0\t300000\n")
    return targets


def _extracted(sample, out, targets, *extra):
    result = extract(sample, out, targets, *extra)
    assert result.exit_code == 0, result.output
    return out / f"{SAMPLE}.bed.gz"


def _record_path(out):
    return out / f"{SAMPLE}.metadata.parquet"


def _gc(out):
    row = pd.read_parquet(_record_path(out)).iloc[0]
    return {
        f: (row[f].item() if hasattr(row[f], "item") else row[f])
        for f in GC_FIELDS
        if f in row.index
    }


def _tool(sample, name, bed, out, targets, *extra):
    """A single tool over the extraction's BED, into ``out``."""
    specific = {
        "fsd": ["-a", write_arms(sample.directory / "arms.validated.bed", sample)],
        "fsc": ["-b", sample.assets.bins],
        "wps": ["--wps-anchors", sample.assets.anchors, "-B", sample.assets.alu],
        "region-entropy": ["--tfbs-regions", sample.assets.regions, "--no-atac"],
    }[name]
    args = [name, "-i", bed, "-o", out, "-s", SAMPLE, *specific, "-T", targets]
    args += ["--skip-pon", "--output-format", "parquet", *extra]
    return CliRunner().invoke(app, [str(a) for a in args])


def _written(out, name):
    stem = {"fsd": "FSD", "fsc": "FSC", "wps": "WPS", "region-entropy": "TFBS"}[name]
    return sorted(p.name for p in out.glob(f"{SAMPLE}.{stem}*"))


def _remove(out, population):
    names = factor_tables(out, population)
    assert names, f"extraction fitted no {population} table to remove"
    for name in names:
        (out / name).unlink()


def test_extraction_records_the_gc_choice_every_tool_then_keeps(dense, chr1, tmp_path):
    """Three extractions, three records (invariant 1): both tables fitted; the dense panel, whose
    off-target fragments are too few to fit; and `--no-gc-correct`. A later `fsd` run as the
    extraction was keeps the record as extraction wrote it."""
    sample, dense_targets = dense
    cases = {
        "fitted": (chr1, ()),
        "unfittable": (dense_targets, ()),
        "disabled": (chr1, ("--no-gc-correct",)),
    }
    records = {}
    for case, (targets, flags) in cases.items():
        out = tmp_path / case
        bed = _extracted(sample, out, targets, *flags)
        records[case] = _gc(out)
        assert set(records[case]) == set(GC_FIELDS), (case, records[case])
        result = _tool(sample, "fsd", bed, out, targets, *flags)
        assert result.exit_code == 0, (case, result.output)
        assert _gc(out) == records[case], case

    fitted, unfittable, disabled = (records[c] for c in cases)
    assert fitted["gc_factors_offtarget"] == fitted["gc_factors_ontarget"] == "parquet"
    assert fitted["gc_correction_applied_offtarget"] is True
    assert unfittable["gc_factors_offtarget"].startswith("none:")
    assert unfittable["gc_correction_applied_offtarget"] is False
    assert unfittable["gc_factors_ontarget"] == "parquet"
    for tag in ("offtarget", "ontarget"):
        assert disabled[f"gc_factors_{tag}"] == "none:gc correction disabled"
        assert disabled[f"gc_correction_applied_{tag}"] is False
    assert len({tuple(sorted(r.items())) for r in records.values()}) == 3


def test_motif_first_does_not_cost_the_sample_its_gc_record(dense, chr1, tmp_path):
    """`motif` corrects nothing and writes no GC field. Run first, it used to leave the record
    without any, every later tool kept that, and the donor was refused as written before the
    GC-lookup fix."""
    from krewlyzer.pon.from_outputs import gc_reason

    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    args = ["motif", "-i", sample.bam, "-r", sample.reference, "-o", out, "-s", SAMPLE]
    args += ["-T", chr1, "-x", sample.assets.blacklist, "--skip-pon"]
    args += ["--output-format", "parquet"]
    result = CliRunner().invoke(app, [str(a) for a in args])
    assert result.exit_code == 0, result.output
    assert _tool(sample, "fsd", bed, out, chr1).exit_code == 0
    assert _gc(out)["gc_factors_offtarget"] == "parquet"
    assert gc_reason(out, SAMPLE) is None


@pytest.mark.parametrize(
    "tool, flags, removed, population",
    [
        ("fsd", ("--no-gc-correct",), None, "off-target"),
        ("region-entropy", ("--no-gc-correct",), None, "off-target"),
        ("wps", (), "off_target", "off-target"),
        ("region-entropy", (), "on_target", "on-target"),
    ],
    ids=[
        "fsd-no-gc",
        "region-entropy-no-gc",
        "wps-table-removed",
        "split-table-removed",
    ],
)
def test_a_tool_that_would_correct_otherwise_stops_before_any_work(
    dense, chr1, tmp_path, tool, flags, removed, population
):
    """A tool whose tables differ from the ones the record names would leave features corrected two
    ways under one record. It stops with a message, before writing anything, as a tool with other
    settings does. region-entropy's split corrects its on-target half with the on-target table.
    """
    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    if removed:
        _remove(out, removed)
    before = _record_path(out).read_bytes()
    result = _tool(sample, tool, bed, out, chr1, *flags)
    assert result.exit_code == 1, result.output
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    assert population in " ".join(result.output.split())
    assert _written(out, tool) == []
    assert _record_path(out).read_bytes() == before


def test_fsd_does_not_compare_an_on_target_table_it_never_reads(dense, chr1, tmp_path):
    """The Rust pass corrects every table it writes, its on-target ones included, with the
    off-target table; FSD reads no other. A missing on-target table is no reason to stop it.
    """
    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    _remove(out, "on_target")
    result = _tool(sample, "fsd", bed, out, chr1)
    assert result.exit_code == 0, result.output
    assert _written(out, "fsd")


@pytest.mark.parametrize("tool", ["fsd", "region-entropy"])
def test_a_table_in_the_output_folder_is_never_used(dense, chr1, tmp_path, tool):
    """Only a table beside the BED corrects it: one in the folder a tool writes to may be another
    extraction's, fitted from other fragments, and the output folder was once searched by sample
    name too. `out/` holds both fitted tables named for the sample, the BED has none beside it, and
    every table written there equals the one written into an empty folder, through both call sites
    that resolve tables: run_features (`fsd`) and `region-entropy` (genome-wide TFBS). Both tools
    correct with the off-target table here; without `--assay` no published table reads the on-target
    one, and the bundled assays' targets lie off the synthetic contigs."""
    sample, _ = dense
    fitted, plain = tmp_path / "fitted", tmp_path / "plain"
    _extracted(sample, fitted, chr1)
    bed = _extracted(sample, plain, chr1, "--no-gc-correct")
    assert factor_tables(plain, "off_target", "on_target") == []
    out, empty = tmp_path / "out", tmp_path / "empty"
    out.mkdir()
    for population in ("off_target", "on_target"):
        (name,) = factor_tables(fitted, population)
        factors = pd.read_parquet(fitted / name)["correction_factor"]
        assert factors.nunique() > 1, f"a flat {population} table"
        (out / name).write_bytes((fitted / name).read_bytes())
    for folder in (out, empty):
        result = _tool(sample, tool, bed, folder, chr1)
        assert result.exit_code == 0, result.output
    written = [name for name in _written(empty, tool) if name.endswith(".parquet")]
    assert (
        written
        and [n for n in _written(out, tool) if n.endswith(".parquet")] == written
    )
    for name in written:
        pd.testing.assert_frame_equal(
            pd.read_parquet(out / name), pd.read_parquet(empty / name), obj=name
        )


def test_two_uncorrected_runs_agree(dense, chr1, tmp_path):
    """Extracted without GC correction, the folder has no table, so a tool with correction on runs
    uncorrected too: nothing differs, nothing stops, and the record keeps extraction's reason.
    """
    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1, "--no-gc-correct")
    result = _tool(sample, "fsd", bed, out, chr1)
    assert result.exit_code == 0, result.output
    assert _gc(out)["gc_factors_offtarget"] == "none:gc correction disabled"


@pytest.mark.parametrize("stamped", [True, False], ids=["stamped", "unstamped"])
@pytest.mark.parametrize("tool", ["fsd", "region-entropy"])
def test_a_single_tool_does_not_vouch_for_a_record_written_before_the_gc_fix(
    dense, chr1, tmp_path, tool, stamped
):
    """A 0.9.x record holds no GC field and three of the six filter columns: 0.9.x run-all stamped
    the PON fields without GC ones, 0.9.x `region-entropy` never stamped, and either may have left
    panel WPS, TFBS, ATAC or gene FSC uncorrected. `build-pon --from-outputs` refuses such a
    folder on the missing column. One current tool run into it must not add the column, or the
    donor is accepted with features that tool never computed."""
    from krewlyzer.core.output_provenance import (
        TABLE_FORMATS,
        stamp_metadata,
    )
    from krewlyzer.pon.from_outputs import gc_reason

    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    frame = pd.read_parquet(_record_path(out))
    newer = (
        *GC_FIELDS,
        "filters_skip_duplicates",
        "filters_require_proper_pair",
        "filters_blacklist_sha256",
        TABLE_FORMATS,
    )
    frame.drop(columns=[c for c in newer if c in frame.columns]).to_parquet(
        _record_path(out), index=False
    )
    if stamped:  # as 0.9.x run-all left it: the PON fields, no GC columns
        assert stamp_metadata(
            out, SAMPLE, pon=None, output_format="parquet", compress=False
        )
    result = _tool(sample, tool, bed, out, chr1)
    assert result.exit_code == 0, result.output
    assert _gc(out) == {}
    assert "written before the GC-lookup fix" in (gc_reason(out, SAMPLE) or "")
