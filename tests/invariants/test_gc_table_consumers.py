"""region-entropy and run-all's mFSD find, check, report and record a sample's GC factor tables
the way run-all's feature pass does (`resolve_feature_gc_tables`, `report_gc_tables`).

region-entropy and run-all's mFSD looked the tables up themselves, with `resolve_gc_factors`, and
each missed part of what the feature pass does:

- region-entropy warned of a missing off-target table without saying whether re-extracting helps
  (for a dense panel it does not). It said nothing when a panel had no on-target table and its
  on-target fragments fell back to the off-target one. And it recorded neither table: into
  extraction's record it wrote the PON fields only, and every later single tool keeps a record
  that has them, so the sample's GC columns never appeared.
- It handed Rust a table it had never read, so an unreadable one ended in a traceback.
- Its `--gc-factors` corrected with a table other than the sample's, which the sample's record,
  holding one choice, cannot describe; it is removed.
- run-all's mFSD looked the tables up a second time, unchecked, and repeated the feature pass's
  warnings without their remedy.
"""

from __future__ import annotations

import logging

import pandas as pd
import pytest
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app

from .pon_donor import dense_panel, extract, factor_tables
from .runall import SAMPLE, outcomes, synth_run_all, write_arms, write_vcf

pytestmark = [pytest.mark.invariant, requires_data]


@pytest.fixture(scope="module")
def dense(tmp_path_factory):
    """(sample, targets): a panel whose targets hold all but a sliver of the fragments, so
    extraction fits an on-target table and no off-target one."""
    return dense_panel(tmp_path_factory.mktemp("consumers_in"))


@pytest.fixture(scope="module")
def chr1(dense):
    """Targets on chr1 only: chr2's fragments are off-target, so extraction fits both tables."""
    sample, _ = dense
    targets = sample.directory.parent / "chr1.bed"
    targets.write_text("chr1\t0\t300000\n")
    return targets


def _extracted(sample, out, targets):
    result = extract(sample, out, targets)
    assert result.exit_code == 0, result.output
    return out / f"{SAMPLE}.bed.gz"


def _region_entropy(sample, bed, out, targets, *extra):
    """TFBS over the synthetic regions (the bundled ones are on hg19 coordinates), into the
    extraction's folder, whose record it stamps."""
    return CliRunner().invoke(
        app,
        [
            "region-entropy",
            "-i",
            str(bed),
            "-o",
            str(out),
            "-s",
            SAMPLE,
            "--tfbs-regions",
            str(sample.assets.regions),
            "--no-atac",
            "--skip-pon",
            "-T",
            str(targets),
            "--output-format",
            "parquet",
            *extra,
        ],
    )


def _warnings(caplog):
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]


def _record(out):
    return pd.read_parquet(out / f"{SAMPLE}.metadata.parquet").iloc[0]


def _tfbs_tables(out):
    return sorted(p.name for p in out.glob(f"{SAMPLE}.TFBS*"))


def test_without_an_off_target_table_region_entropy_warns_with_the_remedy(
    dense, tmp_path, caplog
):
    """The dense panel: extraction fitted no off-target table, and nothing could. The warning
    says so, and how to get one when something can, as run-all's does."""
    sample, targets = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, targets)
    with caplog.at_level(logging.WARNING):
        result = _region_entropy(sample, bed, out, targets)
    assert result.exit_code == 0, result.output
    assert any(
        "GC factors (off-target)" in m
        and "nothing to fit" in m
        and "re-extract with GC correction" in m
        for m in _warnings(caplog)
    ), _warnings(caplog)


def test_a_panel_without_an_on_target_table_warns_of_the_fallback(
    dense, chr1, tmp_path, caplog
):
    """No on-target table beside the BED: the split's on-target fragments are corrected with the
    off-target table, which Rust logged at debug only."""
    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    fitted = factor_tables(out, "on_target")
    assert fitted, "extraction fitted no on-target table to remove"
    for name in fitted:
        (out / name).unlink()
    with caplog.at_level(logging.WARNING):
        result = _region_entropy(sample, bed, out, chr1)
    assert result.exit_code == 0, result.output
    assert any(
        "GC factors (on-target)" in m and "off-target table" in m
        for m in _warnings(caplog)
    ), _warnings(caplog)


@pytest.mark.parametrize("panel", ["fitted", "unfittable"])
def test_region_entropy_records_the_tables_it_corrected_with(
    dense, chr1, tmp_path, caplog, panel
):
    """Into extraction's record, which no run has stamped, region-entropy records the tables it
    used, as the unified-pipeline tools do: the kind when there is one, and why not when there is
    none. Invariant 1: the two panels give two records. With every table there, nothing warns.
    """
    sample, dense_targets = dense
    targets = chr1 if panel == "fitted" else dense_targets
    out = tmp_path / "out"
    bed = _extracted(sample, out, targets)
    with caplog.at_level(logging.WARNING):
        result = _region_entropy(sample, bed, out, targets)
    assert result.exit_code == 0, result.output
    row = _record(out)
    assert row["gc_factors_ontarget"] == "parquet"  # both panels fit one
    assert bool(row["gc_correction_applied_ontarget"])
    if panel == "fitted":
        assert row["gc_factors_offtarget"] == "parquet"
        assert bool(row["gc_correction_applied_offtarget"])
        assert not [m for m in _warnings(caplog) if "GC factors" in m]
    else:
        assert str(row["gc_factors_offtarget"]).startswith("none:")
        assert not bool(row["gc_correction_applied_offtarget"])


def test_an_unreadable_table_is_refused_before_any_table_is_written(
    dense, chr1, tmp_path
):
    """A table beside the BED that cannot be read stops the run with a message, before any work,
    as run-all's feature pass does, rather than in a traceback from Rust's loader."""
    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    (name,) = factor_tables(out, "off_target")
    (out / name).write_bytes(b"not a factor table")
    result = _region_entropy(sample, bed, out, chr1)
    assert result.exit_code == 1, result.output
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    assert "cannot be read" in result.output
    assert _tfbs_tables(out) == []


def test_gc_factors_is_not_an_option(dense, chr1, tmp_path):
    """`--gc-factors` let region-entropy correct with a table other than the sample's, and the
    record, which holds one choice per sample, then said a later tool's uncorrected features were
    corrected. Every tool corrects only with extraction's table beside the BED."""
    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    (name,) = factor_tables(out, "off_target")
    result = _region_entropy(sample, bed, out, chr1, "--gc-factors", str(out / name))
    assert result.exit_code == 2, result.output
    assert _tfbs_tables(out) == []


@pytest.mark.parametrize("tool", ["fsd", "region-entropy"])
def test_a_single_tool_does_not_vouch_for_a_record_written_before_the_gc_fix(
    dense, chr1, tmp_path, tool
):
    """0.9.x run-all stamped the PON fields and no GC columns, and its panel WPS, TFBS, ATAC and
    gene FSC may be uncorrected, so `build-pon --from-outputs` refuses the folder on the missing
    column. One current tool run into it (a region-entropy rerun after a failed TFBS stage) must
    not add the column, or the donor is accepted with features that tool never computed.
    """
    from krewlyzer.core.output_provenance import stamp_metadata
    from krewlyzer.pon.from_outputs import gc_reason

    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    # The record as 0.9.x run-all left it: the PON fields, no GC columns.
    assert stamp_metadata(
        out, SAMPLE, pon=None, output_format="parquet", compress=False
    )
    assert "gc_factors_offtarget" not in _record(out)
    if tool == "fsd":
        arms = write_arms(sample.directory / "arms.validated.bed", sample)
        result = CliRunner().invoke(
            app,
            [
                "fsd",
                "-i",
                str(bed),
                "-o",
                str(out),
                "-s",
                SAMPLE,
                "-a",
                str(arms),
                "-T",
                str(chr1),
                "--skip-pon",
                "--output-format",
                "parquet",
            ],
        )
    else:
        result = _region_entropy(sample, bed, out, chr1)
    assert result.exit_code == 0, result.output
    assert "gc_factors_offtarget" not in _record(out)
    assert "written before the GC-lookup fix" in (gc_reason(out, SAMPLE) or "")


def test_region_entropy_that_ran_nothing_records_no_tables(dense, chr1, tmp_path):
    """Outcome, not intent: with neither TFBS nor ATAC run, region-entropy corrected nothing, so
    it records no GC tables."""
    sample, _ = dense
    out = tmp_path / "out"
    bed = _extracted(sample, out, chr1)
    result = _region_entropy(sample, bed, out, chr1, "--no-tfbs")
    assert result.exit_code == 0, result.output
    row = _record(out)
    assert "pon_applied" in row and "gc_factors_offtarget" not in row


def test_run_all_mfsd_does_not_repeat_the_gc_warning(dense, chr1, tmp_path, caplog):
    """With both tables mFSD takes the on-target one (variants sit in captured regions). On a
    rerun over a folder whose on-target table is gone (the BED is reused, so nothing fits it
    again), the feature pass warns once, with the remedy, and mFSD falls back to the off-target
    table and adds no warning of its own."""
    sample, _ = dense
    out = tmp_path / "out"
    options = dict(
        skip_target_regions=False,
        target_regions=chr1,
        variants=write_vcf(tmp_path / "variants.vcf"),
    )
    with caplog.at_level(logging.INFO):
        assert synth_run_all(sample, out, **options) is None
    assert any(
        r.getMessage().startswith("mFSD: GC factors (on-target)")
        for r in caplog.records
    ), "with both tables, mFSD did not take the on-target one"
    fitted = factor_tables(out, "on_target")
    assert fitted, "extraction fitted no on-target table to remove"
    for name in fitted:
        (out / name).unlink()
    caplog.clear()
    with caplog.at_level(logging.INFO):
        assert synth_run_all(sample, out, **options) is None
    about_on_target = [
        m for m in _warnings(caplog) if "on-target" in m and "GC factors" in m
    ]
    assert len(about_on_target) == 1, about_on_target
    assert any(
        r.getMessage().startswith("mFSD: GC factors (off-target)")
        for r in caplog.records
    ), "mFSD did not say which table it corrects with"
    assert outcomes(out)["mfsd"] == ("ran", "")
