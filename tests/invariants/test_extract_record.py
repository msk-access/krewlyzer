"""What `extract` writes into `{sample}.metadata.*` is what the run did, and it never replaces a
record another stage wrote into.

`{sample}.metadata.parquet` is downstream's completion marker (invariant #2). `extract` writes it
from scratch. Run into a folder run-all had finished, it dropped about thirty columns --
`krewlyzer_version`, `pon_applied`, every `{stage}_status` -- and the sample still looked
complete to a reader of the marker alone. It also recorded filters it had not used (always
20/65/1000) and `gc_correction_computed=True` for a run that fitted no factor and wrote no table.
Every one of those values was present and plausible, so each is pinned here by giving two
different runs and asserting two different records (invariant #1).
"""

from __future__ import annotations

import logging

import pandas as pd
import pysam
import pytest
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app

from .runall import SAMPLE, record, synth_run_all, write_vcf
from .synth import BASELINE, make_sample

pytestmark = [pytest.mark.invariant, requires_data]


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("extract_record_in"), BASELINE)


def _files(out):
    return {p.name: p.stat().st_mtime_ns for p in sorted(out.iterdir()) if p.is_file()}


def _cli(*args):
    return CliRunner().invoke(app, [str(a) for a in args])


def _extract(sample, out, *extra):
    """`extract` with the settings `synth_run_all` uses, so only what a test adds differs."""
    return _cli(
        "extract",
        "-i",
        sample.bam,
        "-r",
        sample.reference,
        "-o",
        out,
        "-s",
        SAMPLE,
        "-x",
        sample.assets.blacklist,
        "--skip-target-regions",
        "-F",
        "parquet",
        *extra,
    )


# ── 1. A record another stage wrote into is not replaced ─────────────────────────────────────


def test_extract_into_a_run_all_folder_is_refused(sample, tmp_path):
    """The same settings as run-all, so only the record collision can stop it."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    kept = record(out)
    assert "krewlyzer_version" in kept and "pon_applied" in kept
    before = _files(out)

    result = _extract(sample, out)

    assert result.exit_code == 1, result.output
    assert "Traceback" not in result.output
    assert _files(out) == before, "a refused extract wrote into the folder"
    assert record(out) == kept


def test_extract_into_a_folder_a_single_tool_recorded_into_is_refused(sample, tmp_path):
    """`fsr` read this BED and recorded its stage; a new BED would leave its table stale, and the
    rewrite would erase what it recorded."""
    out = tmp_path / "out"
    assert _extract(sample, out).exit_code == 0
    fsr = _cli(
        "fsr",
        "-i",
        out / f"{SAMPLE}.bed.gz",
        "-o",
        out,
        "-s",
        SAMPLE,
        "-b",
        sample.assets.bins,
        "--skip-pon",
        "--skip-target-regions",
        "-F",
        "parquet",
    )
    assert fsr.exit_code == 0, fsr.output
    assert "pon_applied" in record(out)
    before = _files(out)

    result = _extract(sample, out)

    assert result.exit_code == 1, result.output
    assert _files(out) == before


def test_extract_into_a_folder_with_only_a_stage_status_is_refused(sample, tmp_path):
    """`mfsd` records `mfsd_status` and no provenance column, so only the status half of the
    check sees it: a refusal reading `pon_applied` alone let this extract through."""
    out = tmp_path / "out"
    assert _extract(sample, out).exit_code == 0
    mfsd = _cli(
        "mfsd",
        "-i",
        sample.bam,
        "-V",
        write_vcf(tmp_path / "v.vcf"),
        "-g",
        sample.reference,
        "-o",
        out,
        "--sample-name",
        SAMPLE,
        "--output-format",
        "parquet",
    )
    assert mfsd.exit_code == 0, mfsd.output
    kept = record(out)
    assert "mfsd_status" in kept and "pon_applied" not in kept
    before = _files(out)

    result = _extract(sample, out)

    assert result.exit_code == 1, result.output
    assert _files(out) == before


def test_extract_rerun_over_its_own_record_is_allowed(sample, tmp_path):
    """The refusal is for a record something else wrote into, not for a rerun of extract."""
    out = tmp_path / "out"
    assert _extract(sample, out).exit_code == 0
    first = record(out)["timestamp"]
    result = _extract(sample, out)
    assert result.exit_code == 0, result.output
    assert record(out)["timestamp"] != first, "the rerun did not rewrite its record"


# ── 2. The record says what the run did ───────────────────────────────────────────────────────


def test_recorded_filters_are_the_ones_the_run_used(sample, tmp_path):
    runs = {
        (20, 65, 1000): tmp_path / "default",
        (30, 80, 600): tmp_path / "strict",
    }
    recorded = {}
    for (mapq, minlen, maxlen), out in runs.items():
        result = _extract(
            sample, out, "--mapq", mapq, "--minlen", minlen, "--maxlen", maxlen
        )
        assert result.exit_code == 0, result.output
        row = record(out)
        recorded[(mapq, minlen, maxlen)] = (
            int(row["filters_mapq"]),
            int(row["filters_min_length"]),
            int(row["filters_max_length"]),
        )
    for given, written in recorded.items():
        assert written == given
    assert len(set(recorded.values())) == 2, "two runs, one record"


def _extraction(directory, gc_observations):
    """An `ExtractionResult` as `extract_sample` returns one, over a one-fragment BED."""
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
    )


def test_gc_correction_computed_says_whether_a_table_was_written(tmp_path):
    """Rust returns 0 and writes nothing when no length bin has enough points to fit; the record
    said True regardless, and nothing reading the marker could tell."""
    from krewlyzer.core.gc_assets import factor_table_base
    from krewlyzer.core.sample_processor import write_extraction_outputs

    observations = {
        # Every (length bin, GC%) well populated: every bin the reference covers can be fitted.
        "dense": {(b, gc): 10_000 for b in range(188) for gc in range(101)},
        # One cell: no bin reaches the ten points a fit needs.
        "sparse": {(40, 50): 1},
    }
    computed = {}
    for name, obs in observations.items():
        out = tmp_path / name
        out.mkdir()
        write_extraction_outputs(
            _extraction(out, obs), out, genome="hg19", output_format="parquet"
        )
        row = pd.read_parquet(out / f"{SAMPLE}.metadata.parquet").iloc[0]
        table = factor_table_base(out, SAMPLE, "off_target").name
        written = any(p.name.startswith(table) for p in out.iterdir())
        assert bool(row["gc_correction_computed"]) == written, name
        computed[name] = bool(row["gc_correction_computed"])
    assert computed == {"dense": True, "sparse": False}


def test_extract_names_no_factor_table_it_did_not_write(sample, tmp_path, caplog):
    """The synthetic sample is too small to fit any GC factor, so no table is written."""
    out = tmp_path / "out"
    with caplog.at_level(logging.INFO):
        result = _extract(sample, out)
    assert result.exit_code == 0, result.output
    assert not any("correction_factors" in p.name for p in out.iterdir())
    claims = [
        r.getMessage() for r in caplog.records if "Correction factors" in r.getMessage()
    ]
    assert claims == [], f"logged a factor table that does not exist: {claims}"
    assert bool(record(out)["gc_correction_computed"]) is False


# ── 3. --chromosomes is not silently ignored ─────────────────────────────────────────────────


def test_extract_refuses_chromosomes_rather_than_ignoring_it(sample, tmp_path):
    """It was parsed and thrown away: `--chromosomes chr1` produced a whole-genome BED, and a
    sample marked complete. Refused as an unknown option, before anything is written."""
    out = tmp_path / "out"
    result = _extract(sample, out, "--chromosomes", "chr1")
    assert result.exit_code == 2, f"not refused as a usage error: {result.output}"
    assert "No such option" in result.output
    assert not (out / f"{SAMPLE}.bed.gz").exists()


def test_a_failed_gc_computation_is_warned_about(tmp_path, caplog):
    """A factor computation that raises (a genome without assets, an unreadable GC reference)
    records `gc_correction_computed=False`, correctly, and was logged at debug only: the run
    looked like one that had fitted nothing, with nothing at a level anyone reads."""
    from krewlyzer.core.sample_processor import write_extraction_outputs

    dense = {(b, gc): 10_000 for b in range(188) for gc in range(101)}
    with caplog.at_level(logging.INFO):
        write_extraction_outputs(
            _extraction(tmp_path, dense),
            tmp_path,
            genome="hg99",
            output_format="parquet",
        )
    row = pd.read_parquet(tmp_path / f"{SAMPLE}.metadata.parquet").iloc[0]
    assert bool(row["gc_correction_computed"]) is False
    warned = [
        r.getMessage()
        for r in caplog.records
        if r.levelno >= logging.WARNING and "GC" in r.getMessage()
    ]
    assert warned, "a failed GC computation left no warning"
