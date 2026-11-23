"""A GC factor table exists exactly when the run says it wrote one.

The tools reading a sample's folder find its factor table by existence (`resolve_gc_factors`), so a
table left behind is used, and a log line naming a table that was not written sends whoever reads
the log after one that does not exist. Two ways that went wrong:

- With `-F both`, Rust wrote the TSV table and then failed on the Parquet one. The record said
  `gc_correction_computed=False` while a `.tsv` table sat in the folder for every tool to use.
- The on-target factor writes in run-all and `process_sample` logged "Wrote 0 ... factors" (and the
  off-target one in `process_sample` "Computed 0") when nothing was fitted and nothing written.
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
    from krewlyzer.core.sample_processor import SampleParams, process_sample

    out = tmp_path / "out"
    out.mkdir()
    with caplog.at_level(logging.DEBUG):
        process_sample(
            input_path=sample.bam,
            output_dir=out,
            sample_name=SAMPLE,
            reference=sample.reference,
            params=SampleParams(threads=1, exclude_regions=sample.assets.blacklist),
            target_regions=sample.assets.targets,
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
    assert not any("correction_factors" in p.name for p in out.iterdir())
    assert _claims(caplog) == []
    assert _warned_none(caplog, "off-target")
    assert _warned_none(caplog, "on-target")


def test_extract_names_no_on_target_table_it_did_not_write(sample, tmp_path, caplog):
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    out = tmp_path / "out"
    with caplog.at_level(logging.DEBUG):
        result = CliRunner().invoke(
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
                str(sample.assets.targets),
                "-F",
                "parquet",
            ],
        )
    assert result.exit_code == 0, result.output
    assert not any("correction_factors" in p.name for p in out.iterdir())
    assert _claims(caplog) == []
    assert _warned_none(caplog, "on-target")
