"""A kept ``build-pon`` cache must be readable by ``build-pon --from-outputs``.

`--keep-sample-outputs` exists so a rebuild does not re-read every BAM: 55-97
minutes per sample, roughly 15 hours for a 47-sample cohort. Every defect found
in the 0.9.0 models has been in *aggregation*, so re-checking one should cost
minutes.

It could not. `process_sample` never called `write_motif_outputs` -- `run-all`
did it at its own call site and `build-pon` did not -- so a kept cache had no
``EndMotif`` or ``MDS`` table, even though the same k-mer counts reached the
model through memory. The one input `--from-outputs` cannot reconstruct was the
one input missing, and `krewlyzer motif` needs a BAM, so it could not be
backfilled either.

Measured against the real 0.9.0 cluster cache before the fix: all 47 xs1 duplex
directories were refused, each with ``missing .EndMotif``.

The resume story is deliberately "re-aggregate the cache", not "skip finished
samples in the extraction loop". The in-process path collects from live
`SampleOutputs` objects rather than from the directory, so skipping extraction
would skip collection too and drop the sample from every baseline without
saying so -- the exact class of failure this release exists to remove.
"""

from __future__ import annotations

import pytest

from krewlyzer.pon import from_outputs


@pytest.fixture
def processed_sample(tmp_path, real_bam, real_reference):
    """One sample through `process_sample`, as `build-pon` invokes it."""
    from krewlyzer.core.sample_processor import SampleParams, process_sample

    out = tmp_path / "DONOR1"
    out.mkdir(parents=True)
    process_sample(
        input_path=real_bam,
        output_dir=out,
        sample_name="DONOR1",
        reference=real_reference,
        params=SampleParams(threads=1),
        enable_fsc=True,
        enable_fsr=False,
        enable_fsd=True,
        enable_wps=True,
        enable_ocf=True,
        pon_mode=True,
        output_format="tsv",
        compress=False,
        write_motif_files=True,
    )
    return out


def test_the_motif_tables_are_written(processed_sample):
    """The two files a kept cache used to lack."""
    stem = "DONOR1"
    for name in ("EndMotif", "MDS"):
        assert from_outputs.resolve_table_path(
            processed_sample / f"{stem}.{name}"
        ), f"{name} table absent; the cache is not re-aggregatable"


def test_the_end_motif_table_holds_real_frequencies(processed_sample):
    """Present is not enough -- it has to carry the k-mer distribution."""
    import pandas as pd

    path = from_outputs.resolve_table_path(processed_sample / "DONOR1.EndMotif")
    frame = pd.read_csv(path, sep="\t")
    assert {"Motif", "Frequency"}.issubset(frame.columns)
    assert len(frame) == 256, "expected the full 4-mer alphabet"
    assert frame["Frequency"].sum() == pytest.approx(1.0, abs=1e-3)
    # ...and it varies. A flat 1/256 everywhere would mean nothing was counted.
    assert frame["Frequency"].nunique() > 1


def test_the_collector_reconstructs_the_motif_record(processed_sample):
    """The round trip that `--from-outputs` depends on."""
    record = from_outputs._motif_record(processed_sample, "DONOR1")
    assert record is not None, "the collector could not read what we just wrote"
    assert len(record["kmers"]) == 256
    assert 0.0 < record["mds"] <= 1.0


def test_writing_them_is_off_by_default(tmp_path, real_bam, real_reference):
    """`run-all` writes them at its own call site; twice would be waste.

    Also pins that the flag is what changed the behaviour, rather than the
    files having appeared for some unrelated reason.
    """
    from krewlyzer.core.sample_processor import SampleParams, process_sample

    out = tmp_path / "DONOR2"
    out.mkdir(parents=True)
    process_sample(
        input_path=real_bam,
        output_dir=out,
        sample_name="DONOR2",
        reference=real_reference,
        params=SampleParams(threads=1),
        enable_fsc=True,
        enable_fsr=False,
        enable_fsd=False,
        enable_wps=False,
        enable_ocf=False,
        pon_mode=True,
        output_format="tsv",
        compress=False,
    )
    assert from_outputs.resolve_table_path(out / "DONOR2.EndMotif") is None
