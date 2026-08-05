"""Two Rust baseline readers parse plain TSV. They must be given plain TSV.

``compute_fsd_baseline`` and ``compute_region_mds_baseline`` read with
``BufReader::lines()`` -- no gzip, no parquet. ``File::open`` succeeds on both
of those anyway, the header parse then finds no usable columns, and every
sample is skipped. Three samples in, zero out, exit 0.

That is not hypothetical. A ``run-all`` directory writes ``.parquet`` and
``.tsv.gz`` and *no* plain ``.tsv``, so it is the normal case for the very
layout ``--from-outputs`` exists to read. Measured against the real 0.8.3
corpus: ``region_mds`` came back empty, logged one warning, and the build
carried on to write a PON with no ``region_mds`` block at all.

Nothing downstream would have said so. ``validate-pon`` iterates the blocks
present in the file and skips the empty ones, so a block that vanished
entirely is indistinguishable from one that was never expected.

The fix normalises inputs where the constraint lives, in the two callers that
know about it, rather than at each call site.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.pon.build import (
    _as_plain_tsv,
    _compute_fsd_baseline,
    _compute_region_mds_baseline,
)

SIZE_BINS = [f"{s}-{s + 4}" for s in range(65, 105, 5)]


def _fsd_frame(scale: float) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"arm": arm, **{b: (i + 1) * scale for i, b in enumerate(SIZE_BINS)}}
            for arm in ("1p", "1q", "2p")
        ]
    )


def _write_every_format(frame: pd.DataFrame, base: Path) -> dict:
    """The same table as ``.tsv``, ``.tsv.gz`` and ``.parquet``.

    Extensions are appended, never `with_suffix`. On a compound base like
    ``s0.MDS.gene`` that call replaces ``.gene`` and yields ``s0.MDS.tsv`` --
    the trap `write_motif_outputs` documents, met here while writing this file.
    """
    plain = Path(f"{base}.tsv")
    frame.to_csv(plain, sep="\t", index=False)
    gz = Path(f"{base}.tsv.gz")
    with gzip.open(gz, "wt") as handle:
        frame.to_csv(handle, sep="\t", index=False)
    parquet = Path(f"{base}.parquet")
    frame.to_parquet(parquet)
    return {"tsv": plain, "tsv.gz": gz, "parquet": parquet}


@pytest.fixture
def fsd_formats(tmp_path):
    written = [
        _write_every_format(_fsd_frame(1.0 + i), tmp_path / f"s{i}.FSD")
        for i in range(3)
    ]
    return {fmt: [str(w[fmt]) for w in written] for fmt in ("tsv", "tsv.gz", "parquet")}


def test_the_premise_gzip_and_parquet_reach_rust_as_nothing(tmp_path, fsd_formats):
    """Guard the assumption this whole module rests on.

    If Rust ever learns to read these formats, the normalisation below becomes
    dead weight and this test says so by failing.
    """
    from krewlyzer import _core

    assert _core.pon_builder.compute_fsd_baseline(fsd_formats["tsv"]), "plain TSV works"
    assert not _core.pon_builder.compute_fsd_baseline(fsd_formats["tsv.gz"])
    assert not _core.pon_builder.compute_fsd_baseline(fsd_formats["parquet"])


@pytest.mark.parametrize("fmt", ["tsv", "tsv.gz", "parquet"])
def test_the_fsd_baseline_is_the_same_whatever_the_format(fmt, fsd_formats):
    """The point of the fix: output must not depend on how the input was stored."""
    reference = _compute_fsd_baseline([], fsd_formats["tsv"])
    baseline = _compute_fsd_baseline([], fsd_formats[fmt])

    assert len(baseline.arms) == 3, f"{fmt} lost arms"
    assert sorted(baseline.arms) == sorted(reference.arms)
    for arm in reference.arms:
        assert baseline.arms[arm]["expected"] == reference.arms[arm]["expected"]
        assert baseline.arms[arm]["std"] == reference.arms[arm]["std"]


def test_an_unreadable_fsd_cohort_says_so_instead_of_returning_empty(tmp_path):
    """The old wording blamed the backend for doing what it was asked.

    "no data returned from Rust" sent a reader hunting a backend bug once
    already, in `_compute_wps_baseline`. The message has to describe the
    inputs.
    """
    empty = tmp_path / "junk.FSD.tsv"
    empty.write_text("not\ta\tfsd\ttable\n")
    with pytest.raises(RuntimeError) as excinfo:
        _compute_fsd_baseline([], [str(empty)])
    message = str(excinfo.value)
    assert "Rust" not in message, "still blames the backend"
    assert "junk.FSD.tsv" in message, "does not name an input"


def test_plain_tsv_is_passed_through_untouched(tmp_path, fsd_formats):
    """The in-process path is the common case and must copy nothing."""
    readable, complaints = _as_plain_tsv(fsd_formats["tsv"], tmp_path / "staging")
    assert readable == fsd_formats["tsv"]
    assert not complaints


def test_other_formats_are_materialised_into_the_staging_directory(tmp_path):
    """And the original is left alone."""
    staging = tmp_path / "staging"
    staging.mkdir()
    source = tmp_path / "s.FSD.parquet"
    _fsd_frame(1.0).to_parquet(source)

    readable, complaints = _as_plain_tsv([str(source)], staging)
    assert not complaints
    assert len(readable) == 1
    written = Path(readable[0])
    assert written.parent == staging and written.suffix == ".tsv"
    assert source.exists(), "the input was consumed"
    assert pd.read_csv(written, sep="\t").equals(pd.read_parquet(source))


def test_an_unreadable_input_is_reported_not_dropped(tmp_path):
    """Silently shrinking the cohort is the failure being fixed, not the fix."""
    good = tmp_path / "good.FSD.tsv"
    _fsd_frame(1.0).to_csv(good, sep="\t", index=False)
    bad = tmp_path / "bad.FSD.parquet"
    bad.write_bytes(b"this is not parquet")

    readable, complaints = _as_plain_tsv([str(good), str(bad)], tmp_path / "st")
    assert len(readable) == 1
    assert len(complaints) == 1 and "bad.FSD.parquet" in complaints[0]


# ---------------------------------------------------------------------------
# region_mds -- the one that failed quietly
# ---------------------------------------------------------------------------


def _mds_gene_frame(offset: float) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "gene": g,
                "mds_mean": 0.84 + offset + i / 100,
                "mds_e1": 0.80 + offset + i / 100,
            }
            for i, g in enumerate(("TP53", "EGFR", "KRAS"))
        ]
    )


@pytest.fixture
def mds_gene_formats(tmp_path):
    written = [
        _write_every_format(_mds_gene_frame(i / 50), tmp_path / f"s{i}.MDS.gene")
        for i in range(3)
    ]
    return {fmt: [str(w[fmt]) for w in written] for fmt in ("tsv", "tsv.gz", "parquet")}


@pytest.mark.parametrize("fmt", ["tsv", "tsv.gz", "parquet"])
def test_region_mds_builds_from_any_format(fmt, mds_gene_formats):
    baseline = _compute_region_mds_baseline(mds_gene_formats[fmt])
    assert baseline is not None, f"{fmt} produced no baseline"
    assert len(baseline.gene_baseline) == 3


def test_an_empty_region_mds_result_raises_rather_than_warning(tmp_path):
    """It used to warn and return None, so the block simply vanished.

    `validate-pon` skips blocks that are absent from the file, so a PON missing
    `region_mds` entirely still exits 0 -- the block cannot be checked once it
    is gone. The build is the only place that can notice.
    """
    junk = tmp_path / "junk.MDS.gene.tsv"
    junk.write_text("gene\tsomething_else\nTP53\t1\n")
    with pytest.raises(RuntimeError, match="region-MDS baseline is empty"):
        _compute_region_mds_baseline([str(junk)])


def test_region_mds_does_not_fall_back_to_a_standard_normal(
    monkeypatch, mds_gene_formats
):
    """`data.get("mds_mean", 0.0), data.get("mds_std", 1.0)`.

    Mean 0 with sigma 1 makes z equal the raw MDS -- about 0.95, an entirely
    ordinary-looking z-score for a statistic that was never computed.

    The backend is stubbed deliberately. Rust always supplies these keys, so
    against real input the fallback never fires and asserting on the output
    would pin nothing -- a first version of this test passed just as happily
    with the `0.0 / 1.0` defaults restored. The only way to test a fallback is
    to make it fire.
    """
    from krewlyzer import _core

    monkeypatch.setattr(
        _core.pon_builder,
        "compute_region_mds_baseline",
        lambda paths: {"TP53": {"n_samples": 4}},
    )
    baseline = _compute_region_mds_baseline(mds_gene_formats["tsv"])
    stats = baseline.gene_baseline["TP53"]

    import math

    assert math.isnan(stats["mds_std"]), f"placeholder sigma: {stats['mds_std']}"
    assert math.isnan(stats["mds_mean"]), f"placeholder mean: {stats['mds_mean']}"
    assert math.isnan(stats["mds_e1_std"]) and math.isnan(stats["mds_e1_mean"])
