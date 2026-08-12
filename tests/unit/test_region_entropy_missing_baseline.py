"""A PON without a matching baseline must not crash -- and must not invent a z.

Two regressions, both about what happens when the baseline is missing.

The original one: `apply_pon_zscore` returns early without writing, then
`load_entropy_tsv` asserted on the missing output and the caller unlinked the
raw file -- losing the entropy results entirely.

The second, found later: the degraded path wrote `z_score = 0.0`. Zero is not
a neutral placeholder here, it is the most confident claim the column can
make -- "this sample sits exactly at the healthy baseline" -- and a reader has
no way to tell it apart from a measured zero. NaN is the honest value.
"""

import numpy as np
import pandas as pd

from krewlyzer.core import region_entropy_processor as rep


def test_missing_baseline_still_writes_the_entropy_results(tmp_path, monkeypatch):
    """The original regression: the output must survive a missing baseline."""
    raw = tmp_path / "sample.TFBS.raw.tsv"
    raw.write_text("label\tcount\tmean_size\tentropy\nCTCF\t1234\t167.2\t5.2300\n")
    out = tmp_path / "sample.TFBS.tsv"

    class _FakeRegionEntropy:
        @staticmethod
        def apply_pon_zscore(raw_path, pon_path, output_path, baseline_table):
            # Mimic the Rust early-return: no output file is created.
            return 0

    fake_core = type("_C", (), {"region_entropy": _FakeRegionEntropy})
    monkeypatch.setitem(__import__("sys").modules, "krewlyzer", __import__("krewlyzer"))
    monkeypatch.setattr("krewlyzer._core", fake_core, raising=False)

    n = rep.process_region_entropy(
        raw_path=raw,
        output_path=out,
        pon_parquet_path=tmp_path / "empty.pon.parquet",
        baseline_table="tfbs_baseline",
    )

    assert n == 0
    assert out.exists(), "output must still be written when the baseline is absent"
    df = pd.read_csv(out, sep="\t")
    assert df["label"].tolist() == ["CTCF"]
    assert df["entropy"].tolist() == [5.23], "the measurement itself must survive"


def test_a_missing_baseline_yields_nan_not_zero(tmp_path, monkeypatch):
    """Zero would claim the sample sits exactly at the healthy baseline.

    That is the most confident statement this column can make, and it would be
    made on the strength of having no baseline at all. It is also
    indistinguishable from a measured zero. NaN says what happened.
    """
    raw = tmp_path / "sample.TFBS.raw.tsv"
    raw.write_text("label\tcount\tmean_size\tentropy\nCTCF\t1234\t167.2\t5.2300\n")
    out = tmp_path / "sample.TFBS.tsv"

    class _FakeRegionEntropy:
        @staticmethod
        def apply_pon_zscore(raw_path, pon_path, output_path, baseline_table):
            return 0

    monkeypatch.setattr(
        "krewlyzer._core",
        type("_C", (), {"region_entropy": _FakeRegionEntropy}),
        raising=False,
    )

    rep.process_region_entropy(
        raw_path=raw,
        output_path=out,
        pon_parquet_path=tmp_path / "empty.pon.parquet",
        baseline_table="tfbs_baseline",
    )
    z = pd.read_csv(out, sep="\t")["z_score"]
    assert z.isna().all(), f"expected NaN, got {z.tolist()}"
    assert not z.eq(0.0).any()


def test_no_pon_at_all_also_yields_nan(tmp_path):
    """Same reasoning, the commoner path: nobody passed --pon-model."""
    raw = tmp_path / "sample.ATAC.raw.tsv"
    raw.write_text("label\tcount\tmean_size\tentropy\nBRCA\t900\t166.0\t4.9000\n")
    out = tmp_path / "sample.ATAC.tsv"

    rep.process_region_entropy(
        raw_path=raw,
        output_path=out,
        pon_parquet_path=None,
        baseline_table="atac_baseline",
    )
    df = pd.read_csv(out, sep="\t")
    assert np.isnan(df["z_score"]).all(), f"expected NaN, got {df['z_score'].tolist()}"
