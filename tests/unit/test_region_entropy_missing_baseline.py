"""Regression test: PON without a matching baseline table must not crash."""

import pandas as pd

from krewlyzer.core import region_entropy_processor as rep


def test_missing_baseline_degrades_to_zero_zscores(tmp_path, monkeypatch):
    """apply_pon_zscore returns early without writing when the PON lacks the table.

    Previously process_region_entropy() then called load_entropy_tsv() on the
    missing output, whose `assert df is not None` fired, and the raw file was
    subsequently unlinked by the caller -- losing the entropy results entirely.
    """
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
    assert "z_score" in df.columns
    assert df["z_score"].eq(0.0).all()
    assert df["label"].tolist() == ["CTCF"]
