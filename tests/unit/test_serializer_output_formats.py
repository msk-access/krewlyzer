"""features.json must be format-independent (.tsv / .tsv.gz / .parquet)."""

from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.core.feature_serializer import FeatureSerializer, _resolve_output_path

SID = "S1"

TABLES = {
    "FSD": pd.DataFrame([{"region": "chr1:0-125000000", "65-69": 1.0, "70-74": 2.0, "total": 3.0}]),
    "FSR": pd.DataFrame([{"region": "chr1", "core_short_long_ratio": 1.4}]),
    "FSC": pd.DataFrame([{"chrom": "1", "start": 0, "end": 100, "core_short": 5.0, "total": 9.0}]),
    "FSC.gene": pd.DataFrame([{"gene": "TP53", "total": 10.0}]),
    "FSC.regions": pd.DataFrame([{"gene": "TP53", "region_name": "e1", "total": 4.0}]),
    "FSC.regions.e1only": pd.DataFrame([{"gene": "TP53", "region_name": "e1", "total": 4.0}]),
    "EndMotif": pd.DataFrame([{"AAAA": 0.01, "CCCC": 0.02}]),
    "MDS": pd.DataFrame([{"Sample": SID, "MDS": 0.97}]),
    "OCF": pd.DataFrame([{"tissue": "Liver", "OCF": 265.3}]),
    "TFBS": pd.DataFrame([{"label": "CTCF", "count": 100.0, "entropy": 6.9, "z_score": 1.2}]),
    "ATAC": pd.DataFrame([{"label": "LIHC", "count": 80.0, "entropy": 6.5, "z_score": 2.1}]),
    "metadata": pd.DataFrame([{"genome": "hg19", "assay": "xs2"}]),
}


def _write(directory: Path, fmt: str) -> None:
    for name, df in TABLES.items():
        if fmt == "tsv":
            df.to_csv(directory / f"{SID}.{name}.tsv", sep="\t", index=False)
        elif fmt == "tsv.gz":
            df.to_csv(directory / f"{SID}.{name}.tsv.gz", sep="\t",
                      index=False, compression="gzip")
        elif fmt == "parquet":
            df.to_parquet(directory / f"{SID}.{name}.parquet", index=False)


@pytest.mark.parametrize("fmt", ["tsv", "tsv.gz", "parquet"])
def test_from_outputs_is_format_independent(tmp_path, fmt):
    """Regression: compressed and Parquet runs produced an EMPTY features.json.

    Every probe in from_outputs() was `Path(f"{sample}.FOO.tsv").exists()`, and
    that gate ran *before* read_table() (which does understand .gz/.parquet).
    So with --compress or --output-format parquet every feature was silently
    skipped and --generate-json emitted an empty payload.
    """
    out = tmp_path / fmt.replace(".", "_")
    out.mkdir()
    _write(out, fmt)

    serializer = FeatureSerializer.from_outputs(SID, out)

    expected = {
        "fsd", "fsr", "fsc", "fsc_gene", "fsc_region", "fsc_region_e1",
        "motif", "ocf", "tfbs", "atac",
    }
    assert expected.issubset(set(serializer.features)), (
        f"{fmt}: missing {expected - set(serializer.features)}"
    )
    assert serializer.metadata.get("assay") == "xs2", f"{fmt}: metadata not loaded"


def test_all_formats_agree(tmp_path):
    """The three formats must produce the same feature set."""
    seen = {}
    for fmt in ("tsv", "tsv.gz", "parquet"):
        out = tmp_path / fmt.replace(".", "_")
        out.mkdir()
        _write(out, fmt)
        seen[fmt] = set(FeatureSerializer.from_outputs(SID, out).features)
    assert seen["tsv"] == seen["tsv.gz"] == seen["parquet"], seen


def test_resolver_prefers_tsv_then_gz_then_parquet(tmp_path):
    assert _resolve_output_path(tmp_path, SID, "FSD") is None

    gz = tmp_path / f"{SID}.FSD.tsv.gz"
    gz.write_bytes(b"")
    assert _resolve_output_path(tmp_path, SID, "FSD") == gz

    plain = tmp_path / f"{SID}.FSD.tsv"
    plain.write_text("")
    assert _resolve_output_path(tmp_path, SID, "FSD") == plain


# ---------------------------------------------------------------------------
# E1-only FSC naming across output formats
# ---------------------------------------------------------------------------

FSC_REGIONS = pd.DataFrame(
    [
        {"chrom": "1", "start": 100, "end": 200, "gene": "TP53", "region_name": "t1",
         "region_bp": 100, "ultra_short": 1.0, "core_short": 2.0, "mono_nucl": 3.0,
         "di_nucl": 1.0, "long": 1.0, "total": 8.0, "normalized_depth": 10.0},
        {"chrom": "1", "start": 300, "end": 400, "gene": "TP53", "region_name": "t2",
         "region_bp": 100, "ultra_short": 1.0, "core_short": 2.0, "mono_nucl": 3.0,
         "di_nucl": 1.0, "long": 1.0, "total": 8.0, "normalized_depth": 10.0},
    ]
)


@pytest.mark.parametrize(
    "fmt,compress,in_name,out_name",
    [
        ("tsv", False, "S1.FSC.regions.tsv", "S1.FSC.regions.e1only.tsv"),
        ("tsv", True, "S1.FSC.regions.tsv.gz", "S1.FSC.regions.e1only.tsv.gz"),
        ("parquet", False, "S1.FSC.regions.parquet", "S1.FSC.regions.e1only.parquet"),
    ],
)
def test_e1_output_name_is_correct_for_every_format(
    tmp_path, fmt, compress, in_name, out_name
):
    """Regression: with .tsv.gz input the E1 file was named
    'S1.FSC.regions.tsv.e1only.tsv.gz'.

    Path.stem strips only the last dot-segment, so for 'S1.FSC.regions.tsv.gz'
    it returns 'S1.FSC.regions.tsv'; the '.regions' test then missed and
    '.e1only' was appended to a stem still containing '.tsv'.
    """
    from krewlyzer.core.fsc_processor import filter_fsc_to_e1

    src = tmp_path / in_name
    if fmt == "parquet":
        FSC_REGIONS.to_parquet(src, index=False)
    elif compress:
        FSC_REGIONS.to_csv(src, sep="\t", index=False, compression="gzip")
    else:
        FSC_REGIONS.to_csv(src, sep="\t", index=False)

    produced = filter_fsc_to_e1(src, output_format=fmt, compress=compress)

    assert produced is not None
    assert produced.name == out_name, f"got {produced.name}, want {out_name}"
    assert produced.exists()


def test_strip_table_extension_handles_compound_suffixes():
    from krewlyzer.core.output_utils import strip_table_extension

    assert strip_table_extension("s.FSC.regions.tsv") == "s.FSC.regions"
    assert strip_table_extension("s.FSC.regions.tsv.gz") == "s.FSC.regions"
    assert strip_table_extension("s.FSC.regions.parquet") == "s.FSC.regions"
    assert strip_table_extension("s.FSC.regions") == "s.FSC.regions"
