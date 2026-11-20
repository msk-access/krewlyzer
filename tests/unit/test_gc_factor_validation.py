"""User-supplied GC factor tables are validated by the loader that will read them (redesign item 0).

`FileSchema.GC_FACTORS_TSV` described a three-column text file (`length_bin gc_pct factor`) that
no writer has produced for a long time, and it opened every file as text -- so the single-tool
`region-entropy` and `validate` commands rejected the six-column `.parquet` tables run-all
writes, and accepted the three-column shape the Rust loader skipped row by row. Validation now
asks the Rust loader itself, so "valid" means "loadable", in every format.
"""

import gzip

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from krewlyzer.core.asset_validation import FileSchema, validate_file

COLUMNS = [
    "length_bin_min",
    "length_bin_max",
    "gc_percent",
    "observed",
    "expected",
    "correction_factor",
]
ROWS = [(170, 175, 50, 100, 90, 1.25), (170, 175, 30, 80, 95, 0.75)]


def _tsv(path, gz=False):
    text = pd.DataFrame(ROWS, columns=COLUMNS).to_csv(sep="\t", index=False)
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    return path


def _parquet(path):
    schema = pa.schema(
        [
            ("length_bin_min", pa.uint64()),
            ("length_bin_max", pa.uint64()),
            ("gc_percent", pa.uint8()),
            ("observed", pa.uint64()),
            ("expected", pa.uint64()),
            ("correction_factor", pa.float64()),
        ]
    )
    pq.write_table(
        pa.Table.from_pandas(
            pd.DataFrame(ROWS, columns=COLUMNS), schema=schema, preserve_index=False
        ),
        path,
    )
    return path


@pytest.mark.unit
@pytest.mark.rust
@pytest.mark.parametrize("make", ["tsv", "tsv.gz", "parquet"])
def test_every_format_run_all_writes_validates(tmp_path, make):
    path = tmp_path / f"s.correction_factors.{make}"
    {
        "tsv": lambda: _tsv(path),
        "tsv.gz": lambda: _tsv(path, gz=True),
        "parquet": lambda: _parquet(path),
    }[make]()
    assert validate_file(path, FileSchema.GC_FACTORS_TSV) == len(ROWS)


@pytest.mark.unit
@pytest.mark.rust
def test_the_old_three_column_shape_is_refused(tmp_path):
    path = tmp_path / "old.tsv"
    path.write_text("length_bin\tgc_pct\tfactor\n22\t50\t2.0\n")
    with pytest.raises(ValueError, match="(?i)correction factor"):
        validate_file(path, FileSchema.GC_FACTORS_TSV)
