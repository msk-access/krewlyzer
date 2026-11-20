"""Gene FSC must actually be GC-corrected when given a factor table (redesign item 0).

Before this fix, run-all loaded the on-target table into a Python dict and `aggregate_by_gene`
wrote it back out as a three-column TSV. The Rust loader skips rows with fewer than six fields,
so it loaded zero factors and every fragment counted 1.0 -- in every sample and every PON donor,
while the log said "Using on-target GC factors". These tests pin that a table changes the counts,
in either format, and that a table the loader cannot read is an error, not a silent no-op.
"""

import gzip
from types import SimpleNamespace

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from krewlyzer.core.fsc_processor import aggregate_by_gene

COLUMNS = [
    "length_bin_min",
    "length_bin_max",
    "gc_percent",
    "observed",
    "expected",
    "correction_factor",
]
# 170-bp fragments at GC 0.50 fall in LengthBin 22 (170-174 bp), GC key 50.
ROWS = [(170, 175, 50, 100, 100, 2.0), (170, 175, 30, 100, 100, 0.5)]


def _bed(tmp_path):
    bed = tmp_path / "s.bed.gz"
    with gzip.open(bed, "wt") as f:
        f.write("chr1\t100\t270\t0.5\n")
        f.write("chr1\t110\t280\t0.5\n")
    genes = {
        "Gene1": [SimpleNamespace(chrom="chr1", start=0, end=500, name="Gene1_t1")]
    }
    return bed, genes


def _total(path):
    return float(pd.read_csv(path, sep="\t").set_index("gene").loc["Gene1", "total"])


def _write_tsv(path, rows):
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, sep="\t", index=False)
    return path


def _write_parquet(path, rows):
    df = pd.DataFrame(rows, columns=COLUMNS)
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
    pq.write_table(pa.Table.from_pandas(df, schema=schema, preserve_index=False), path)
    return path


@pytest.mark.unit
@pytest.mark.rust
@pytest.mark.parametrize("fmt", ["tsv", "parquet"])
def test_a_factor_table_changes_gene_counts(tmp_path, fmt):
    bed, genes = _bed(tmp_path)
    raw = aggregate_by_gene(
        bed, genes, tmp_path / "raw.FSC.gene.tsv", aggregate_by="gene"
    )
    write = _write_tsv if fmt == "tsv" else _write_parquet
    table = write(tmp_path / f"s.correction_factors.ontarget.{fmt}", ROWS)
    weighted = aggregate_by_gene(
        bed,
        genes,
        tmp_path / "gc.FSC.gene.tsv",
        correction_factors=table,
        aggregate_by="gene",
    )
    assert _total(raw) == 2.0
    # Two inputs, two outputs (invariant #1): the GC-50 factor is above 1, so the total rises.
    assert _total(weighted) > 2.0


@pytest.mark.unit
@pytest.mark.rust
def test_a_table_the_loader_cannot_read_is_an_error(tmp_path):
    # The three-column shape the old dict path wrote.
    bed, genes = _bed(tmp_path)
    bad = tmp_path / "s.correction_factors.ontarget.tsv"
    bad.write_text("length_bin\tgc_pct\tfactor\n22\t50\t2.0\n")
    with pytest.raises(Exception, match="(?i)correction factor"):
        aggregate_by_gene(
            bed,
            genes,
            tmp_path / "x.FSC.gene.tsv",
            correction_factors=bad,
            aggregate_by="gene",
        )


@pytest.mark.unit
def test_a_dict_of_factors_is_no_longer_accepted(tmp_path):
    bed, genes = _bed(tmp_path)
    with pytest.raises(TypeError, match="path"):
        aggregate_by_gene(
            bed,
            genes,
            tmp_path / "x.FSC.gene.tsv",
            correction_factors={(22, 50): 2.0},
            aggregate_by="gene",
        )
