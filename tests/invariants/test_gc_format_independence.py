"""GC correction must not depend on the output format (redesign item 0; invariant #6).

Production ran `--output-format both --compress`, which writes `.correction_factors.tsv.gz` and
`.parquet` but never a plain `.tsv`. Three Python consumers looked only for the plain `.tsv`, so in
0.9.x cohorts panel WPS, every TFBS and ATAC table, and gene FSC ran **uncorrected** -- while PON
donors, built with TSV forced, were corrected. Nothing failed and nothing said so.

This drives the real Rust writer on a GC-biased synthetic sample in each output layout and
asserts the run resolves, in every layout, a table carrying *identical* factors -- and factors that
are not all 1, so the test cannot pass because nothing is corrected (invariant #1). The writer
keeps full precision in TSV, so the layouts agree exactly.
"""

import gzip
import math

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from krewlyzer import _core
from krewlyzer.core.gc_assets import resolve_feature_gc_tables

pytestmark = [pytest.mark.invariant, pytest.mark.rust]

LAYOUTS = {  # name -> (output_format, compress)
    "tsv": ("tsv", False),
    "parquet": ("parquet", False),
    "production": ("both", True),  # what the 0.9.x cohorts ran
}


def _reference(path):
    """Genome-shaped expected counts for one nucleosome-length bin (LengthBin 22, 170-174 bp)."""
    rows = []
    for gc in range(20, 81):
        z = (gc / 100 - 0.43) / 0.08
        rows.append((170, 175, gc, int(round(2e6 * math.exp(-0.5 * z * z))) + 100))
    table = pa.table(
        {
            "length_bin_min": pa.array([r[0] for r in rows], pa.uint32()),
            "length_bin_max": pa.array([r[1] for r in rows], pa.uint32()),
            "gc_percent": pa.array([r[2] for r in rows], pa.uint8()),
            "expected_count": pa.array([r[3] for r in rows], pa.uint64()),
        }
    )
    pq.write_table(table, path)
    return path, {r[2]: r[3] for r in rows}


def _observations(expected):
    """A sample whose coverage rises with GC: observed = expected x exp(6 (gc - 0.5)) / 50."""
    return {
        (22, gc): max(60, int(e * math.exp(6 * (gc / 100 - 0.5)) / 50))
        for gc, e in expected.items()
    }


def _read(path):
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(
        path, sep="\t", compression="gzip" if path.suffix == ".gz" else None
    )


def _run_layout(tmp_path, name, ref, obs, valid):
    fmt, compress = LAYOUTS[name]
    d = tmp_path / name
    d.mkdir()
    bed = d / "s.bed.gz"
    with gzip.open(bed, "wt") as f:
        f.write("chr1\t100\t270\t0.5\n")
    n = _core.gc.compute_and_write_gc_factors(
        obs,
        str(ref),
        str(valid),
        str(d / "s.correction_factors.tsv"),
        output_format=fmt,
        compress=compress,
    )
    assert n > 0
    return d, bed


def test_every_output_layout_resolves_the_same_factors(tmp_path):
    ref, expected = _reference(tmp_path / "ref.parquet")
    valid = tmp_path / "valid.bed"
    valid.write_text("chr1\t0\t1000\n")
    obs = _observations(expected)

    tables = {}
    for name in LAYOUTS:
        d, bed = _run_layout(tmp_path, name, ref, obs, valid)
        off, _ = resolve_feature_gc_tables(bed, gc_correct=True, panel_mode=False)
        assert off.found, f"{name}: no factor table resolved ({off.reason})"
        frame = (
            _read(off.path)
            .sort_values(["length_bin_min", "gc_percent"])
            .reset_index(drop=True)
        )
        tables[name] = frame

    # The production layout writes no plain .tsv -- the file the old lookups required.
    assert not (tmp_path / "production" / "s.correction_factors.tsv").exists()

    base = tables["tsv"]
    for name, frame in tables.items():
        pd.testing.assert_frame_equal(
            frame.astype({"correction_factor": float}),
            base.astype({"correction_factor": float}),
            check_dtype=False,
            obj=f"{name} vs tsv",
        )
    # Two inputs, two outputs: a GC-biased sample must get factors that are not all 1.
    assert (base["correction_factor"] - 1).abs().max() > 0.1


def test_gc_correction_off_resolves_nothing(tmp_path):
    ref, expected = _reference(tmp_path / "ref.parquet")
    valid = tmp_path / "valid.bed"
    valid.write_text("chr1\t0\t1000\n")
    d, bed = _run_layout(tmp_path, "production", ref, _observations(expected), valid)
    off, on = resolve_feature_gc_tables(bed, gc_correct=False, panel_mode=True)
    assert not off.found and not on.found
    assert off.record() == "none:gc correction disabled"


@pytest.mark.parametrize("population", ["off_target", "on_target"])
def test_an_unreadable_table_stops_the_run_instead_of_being_recorded_as_applied(
    tmp_path, population
):
    """A table that exists but cannot be read used to count as found: gene FSC, panel WPS, TFBS
    and ATAC each failed on it, each failure was caught and logged, and the metadata still
    recorded `gc_correction_applied_*=True`. Review of PR #89."""
    ref, expected = _reference(tmp_path / "ref.parquet")
    valid = tmp_path / "valid.bed"
    valid.write_text("chr1\t0\t1000\n")
    d, bed = _run_layout(tmp_path, "production", ref, _observations(expected), valid)
    ontarget = d / "s.correction_factors.ontarget.parquet"
    ontarget.write_bytes((d / "s.correction_factors.parquet").read_bytes())
    resolve_feature_gc_tables(bed, gc_correct=True, panel_mode=True)

    broken = (
        ontarget if population == "on_target" else d / "s.correction_factors.parquet"
    )
    broken.write_bytes(b"not a parquet file")
    label = population.replace("_", "-")
    with pytest.raises(
        RuntimeError, match=f"{label} GC factor table .* cannot be read"
    ):
        resolve_feature_gc_tables(bed, gc_correct=True, panel_mode=True)
