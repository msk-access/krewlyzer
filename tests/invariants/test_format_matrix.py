"""The output format must not change the numbers.

`test_serializer_output_formats.py` covers the serializer's *discovery* across
formats, and the CLI parity test covers whether the flags exist. Neither runs
the pipeline twice and compares values, which is the assumption everything else
rests on: that `--output-format parquet` is a serialisation choice and not a
different computation.

Worth asserting because the two paths genuinely differ -- Parquet is written
directly, TSV goes through a text round-trip with `float_format="%.4f"`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import pandas as pd
import pytest

from krewlyzer.core.output_utils import read_table

from .pipeline import run_pipeline
from .synth import BASELINE

pytestmark = [pytest.mark.invariant, pytest.mark.rust, pytest.mark.slow]

# %.4f on the TSV path, so values agree to four decimals and no further.
TSV_DECIMALS = 4
TOLERANCE = 0.5 * 10**-TSV_DECIMALS


@pytest.fixture(scope="module")
def by_format(tmp_path_factory) -> Dict[str, Path]:
    """Run the same input through each output format."""
    root = tmp_path_factory.mktemp("formats")
    produced = {}
    for fmt, compress in (("tsv", False), ("tsv", True), ("parquet", False)):
        label = f"{fmt}{'_gz' if compress else ''}"
        run = run_pipeline(root / label, BASELINE, output_format=fmt, compress=compress)
        produced[label] = run.output_dir
    return produced


def _load(directory: Path, stem: str) -> pd.DataFrame:
    for suffix in (".tsv", ".tsv.gz", ".parquet"):
        candidate = directory / f"{stem}{suffix}"
        if candidate.exists():
            frame = read_table(candidate)
            assert frame is not None, f"could not read {candidate}"
            return frame
    raise AssertionError(
        f"no {stem} table in {directory} (found: "
        f"{sorted(p.name for p in directory.iterdir())})"
    )


@pytest.mark.parametrize("stem", ["fsc_counts", "FSD"])
def test_values_are_identical_across_formats(by_format, stem):
    """Same input, three serialisations, one set of numbers."""
    name = f"{BASELINE.name}.{stem}"
    frames = {label: _load(d, name) for label, d in by_format.items()}

    reference_label = "tsv"
    reference = frames[reference_label]

    for label, frame in frames.items():
        if label == reference_label:
            continue
        assert list(frame.columns) == list(
            reference.columns
        ), f"{stem}: {label} has different columns from {reference_label}"
        assert len(frame) == len(
            reference
        ), f"{stem}: {label} has a different row count"

        for column in reference.select_dtypes(include="number").columns:
            delta = (frame[column] - reference[column]).abs().max()
            assert delta <= TOLERANCE, (
                f"{stem}.{column}: {label} differs from {reference_label} by "
                f"{delta:g}, more than the {TSV_DECIMALS}-decimal TSV rounding "
                f"allows ({TOLERANCE:g})"
            )


def test_each_format_actually_produced_its_own_extension(by_format):
    """Guards the fixture: if every run silently wrote TSV, the comparison
    above would pass while testing nothing."""
    expected = {"tsv": ".tsv", "tsv_gz": ".tsv.gz", "parquet": ".parquet"}
    for label, directory in by_format.items():
        names = [p.name for p in directory.iterdir()]
        suffix = expected[label]
        assert any(
            n.endswith(suffix) for n in names
        ), f"{label}: no {suffix} file was written (found {sorted(names)})"
