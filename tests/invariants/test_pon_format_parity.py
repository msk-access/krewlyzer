"""A PON-derived column must not depend on `--output-format`.

Three defects of this shape were live at once, and each was silent:

- **FSD** produced log2 ratios under `tsv` and raw counts under `parquet`.
  Two causes stacked: the caller guarded on `outputs.fsd.exists()` against a
  hardcoded `.tsv` name the writer no longer produced, and `read_table` is
  parquet-first, so the write-back re-read the stale pre-normalisation table.
- **`apply_fsc_gene_pon`** returned `None` with "file not found" under parquet
  and added no `depth_zscore`.
- **`apply_fsc_region_pon`** and the e1-only filter, identically.

None of them raised, and every column kept its name. 0.9.0 makes parquet the
Nextflow default, so all three were about to go from affecting nobody to
affecting every run.

The generalisable check is here rather than in each feature's own test file:
the defect is not in FSD or FSC, it is in the pattern of naming one extension
while the writer honours a flag.
"""

from __future__ import annotations

from pathlib import Path

from conftest import is_hydrated

import numpy as np
import pandas as pd
import pytest

from krewlyzer.core.fsc_processor import apply_fsc_gene_pon, apply_fsc_region_pon
from krewlyzer.core.fsd_processor import process_fsd
from krewlyzer.core.output_utils import read_table
from krewlyzer.pon.model import PonModel

_ROOT = Path(__file__).resolve().parents[2]
_PON = _ROOT / "src/krewlyzer/data/pon/GRCh37/all_unique/xs1.all_unique.pon.parquet"

#: Every way a caller can ask for output. `both` matters as much as the others
#: -- it leaves a TSV *and* a parquet on disk, which is what let the
#: parquet-first reader pick the wrong one.
FORMATS = ["tsv", "parquet", "both"]

pytestmark = pytest.mark.invariant


@pytest.fixture(scope="module")
def pon():
    # `is_hydrated`, not `.exists()`: a checkout without LFS leaves a 130-byte
    # pointer here, which exists happily and then dies inside pyarrow with
    # "Parquet magic bytes not found in footer". That is precisely the
    # guarded-on-exists() defect this module was written about.
    if not is_hydrated(_PON):
        pytest.skip("bundled PON not hydrated (git lfs pull)")
    model = PonModel.load(_PON)
    if not model.fsc_gene_baseline:
        pytest.skip("bundled PON has no fsc_gene_baseline")
    return model


def _written(directory: Path, stem: str) -> pd.DataFrame:
    """Read whichever file the format actually produced."""
    frame = read_table(directory / f"{stem}.tsv")
    assert frame is not None, f"nothing readable at {stem}"
    return frame


# ---------------------------------------------------------------------------
# FSD -- the PON adds `{bin}_logR` columns beside the raw counts
# ---------------------------------------------------------------------------


def _fsd_frame(pon) -> pd.DataFrame:
    arms = list(pon.fsd_baseline.arms)[:4] if pon.fsd_baseline else []
    if not arms:
        pytest.skip("bundled PON has no fsd_baseline")
    bins = [f"{s}-{s + 4}" for s in range(65, 400, 5)]
    rng = np.random.default_rng(0)
    data: dict = {"region": arms}
    for i, b in enumerate(bins):
        data[b] = rng.integers(50, 500, size=len(arms)).astype(float) + i
    frame = pd.DataFrame(data)
    frame["total"] = frame[bins].sum(axis=1)
    return frame


def test_fsd_logratios_are_identical_across_formats(tmp_path, pon):
    """`{bin}_logR` must be present, populated and equal whichever format ran.

    Measured on a real run: 0 `_logR` columns under parquet before the fix,
    67 after, and 67 under tsv throughout.
    """
    results = {}
    for fmt in FORMATS:
        d = tmp_path / fmt
        d.mkdir()
        _fsd_frame(pon).to_csv(d / "s.FSD.tsv", sep="\t", index=False)
        process_fsd(d / "s.FSD.tsv", pon_parquet_path=_PON, output_format=fmt)
        frame = _written(d, "s.FSD")
        logr = sorted(c for c in frame.columns if c.endswith("_logR"))
        assert logr, (
            f"--output-format {fmt} produced no _logR columns: the PON was "
            "silently not applied."
        )
        results[fmt] = (logr, frame[logr].to_numpy(float))

    names, baseline = results["tsv"]
    assert np.isfinite(baseline).any(), "all log-ratios are NaN; parity is vacuous"
    for fmt, (fmt_names, values) in results.items():
        assert (
            fmt_names == names
        ), f"--output-format {fmt} produced a different set of _logR columns"
        assert np.allclose(
            baseline, values, equal_nan=True
        ), f"FSD log-ratios under --output-format {fmt} differ from tsv."


# ---------------------------------------------------------------------------
# FSC gene and region -- a column is added, so parity is about its presence
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS)
def test_fsc_gene_zscore_appears_under_every_format(tmp_path, pon, fmt):
    genes = list(pon.fsc_gene_baseline.data)[:5]
    d = tmp_path / fmt
    d.mkdir()
    frame = pd.DataFrame(
        {"gene": genes, "normalized_depth": [100.0, 200.0, 300.0, 400.0, 500.0]}
    )
    # Write only what this format would have produced -- reproducing the
    # condition, not just the flag.
    if fmt == "tsv":
        frame.to_csv(d / "s.FSC.gene.tsv", sep="\t", index=False)
    elif fmt == "parquet":
        frame.to_parquet(d / "s.FSC.gene.parquet")
    else:
        frame.to_csv(d / "s.FSC.gene.tsv", sep="\t", index=False)
        frame.to_parquet(d / "s.FSC.gene.parquet")

    apply_fsc_gene_pon(d / "s.FSC.gene.tsv", pon, output_format=fmt)
    out = _written(d, "s.FSC.gene")
    assert "depth_zscore" in out.columns, f"no depth_zscore under {fmt}"
    assert out["depth_zscore"].notna().any(), f"depth_zscore all NaN under {fmt}"


@pytest.mark.parametrize("fmt", FORMATS)
def test_fsc_region_zscore_appears_under_every_format(tmp_path, pon, fmt):
    if not pon.fsc_region_baseline:
        pytest.skip("bundled PON has no fsc_region_baseline")
    keys = list(pon.fsc_region_baseline.data)[:5]
    parts = [k.replace(":", "-").split("-") for k in keys]
    frame = pd.DataFrame(
        {
            "chrom": [p[0] for p in parts],
            "start": [int(p[1]) for p in parts],
            "end": [int(p[2]) for p in parts],
            "normalized_depth": [100.0, 200.0, 300.0, 400.0, 500.0],
        }
    )
    d = tmp_path / fmt
    d.mkdir()
    if fmt == "tsv":
        frame.to_csv(d / "s.FSC.regions.tsv", sep="\t", index=False)
    elif fmt == "parquet":
        frame.to_parquet(d / "s.FSC.regions.parquet")
    else:
        frame.to_csv(d / "s.FSC.regions.tsv", sep="\t", index=False)
        frame.to_parquet(d / "s.FSC.regions.parquet")

    apply_fsc_region_pon(d / "s.FSC.regions.tsv", pon, output_format=fmt)
    out = _written(d, "s.FSC.regions")
    assert "depth_zscore" in out.columns, f"no depth_zscore under {fmt}"
    assert out["depth_zscore"].notna().any(), f"depth_zscore all NaN under {fmt}"
