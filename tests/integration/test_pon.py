"""
Integration tests for PON (Panel of Normals) functionality.

Tests PON model loading, application, and z-score calculations.
"""

import pytest
import pandas as pd
from pathlib import Path


@pytest.fixture
def sample_pon_parquet(tmp_path):
    """A PON in the schema `build-pon` actually writes.

    The previous version used `version` / `genome` / `sample_count`, which no
    krewlyzer has ever produced -- the real metadata row carries
    `schema_version` / `assay` / `n_samples` / `reference`. `PonModel.load`
    reads through `meta.get(key, default)`, so every field fell back to its
    default and the fixture loaded as a completely empty model: `assay=''`,
    `n_samples=0`, every baseline `None`.

    Nothing caught it because the tests using it only asserted
    `model is not None`, and the one test that checked a real field was
    marked skip with "PON parquet schema requires production format" rather
    than the fixture being corrected.

    Column names verified against
    `src/krewlyzer/data/pon/GRCh37/all_unique/xs1.all_unique.pon.parquet`.
    """
    metadata = pd.DataFrame(
        {
            "table": ["metadata"],
            "schema_version": ["1.0"],
            "assay": ["xs1"],
            "build_date": ["2026-01-01"],
            "n_samples": [10.0],
            "reference": ["Homo_sapiens_assembly19"],
            "panel_mode": [True],
            "target_regions_file": ["xs1.targets.bed.gz"],
            # `load_pon_model` refuses anything below MIN_PON_VERSION, so a
            # fixture standing in for a *valid* PON has to carry a version.
            # Recording one is now part of what valid means.
            "krewlyzer_version": ["0.9.0"],
        }
    )

    gc_bias = pd.DataFrame(
        {
            "table": ["gc_bias"] * 3,
            "gc_bin": [0.3, 0.4, 0.5],
            "short_expected": [100.0, 120.0, 110.0],
            "short_std": [10.0, 12.0, 11.0],
            "intermediate_expected": [200.0, 220.0, 210.0],
            "intermediate_std": [20.0, 22.0, 21.0],
            "long_expected": [50.0, 55.0, 52.0],
            "long_std": [5.0, 5.5, 5.2],
        }
    )

    pon_file = tmp_path / "test.pon.parquet"
    pd.concat([metadata, gc_bias], ignore_index=True).to_parquet(pon_file, index=False)
    return pon_file


@pytest.mark.integration
def test_pon_model_loading(tmp_path, sample_pon_parquet):
    """The metadata row is parsed into the model, field by field.

    Unskipped: it was marked "PON parquet schema requires production format",
    which was true of the fixture, not of the loader.
    """
    from krewlyzer.pon.model import PonModel

    # Load model
    model = PonModel.load(sample_pon_parquet)

    assert model is not None
    assert model.assay == "xs1"
    assert model.n_samples == 10
    assert model.reference == "Homo_sapiens_assembly19"
    assert model.panel_mode is True
    assert model.gc_bias is not None, "the gc_bias block was not parsed"


@pytest.mark.integration
def test_pon_zscore_calculation():
    """Test z-score calculation formula."""
    # Z-score = (value - mean) / std
    value = 120.0
    mean = 100.0
    std = 10.0

    zscore = (value - mean) / std

    assert zscore == 2.0
    assert abs(zscore) > 1.96  # Significant at 95% CI


@pytest.mark.integration
def test_pon_integration_load_model(tmp_path, sample_pon_parquet):
    """`load_pon_model` returns a populated model, not merely a non-None object.

    This asserted only `model is not None`, which the absence of an exception
    had already established. A loader that silently produced an empty model --
    every baseline `None`, `n_samples` 0 -- would have passed, and an empty
    model is exactly what a schema change or a renamed table produces.
    """
    from krewlyzer.core.pon_integration import load_pon_model

    model = load_pon_model(sample_pon_parquet)

    assert model is not None
    assert model.assay, "assay is unset; the metadata row was not parsed"
    assert model.n_samples > 0, f"n_samples is {model.n_samples}"
    populated = [
        name
        for name in ("gc_bias", "fsd_baseline", "wps_baseline", "mds_baseline")
        if getattr(model, name, None) is not None
    ]
    assert populated, (
        "every baseline is None -- the model loaded but carries no data, which "
        "is what a renamed table or a schema change looks like"
    )


@pytest.mark.integration
def test_pon_missing_file():
    """Test PON loading handles missing file gracefully."""
    from krewlyzer.core.pon_integration import load_pon_model

    result = load_pon_model(Path("/nonexistent/path/pon.parquet"))

    assert result is None


@pytest.mark.unit
def test_pon_zscore_edge_cases():
    """Test z-score edge cases."""
    # Zero std should be handled
    value = 100.0
    mean = 100.0
    std = 0.0

    if std > 0:
        zscore = (value - mean) / std
    else:
        zscore = 0.0  # Default for zero std

    assert zscore == 0.0
