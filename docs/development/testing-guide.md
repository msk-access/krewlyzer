# Testing Guide

## Overview

Krewlyzer has **245 tests** covering all features via pytest.

| Category | Tests | Speed | Location |
|----------|:-----:|-------|----------|
| **Unit** | 155 | <1s | `tests/unit/` |
| **Integration** | 52 | 5-30s | `tests/integration/` |
| **CLI** | 10 | 2-5s | `tests/cli/` |
| **E2E** | 3 | 30-60s | `tests/e2e/` |
| **Asset Resolution** | 19 | <1s | `tests/test_asset_resolution.py` |

!!! note
    **Rust code is tested via Python.** The `test_rust_python_equivalence.py` suite verifies Rust output matches Python implementations.

---

<a id="feature--test-map"></a>
## Feature → Test Map

**When modifying a feature, update the corresponding test file:**

| Feature | Test File(s) | Tests |
|---------|--------------|:-----:|
| **FSC** | `unit/test_fsc.py` | 7 |
| **FSD** | `unit/test_fsd.py`, `integration/test_fsd_cli.py` | 8 |
| **FSR** | `integration/test_fsr_cli.py` | 3 |
| **WPS** | `unit/test_wps.py`, `integration/test_wps_cli.py` | 19 |
| **Motif** | `integration/test_motif.py` | 3 |
| **OCF** | `integration/test_ocf.py` | 1 |
| **mFSD** | `integration/test_mfsd.py` | 9 |
| **UXM** | `integration/test_uxm.py` | 1 |
| **Region Entropy** | `integration/test_region_entropy.py` | 10 |
| **Extract** | `integration/test_extract.py` | 4 |
| **PON Model** | `unit/test_pon_model.py` | 51 |
| **PON Validation** | `unit/test_pon_validation.py` | 14 |
| **PON Build** | `integration/test_pon.py`, `test_pon_dual_gc.py` | 8 |
| **Gene BED** | `unit/test_gene_bed.py` | 25 |
| **Asset Manager** | `unit/test_asset_manager.py` | 10 |
| **External Data Dir** | `unit/test_external_data_dir.py` | 6 |
| **Asset Resolution** | `test_asset_resolution.py` | 19 |
| **BGZF Reader** | `unit/test_bgzf_reader.py` | 5 |
| **Normalization** | `unit/test_normalization.py` | 6 |
| **run-all flags** | `unit/test_run_all_flags.py` | 6 |
| **Rust/Python** | `unit/test_rust_python_equivalence.py` | 11 |
| **Real Data** | `integration/test_real_data.py` | 5 |

---

## Running Tests

```bash
# All tests
pytest tests/

# By category
pytest tests/unit/           # Fast unit tests
pytest tests/integration/    # Tool integration
pytest tests/e2e/            # Full pipeline

# Specific feature
pytest tests/unit/test_fsc.py
pytest tests/unit/test_wps.py -v

# Stop on first failure
pytest -x

# With coverage
pytest tests/ --cov=krewlyzer --cov-report=html
```

---

## Test Markers

```bash
pytest -m unit           # Unit tests only
pytest -m integration    # Integration tests
pytest -m "not slow"     # Skip slow tests
```

---

## Data Availability

!!! important
    The **entire `src/krewlyzer/data/` folder is EXCLUDED from PyPI wheels** to keep size <100MB.
    Tests that require bundled data files are **automatically skipped** in CI/PyPI installs.

### Data by Install Method

| Install Method | Data Files | Test Coverage |
|----------------|:----------:|:-------------:|
| `pip install krewlyzer` (PyPI) | ❌ None | ~85% (skips asset tests) |
| `pip install -e .` (git clone) | ✅ All | 100% |
| Docker image | ✅ All | 100% |

### How It Works

Tests that verify bundled assets use the `@requires_data` decorator from `conftest.py`:

```python
from conftest import requires_data

@requires_data
class TestAssetManager:
    def test_gene_bed_exists(self):
        # Skipped if data not available
        ...
```

### Running Full Tests Locally

```bash
# Clone the repository (includes data/)
git clone https://github.com/msk-access/krewlyzer.git
cd krewlyzer

# Development install (uses source data directly)
pip install -e ".[test]"

# Run all tests - data-dependent tests will pass
pytest tests/ -v
```

### External Data Directory

For PyPI installs, you can provide data via environment variable:

```bash
export KREWLYZER_DATA_DIR=/path/to/krewlyzer/src/krewlyzer/data
pytest tests/ -v
```

---

## Fixtures

Key fixtures from `tests/conftest.py`:

| Fixture | Description |
|---------|-------------|
| `temp_bam` | Minimal BAM with proper pair |
| `temp_bedgz` | Minimal BED.gz (3 fragments) |
| `temp_reference` | FASTA reference (12kb chr1) |
| `temp_bins` | FSC bins file |
| `temp_arms` | FSD arms file |
| `temp_transcripts` | WPS transcripts file |
| `temp_ocr` | OCF regions file |
| `temp_vcf` | mFSD variants file |
| `full_test_data` | Bundle for run-all |
| `real_bam` | Production data (3144 reads) |
| `real_pon` | MSK-ACCESS v1 PON |

---

## Adding Tests

### Naming Convention
- **File**: `test_<feature>.py`
- **Function**: `test_<action>_<expected>()`

### Example Unit Test
```python
@pytest.mark.unit
def test_fsc_counts_fragments_by_size(temp_bedgz, temp_bins):
    """FSC should categorize fragments into size bins."""
    result = fsc.count_fragments(temp_bedgz, temp_bins)
    
    assert "core_short" in result.columns
    assert result["total"].sum() > 0
```

### Example Integration Test
```python
@pytest.mark.integration
def test_fsc_cli_produces_output(temp_bedgz, temp_bins, tmp_path):
    """FSC CLI should write TSV output."""
    from typer.testing import CliRunner
    from krewlyzer.cli import app
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "fsc", "-i", str(temp_bedgz), 
        "-b", str(temp_bins), "-o", str(tmp_path)
    ])
    
    assert result.exit_code == 0
    assert (tmp_path / "test.FSC.tsv").exists()
```

### Example Data-Dependent Test
```python
from conftest import requires_data

@requires_data
def test_bundled_gene_bed_loads(manager):
    """Test bundled gene BED (only runs with source data)."""
    path = manager.get_gene_bed("xs1")
    assert path.exists()
```

---

## See Also

- [Developer Guide](developer-guide.md) - Setup and architecture
- [Contributing](contributing.md) - PR process
