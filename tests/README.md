# Krewlyzer Test Suite

## Overview

**239 tests** covering all features:

| Category | Tests | Speed | Location |
|----------|:-----:|-------|----------|
| **Unit** | 155 | <1s | `tests/unit/` |
| **Integration** | 52 | 5-30s | `tests/integration/` |
| **CLI** | 10 | 2-5s | `tests/cli/` |
| **E2E** | 3 | 30-60s | `tests/e2e/` |
| **Asset Resolution** | 19 | <1s | `tests/test_asset_resolution.py` |

## Structure

```
tests/
├── conftest.py             # Shared fixtures (15+ fixtures)
├── README.md               # This file
├── test_asset_resolution.py # Asset resolution (19 tests)
├── cli/                    # CLI help tests (10)
├── unit/                   # Unit tests (155)
│   ├── test_fsc.py         # FSC fragment counting
│   ├── test_fsd.py         # FSD ratio calculations
│   ├── test_wps.py         # WPS score formulas
│   ├── test_pon_model.py   # PON model loading (51)
│   ├── test_gene_bed.py    # Gene BED parsing (25)
│   └── ...
├── integration/            # Integration tests (52)
│   ├── test_extract.py     # Extract pipeline
│   ├── test_mfsd.py        # mFSD variant classification
│   ├── test_region_entropy.py # TFBS/ATAC entropy
│   └── ...
├── e2e/                    # End-to-end tests (3)
│   └── test_run_all.py     # Full pipeline
└── data/                   # Test data
    └── fixtures/           # Pre-generated fixtures
```

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

# Stop on first failure
pytest -x

# With coverage
pytest tests/ --cov=krewlyzer --cov-report=html
```

## Test Markers

| Marker | Description |
|--------|-------------|
| `unit` | Fast unit tests, no I/O |
| `integration` | Tool integration tests |
| `e2e` | End-to-end workflow tests |
| `slow` | Tests taking >10s |

## Fixtures

Key fixtures in `conftest.py`:

| Fixture | Description |
|---------|-------------|
| `temp_bam` | Minimal BAM with proper pair |
| `temp_bedgz` | Minimal BED.gz (3 fragments) |
| `temp_reference` | FASTA reference (12kb) |
| `temp_bins`, `temp_arms` | FSC/FSD resources |
| `full_test_data` | Complete bundle for run-all |
| `real_bam`, `real_pon` | Production data fixtures |

## Documentation

See [Testing Guide](../docs/development/testing-guide.md) for:
- Feature → test file mapping
- Test writing examples
- Fixture reference
