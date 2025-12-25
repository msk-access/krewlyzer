# Krewlyzer Test Suite

## Structure

```
tests/
├── conftest.py                 # Shared fixtures (15+ fixtures)
├── README.md                   # This file
├── cli/                        # CLI tests (10 tests)
│   └── test_cli.py             # --help output validation
├── unit/                       # Unit tests (20 tests)
│   ├── test_fsc.py             # FSC fragment counting (3)
│   ├── test_fsd.py             # FSD ratio calculations (5)
│   ├── test_wps.py             # WPS score formulas (6)
│   └── test_normalization.py   # Depth/GC normalization (6)
├── integration/                # Integration tests (14 tests)
│   ├── test_extract.py         # Extract blacklist/filtering (1)
│   ├── test_motif.py           # Motif extraction (4)
│   ├── test_mfsd.py            # mFSD variant classification (8)
│   ├── test_ocf.py             # OCF processing (1)
│   └── test_uxm.py             # UXM methylation (1)
├── e2e/                        # End-to-end tests (1 test)
│   └── test_run_all.py         # Full pipeline
└── data/                       # Test data
    ├── fixtures/               # Pre-generated fixtures
    └── create_test_data.py     # Data generator script
```

## Running Tests

```bash
# All tests
pytest tests/

# By category
pytest tests/unit/           # Fast unit tests (~0.5s)
pytest tests/integration/    # Integration tests
pytest tests/e2e/            # End-to-end tests

# By marker
pytest -m unit               # Unit tests only
pytest -m integration        # Integration tests
pytest -m rust               # Tests requiring Rust backend

# With coverage
pytest tests/ --cov=krewlyzer
```

## Test Markers

| Marker | Description |
|--------|-------------|
| `unit` | Fast unit tests, no I/O |
| `integration` | Tool integration tests |
| `e2e` | End-to-end workflow tests |
| `slow` | Tests taking >10s |
| `rust` | Tests requiring Rust backend |

## Shared Fixtures (conftest.py)

| Fixture | Description |
|---------|-------------|
| `temp_bam` | Minimal BAM with proper pair |
| `temp_bedgz` | Minimal BED.gz with fragments |
| `temp_reference` | FASTA reference (2kb) |
| `temp_bins` | FSC bins file |
| `temp_arms` | FSD arms file |
| `temp_transcripts` | WPS transcripts file |
| `temp_ocr` | OCF regions file |
| `temp_markers` | UXM markers file |
| `temp_vcf` | mFSD variants file |
| `full_test_data` | Complete bundle for run-all |

## Current Status

| Category | Passed | Skipped | Total |
|----------|--------|---------|-------|
| CLI | 10 | 0 | 10 |
| Unit | 20 | 0 | 20 |
| Integration | 32 | 4 | 36 |
| E2E | 1 | 0 | 1 |
| **Total** | **63** | **4** | **67** |

**Coverage: 46%**

✅ All executing tests pass

## Real Data Fixtures

Test fixtures in `tests/data/fixtures/`:
- `test_sample.bam` - 3,144 reads from chr1:1-2000000
- `test.pon.parquet` - MSK-ACCESS v1 PON model
- `test_genome.fa` - 2Mb chr1 reference
- `test_bins.bed` - 20 x 100kb bins
- `test_arms.bed` - Chr1 p/q arms
- `test_transcripts.tsv` - 3 test genes
