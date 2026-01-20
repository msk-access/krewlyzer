# Developer Guide

This guide covers the Krewlyzer codebase architecture for contributors.

## Repository Structure

```
krewlyzer/
├── src/krewlyzer/          # Python package
│   ├── cli.py              # Typer CLI entry point
│   ├── wrapper.py          # run-all orchestration (750 lines)
│   ├── assets.py           # AssetManager for bundled data
│   ├── extract.py          # BAM → BED extraction
│   ├── fsc.py, fsd.py, ... # Standalone tools
│   ├── core/               # Shared processors
│   │   ├── gc_assets.py    # GC resolution helper
│   │   ├── fsc_processor.py
│   │   ├── wps_processor.py
│   │   └── feature_serializer.py  # JSON output
│   ├── pon/                # PON model code
│   │   ├── model.py        # PonModel dataclass
│   │   └── build.py        # PON building logic
│   └── data/               # Bundled assets
├── rust/                   # Rust backend
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs          # PyO3 module exports
│       ├── pipeline.rs     # Unified pipeline entry
│       ├── fsc.rs, wps.rs  # Feature modules
│       └── gc_correction.rs
├── tests/                  # Test suite
├── docs/                   # MkDocs documentation
└── nextflow/               # Nextflow pipeline
```

---

## Rust Backend Architecture

### Module Structure

```
lib.rs
├── extract_motif     # BAM extraction module
├── gc                # GC correction module  
├── run_unified_pipeline  # Main pipeline function
└── configure_threads # Thread pool setup
```

### Key Rust Modules

| Module | Lines | Purpose |
|--------|------:|---------|
| `wps.rs` | 2000+ | Dual-stream WPS, FFT, smoothing |
| `gc_correction.rs` | 590 | LOESS GC bias correction |
| `gc_reference.rs` | 658 | Asset generation (once per genome) |
| `fsc.rs` | 500+ | 5-bin fragment counting |
| `fsd.rs` | 600+ | Per-arm size distribution |
| `pipeline.rs` | 400+ | Unified pipeline coordination |
| `extract_motif.rs` | 1200+ | BAM reading, motif extraction |

### Unified Pipeline

All feature computation goes through `run_unified_pipeline`:

```rust
pub fn run_unified_pipeline(
    bedgz_path: &str,
    gc_ref_path: Option<&str>,
    valid_regions_path: Option<&str>,
    gc_factors_out: Option<&str>,
    gc_factors_in: Option<&str>,
    fsc_bins_path: Option<&str>,
    fsc_out_path: Option<&str>,
    wps_anchors_path: Option<&str>,
    wps_out_path: Option<&str>,
    wps_bg_path: Option<&str>,
    wps_bg_out_path: Option<&str>,
    include_empty: bool,
    fsd_arms_path: Option<&str>,
    fsd_out_path: Option<&str>,
    ocf_regions_path: Option<&str>,
    ocf_out_path: Option<&str>,
    target_regions_path: Option<&str>,
    bait_padding: u32,
    silent: bool,
) -> PyResult<()>
```

### Adding a New Feature

1. **Create Rust module** (`rust/src/new_feature.rs`)
2. **Add to lib.rs** module exports
3. **Create Python wrapper** (`src/krewlyzer/new_feature.py`)
4. **Add to CLI** in `cli.py`
5. **Integrate with run_unified_pipeline** if applicable
6. **Add to FeatureSerializer** for JSON output
7. **Write tests** in `tests/unit/` and `tests/integration/`

---

## Python Architecture

### Tool Pattern

All standalone tools follow this pattern:

```python
def tool_name(
    input: Path = typer.Option(...),
    output: Path = typer.Option(...),
    # ... options
):
    # 1. Initialize AssetManager
    assets = AssetManager(genome)
    
    # 2. Resolve GC assets (use helper)
    from .core.gc_assets import resolve_gc_assets
    gc = resolve_gc_assets(assets, output, sample_name, input, gc_correct, genome)
    
    # 3. Call Rust backend
    _core.run_unified_pipeline(...)
    
    # 4. Post-process (Python)
    from .core.tool_processor import process_tool
    process_tool(output_file, ...)
```

### AssetManager

Centralizes access to bundled data:

```python
assets = AssetManager("hg19")

# Properties
assets.gc_reference      # GC ref parquet
assets.valid_regions     # Valid regions BED
assets.bins_100kb        # Default bins
assets.wps_anchors       # WPS anchors BED

# Assay-aware methods
assets.get_gene_bed("xs2")     # Panel gene BED
assets.get_wps_anchors("xs2")  # Panel WPS anchors
assets.get_pon("xs2")          # Bundled PON
```

### FeatureSerializer

Collects all features into unified JSON:

```python
serializer = FeatureSerializer(sample_id)
serializer.add_fsc(fsc_df)
serializer.add_wps(wps_df)
serializer.add_motif(edm_dict, bpm_dict, mds)
# ...
serializer.save(output_file)

# Or load from existing outputs
serializer = FeatureSerializer.from_outputs(output_dir, sample_id)
```

---

## Testing Locally

### Setup

```bash
# Clone and install
git clone https://github.com/msk-access/krewlyzer.git
cd krewlyzer

# Create venv
uv venv .venv
source .venv/bin/activate

# Install with test dependencies
uv pip install -e ".[test]"
```

### Running Tests

```bash
# All tests
pytest tests/

# By category
pytest tests/unit/           # Fast (~0.5s)
pytest tests/integration/    # Full tool tests
pytest tests/e2e/            # End-to-end

# With coverage
pytest tests/ --cov=krewlyzer --cov-report=html
```

### Test Markers

```bash
pytest -m unit           # Unit tests only
pytest -m integration    # Integration tests
pytest -m rust           # Tests requiring Rust backend
pytest -m "not slow"     # Skip slow tests
```

### Test Fixtures

Key fixtures in `tests/conftest.py`:

| Fixture | Description |
|---------|-------------|
| `temp_bam` | Minimal BAM with reads |
| `temp_bedgz` | Pre-extracted BED.gz |
| `temp_reference` | 2kb FASTA reference |
| `full_test_data` | Complete bundle for run-all |

---

## Rust Development

### Building

```bash
cd rust

# Debug build (fast compile, slow run)
maturin develop

# Release build (slow compile, fast run)
maturin develop --release
```

### Testing Rust

```bash
# Run Rust tests (currently minimal)
cargo test

# Check before committing
cargo clippy
cargo fmt --check
```

### Debugging

```bash
# Enable verbose logging
RUST_LOG=debug krewlyzer run-all -i sample.bam ...

# Or in Python
python -c "
import logging
logging.basicConfig(level=logging.DEBUG)
from krewlyzer import _core
# ... your test code
"
```

---

## Code Style

### Python

- Use `typer` for CLIs
- Use `rich` for logging/progress
- Type hints for all public functions
- Docstrings (Google style)

### Rust

- Run `cargo fmt` before committing
- Address `cargo clippy` warnings
- Use `anyhow::Result` for error handling
- Log with `log::info!`, `log::debug!`

---

## Contributing Checklist

- [ ] Code follows existing patterns
- [ ] Added/updated tests
- [ ] Updated documentation
- [ ] Ran `pytest tests/`
- [ ] Ran `cargo fmt && cargo clippy`
- [ ] Updated CHANGELOG.md

See [CONTRIBUTING.md](../contributing.md) for full guidelines.
