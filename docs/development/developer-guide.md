# Developer Guide

This guide covers the Krewlyzer codebase architecture for contributors.

## Repository Structure

```
krewlyzer/
├── src/krewlyzer/          # Python package
│   ├── cli.py              # Typer CLI entry point
│   ├── wrapper.py          # run-all orchestration
│   ├── assets.py           # AssetManager for bundled data
│   ├── extract.py          # BAM → BED extraction
│   ├── fsc.py              # Fragment size coverage
│   ├── fsd.py              # Fragment size distribution
│   ├── fsr.py              # Fragment size ratio
│   ├── wps.py              # Windowed protection score
│   ├── ocf.py              # Orientation-aware fragmentation
│   ├── motif.py            # End motif analysis
│   ├── mfsd.py             # Mutant fragment size distribution
│   ├── region_entropy.py   # TFBS/ATAC region entropy
│   ├── region_mds.py       # Per-gene MDS
│   ├── uxm.py              # Fragment-level methylation
│   ├── build_gc_reference.py # GC reference generation
│   ├── core/               # Shared processors
│   │   ├── asset_resolution.py  # Target/PON resolution logic
│   │   ├── asset_validation.py  # Asset validation checks
│   │   ├── bam_utils.py         # BAM utilities
│   │   ├── feature_serializer.py # JSON output
│   │   ├── fsc_processor.py     # FSC post-processing
│   │   ├── fsd_processor.py     # FSD post-processing
│   │   ├── fsr_processor.py     # FSR post-processing
│   │   ├── gc_assets.py         # GC resolution helper
│   │   ├── gene_bed.py          # Gene BED parsing
│   │   ├── logging.py           # Startup banner and logging
│   │   ├── motif_processor.py   # Motif post-processing
│   │   ├── ocf_processor.py     # OCF post-processing
│   │   ├── pon_integration.py   # PON post-processing
│   │   ├── region_entropy_processor.py # TFBS/ATAC processor
│   │   ├── sample_processor.py  # Per-sample orchestration
│   │   ├── unified_processor.py # Unified pipeline Python layer
│   │   ├── wps_processor.py     # WPS post-processing
│   │   └── utils.py, resource_utils.py
│   ├── pon/                # PON model code
│   │   ├── model.py        # PonModel dataclass
│   │   └── build.py        # PON building logic
│   └── data/               # Bundled assets (Git LFS)
├── rust/                   # Rust backend (19 modules)
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs          # PyO3 module exports
│       ├── pipeline.rs     # Unified pipeline entry
│       ├── engine.rs       # Core engine utilities
│       ├── bed.rs          # BGZF/gzip BED reader
│       ├── filters.rs      # Fragment filtering logic
│       └── (feature modules: see table below)
├── tests/                  # Test suite (248 tests)
├── docs/                   # MkDocs documentation
└── nextflow/               # Nextflow pipeline
    ├── main.nf
    ├── nextflow.config
    ├── modules/local/       # Per-tool NF modules
    └── subworkflows/local/  # Subworkflows
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

| Module | Purpose |
|--------|---------||
| `lib.rs` | PyO3 module exports and thread config |
| `pipeline.rs` | Unified pipeline coordination |
| `engine.rs` | Core engine utilities |
| `bed.rs` | BGZF/gzip BED reader |
| `extract_motif.rs` | BAM parsing, fragment + motif extraction |
| `motif_utils.rs` | Shared 4-mer encoding, MDS, GC utils |
| `fsc.rs` | Fragment size coverage + gene aggregation |
| `fsd.rs` | Per-arm size distribution + PON log-ratio |
| `wps.rs` | Dual-stream WPS, FFT, smoothing |
| `ocf.rs` | Orientation-aware fragmentation + PON z-score |
| `mfsd.rs` | Mutant fragment size distribution |
| `region_entropy.rs` | TFBS/ATAC entropy + PON z-score |
| `region_mds.rs` | Per-gene MDS at exon boundaries |
| `uxm.rs` | Fragment-level methylation (UXM) |
| `gc_correction.rs` | LOESS GC bias correction |
| `gc_reference.rs` | Pre-computed GC reference generation |
| `pon_model.rs` | PON model loading and hybrid correction |
| `pon_builder.rs` | PON model construction |
| `filters.rs` | Fragment filtering logic |


### Unified Pipeline

All feature computation goes through `run_unified_pipeline`:

```rust
pub fn run_unified_pipeline(
    _py: Python,
    bed_path: PathBuf,
    // GC Correction
    gc_ref_path: Option<PathBuf>,
    valid_regions_path: Option<PathBuf>,
    correction_out_path: Option<PathBuf>,
    correction_input_path: Option<PathBuf>,
    // FSC
    fsc_bins: Option<PathBuf>, fsc_output: Option<PathBuf>,
    // WPS Foreground (TSS/CTCF anchors)
    wps_regions: Option<PathBuf>, wps_output: Option<PathBuf>,
    // WPS Background (Alu stacking)
    wps_background_regions: Option<PathBuf>, wps_background_output: Option<PathBuf>,
    wps_empty: bool,
    // FSD
    fsd_arms: Option<PathBuf>, fsd_output: Option<PathBuf>,
    // OCF
    ocf_regions: Option<PathBuf>, ocf_output: Option<PathBuf>,
    // Target regions for on/off-target split (panel mode)
    target_regions_path: Option<PathBuf>,
    bait_padding: u64,
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
serializer = FeatureSerializer(sample_id, version="0.5.3")
serializer.add_fsc(fsc_df)
serializer.add_fsc_e1(fsc_e1_df)
serializer.add_wps(wps_df)
serializer.add_motif(edm_df, bpm_df, mds, mds_z)
serializer.add_ocf(ocf_df)
serializer.add_mfsd(mfsd_df)
serializer.add_uxm(uxm_df)
serializer.add_qc("total_fragments", 1234567)
serializer.save(output_dir)

# Or load from existing outputs
serializer = FeatureSerializer.from_outputs(
    sample_id=sample_id,
    output_dir=output_dir,
    version="0.5.3"
)
```

---

## Testing

Krewlyzer has **248 tests** across unit, integration, and e2e categories.

**→ [Testing Guide](testing-guide.md)** for complete documentation including:
- Feature → test file mapping
- Fixture reference
- Test writing examples

### Quick Commands

```bash
# All tests
pytest tests/

# Fast unit tests only
pytest tests/unit/

# Specific feature
pytest tests/unit/test_fsc.py

# With coverage
pytest tests/ --cov=krewlyzer --cov-report=html
```

!!! tip
    **Modifying a feature?** Check the [Feature → Test Map](testing-guide.md#feature--test-map) to find which test file to update.

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

See [CONTRIBUTING.md](contributing.md) for full guidelines.
