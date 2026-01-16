# Panel of Normals (PON)

The Panel of Normals (PON) is a unified model built from healthy plasma samples that enables:

1. **GC bias correction** - Per-fragment correction for GC content bias
2. **Z-score normalization** - Detect deviations from healthy baseline for all features
3. **Panel mode support** - Dual on/off-target baselines for capture panels

## Quick Start

```bash
# Build PON from healthy samples
krewlyzer build-pon samples.txt --assay msk-access-v2 -r hg19.fa -o pon.parquet

# Use PON for sample processing
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ -P pon.parquet
```

## PON Components

| Component | Description | Used By |
|-----------|-------------|---------|
| **GC Bias Model** | Expected coverage by GC content per fragment type | FSC, FSR, WPS |
| **FSD Baseline** | Size distribution per chromosome arm | FSD |
| **WPS Baseline** | WPS mean/std per transcript region | WPS |
| **OCF Baseline** | Open chromatin scores per region | OCF |
| **MDS Baseline** | k-mer frequencies and motif diversity | Motif |

## Panel Mode

For capture panels (like MSK-ACCESS), use `--target-regions` when building the PON:

```bash
krewlyzer build-pon samples.txt -a msk-access-v2 -r hg19.fa -T targets.bed -o pon.parquet
```

This enables:

- **GC model trained on off-target only** - Avoids capture bias
- **Separate on/off-target baselines** - For features that differ in captured regions
- **Panel mode detection** - Sample processing auto-detects matching PON mode

## Building a PON

See [build-pon CLI](../features/build-pon.md) for detailed options.

**Requirements:**
- Minimum 10 healthy samples recommended
- Same assay/panel as samples to be processed
- Same reference genome

## Using PON in Processing

When `--pon-model` is provided to `run-all`:

1. PON is loaded once and passed to all processors
2. Each feature computes z-scores against healthy baseline
3. Output includes both raw values and PON-normalized columns

## Output Columns

With PON, additional columns are added to outputs:

| Feature | PON Column | Description |
|---------|------------|-------------|
| FSC | `*_log2` | Log2 ratio vs PON expected |
| FSD | `ratio_log2` | Size distribution log ratio |
| WPS | `wps_zscore` | Z-score vs region baseline |
| OCF | `ocf_zscore` | Z-score vs OCF baseline |

## API Reference

```python
from krewlyzer.pon.model import PonModel

# Load existing PON
pon = PonModel.load("path/to/pon.parquet")

# Access components
gc_expected = pon.get_mean("short")  # Expected at median GC
variance = pon.get_variance("short")  # For reliability scoring

# Check panel mode
if pon.panel_mode:
    print(f"Built with: {pon.target_regions_file}")
```
