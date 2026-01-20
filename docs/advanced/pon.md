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

---

## PON Baselines in Detail

### GC Bias Model (`gc_bias`)

Stores expected fragment coverage per GC content (0-100%) for each fragment type:

| Fragment Type | Size Range | Purpose |
|---------------|------------|---------|
| `short` | 65-149bp | Short fragment correction |
| `intermediate` | 150-259bp | Mono-nucleosomal |
| `long` | 260-400bp | Di-nucleosomal |
| `wps_long` | 120-180bp | WPS nucleosomal |
| `wps_short` | 35-80bp | WPS TF footprint |

### FSD Baseline (`fsd_baseline`)

Size distribution per chromosome arm (46 arms):
- `expected`: Mean proportion at each size bin
- `std`: Standard deviation across PON samples

### WPS Baseline (`wps_baseline`)

Per-region nucleosome metrics:
- `wps_long_mean/std`: Nucleosomal WPS baseline
- `wps_short_mean/std`: TF footprint baseline
- Schema v2.0: Full 200-bin vectors for position-specific z-scores

### OCF Baseline (`ocf_baseline`)

Per-region open chromatin footprint:
- `ocf_mean/std`: OCF score baseline
- `sync_mean/std`: Synchronization score baseline

### MDS Baseline (`mds_baseline`)

Motif diversity expectations:
- `kmer_expected`: 256 4-mer frequencies from healthy samples
- `kmer_std`: Variability per k-mer
- `mds_mean/std`: Expected Motif Diversity Score

---

## Interpreting Z-Scores

Z-scores measure how many standard deviations a sample differs from the healthy PON baseline:

$$
z = \frac{x_{\text{sample}} - \mu_{\text{PON}}}{\sigma_{\text{PON}}}
$$

### Clinical Interpretation

| Z-Score Range | Interpretation | Action |
|---------------|----------------|--------|
| **-2 to +2** | Normal range | Within healthy variation |
| **|z| = 2-3** | Mild deviation | Monitor, may be noise |
| **|z| > 3** | Significant | Investigate for ctDNA |
| **|z| > 5** | Extreme | High tumor burden likely |

### Per-Feature Z-Score Meaning

| Feature | Z-Score Column | Positive Z Means | Negative Z Means |
|---------|----------------|------------------|------------------|
| **FSC** | `z_core_short` | More short fragments | Fewer short fragments |
| **FSD** | - | Shifted size distribution | - |
| **WPS** | `wps_nuc_z` | Stronger nucleosome signal | Disrupted nucleosomes |
| **OCF** | `ocf_z` | More open chromatin | Less accessible |
| **MDS** | `mds_z` | More diverse motifs | Less diverse |

### ML Feature Usage

```python
# Extract z-score features for classification
features = {
    "fsc_short_z": sample_fsc["z_core_short"].mean(),
    "wps_nuc_z": sample_wps["wps_nuc_z"].mean(),
    "mds_z": sample_motif["mds_z"],
}

# Higher |z| = more likely to be tumor
combined_signal = sum(abs(z) for z in features.values())
```

> [!TIP]
> **Combine z-scores across features** - Single extreme values may be noise, but consistent deviations across FSC, WPS, and MDS are highly indicative of ctDNA.

