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

## Auto-PON Loading

When you specify an assay with `-A`, krewlyzer automatically loads the bundled PON:

```bash
# Auto-loads bundled PON for xs2 assay
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ -A xs2 -G hg19
```

This is equivalent to explicitly passing `-P` with the bundled PON path.

## Skipping Z-Score Normalization (`--skip-pon`)

For **ML training workflows** where PON samples are used as true negatives, use `--skip-pon` to output raw features without z-score normalization:

```bash
# Process PON samples as ML negatives (no z-scores)
krewlyzer run-all -i pon_sample.bam -r hg19.fa -o out/ -A xs2 --skip-pon
```

> [!WARNING]
> `-P` and `--skip-pon` are **mutually exclusive**. If you specify an explicit PON model, you want z-scores applied. Use `--skip-pon` only with `-A` (assay) for the ML negatives workflow.

The `--skip-pon` flag:
- Works with `-A/--assay` (auto-loads bundled PON but skips z-scores)
- Available on all tools: `run-all`, `fsc`, `fsd`, `fsr`, `wps`, `ocf`, `region-entropy`, `motif`
- Logs which tools are skipping normalization

## PON Components

| Component | Description | Used By |
|-----------|-------------|---------|
| **GC Bias Model** | Expected coverage by GC content per fragment type | FSC, FSR, WPS |
| **FSD Baseline** | Size distribution per chromosome arm | FSD |
| **WPS Baseline** | WPS mean/std per transcript region | WPS |
| **OCF Baseline** | Open chromatin scores per region | OCF |
| **MDS Baseline** | k-mer frequencies and motif diversity | Motif |
| **TFBS Baseline** | Per-TF entropy mean/std | Region Entropy |
| **ATAC Baseline** | Per-cancer-type entropy mean/std | Region Entropy |
| **Region MDS Baseline** | Per-gene MDS mean/std for E1 | Region MDS |
| **FSC Gene Baseline** | Per-gene normalized depth mean/std | FSC Gene |
| **FSC Region Baseline** | Per-exon normalized depth mean/std | FSC Region |

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

| Feature | PON Column(s) | Description |
|---------|--------------|-------------|
| FSC | `*_log2` | Log2 ratio vs PON expected |
| FSC Gene | `depth_zscore` | Gene-level depth z-score |
| FSC Region | `depth_zscore` | Exon-level depth z-score |
| FSD | `ratio_log2` | Size distribution log ratio |
| WPS | `wps_zscore` | Z-score vs region baseline |
| OCF | `ocf_zscore` | Z-score vs OCF baseline |
| Motif | `mds_z` | Z-score for MDS |
| TFBS | `entropy_z` | Z-score per TF |
| ATAC | `entropy_z` | Z-score per cancer type |
| Region MDS | `mds_z`, `mds_e1_z` | Gene-level and E1 z-scores |

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

Per-region nucleosome positioning metrics.

**Schema v1.0 (Scalar):**
- `wps_long_mean/std`: Single nucleosomal WPS value per region
- `wps_short_mean/std`: Single TF footprint value per region

**Schema v2.0 (Vector):**
- `wps_nuc_mean/std`: 200-element vector (nucleosomal footprint)
- `wps_tf_mean/std`: 200-element vector (TF footprint)

> [!TIP]
> v2.0 enables position-specific z-scores and **Shape Correlation Score** for cancer detection.

#### Shape Score Interpretation

| Score | Interpretation |
|-------|---------------|
| 0.9-1.0 | Healthy nucleosome positioning |
| 0.5-0.9 | Mild chromatin disorganization |
| <0.5 | Significant disruption (cancer signal) |

See [WPS Features](../features/wps.md) for output column details.

### OCF Baseline (`ocf_baseline`)

Per-region open chromatin footprint:
- `ocf_mean/std`: OCF score baseline
- `sync_mean/std`: Synchronization score baseline

### MDS Baseline (`mds_baseline`)

Motif diversity expectations:
- `kmer_expected`: 256 4-mer frequencies from healthy samples
- `kmer_std`: Variability per k-mer
- `mds_mean/std`: Expected Motif Diversity Score

### TFBS Baseline (`tfbs_baseline`)

Per-TF size entropy:
- `label_stats`: Mean/std entropy per TF (808 transcription factors)
- Enables z-score per TF for detailed regulatory analysis

### ATAC Baseline (`atac_baseline`)

Per-cancer-type size entropy:
- `label_stats`: Mean/std entropy per cancer type (23 types)
- Enables tissue-of-origin scoring

### Region MDS Baseline (`region_mds`)

Per-gene MDS expectations:
- `gene_baseline`: Dict of gene → {mds_mean, mds_std, mds_e1_mean, mds_e1_std}
- Enables gene-level anomaly detection
- E1 (first exon) tracked separately for promoter-proximal sensitivity

### FSC Gene Baseline (`fsc_gene_baseline`)

Per-gene normalized depth baseline (panel mode only):
- `data`: Dict of gene → (mean_depth, std_depth, n_samples)
- Requires minimum 3 samples for reliable statistics
- Clinical use: z-score >> 0 = amplification, z-score << 0 = deletion

### FSC Region Baseline (`fsc_region_baseline`)

Per-exon/probe normalized depth baseline (panel mode only):
- `data`: Dict of region_id → (mean_depth, std_depth, n_samples)
- Region IDs formatted as "chrom:start-end"
- Covers all exons (no filtering by variance)
- Enables detection of focal copy number changes affecting single exons

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
| **TFBS** | `entropy_z` | Higher entropy (diverse sizes) | Lower entropy (restricted) |
| **ATAC** | `entropy_z` | Higher entropy | Lower entropy |
| **Region MDS** | `mds_z`, `mds_e1_z` | More diverse at gene | Restricted motifs (aberrant) |

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

