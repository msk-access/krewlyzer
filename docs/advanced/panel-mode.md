# Panel Mode

Panel mode enables accurate cfDNA analysis for capture-based sequencing panels like MSK-ACCESS.

## Overview

When using targeted capture panels, two key issues affect cfDNA analysis:

1. **GC Bias**: Capture probes introduce additional GC bias on top of sequencing bias
2. **Coverage Splitting**: On-target fragments behave differently than off-target fragments

Panel mode addresses both by:
- Training the GC model on **off-target fragments only** (unbiased by capture)
- Computing **dual baselines** for on-target and off-target regions separately

## Enabling Panel Mode

### Building GC Reference Assets (One-time)

For panel mode, generate panel-specific GC reference assets:

```bash
krewlyzer build-gc-reference hg19.fa -o data/gc/ \
    --target-regions msk_access_baits.bed
```

This generates both standard and on-target GC reference files.

### At PON Build Time

```bash
krewlyzer build-pon samples.txt \
    --assay msk-access-v2 \
    --reference hg19.fa \
    --target-regions msk_access_baits.bed \
    --output msk-access.pon.parquet
```

### At Sample Processing Time

The `--assay` flag enables panel-specific optimizations:

```bash
# MSK-ACCESS v2 with all panel features
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ \
    --assay xs2 \
    --target-regions msk_access_baits.bed \
    --pon-model msk-access.pon.parquet
```

### What `--assay` Enables

| Feature | Without --assay | With --assay |
|---------|-----------------|--------------|
| **Gene FSC** | Window-based only | + Gene-level aggregation (`FSC.gene.tsv`) |
| **WPS Anchors** | Genome-wide (~15k) | Panel-specific (~2k + genome-wide) |
| **WPS Output** | Single `WPS.parquet` | Dual: `WPS.parquet` + `WPS.panel.parquet` |
| **JSON Output** | Standard features | + `fsc_gene`, `wps_panel` |

### Dual WPS Output

With `--assay`, Krewlyzer generates **two** WPS files:

| File | Anchors | Use Case |
|------|---------|----------|
| `{sample}.WPS.parquet` | Genome-wide TSS+CTCF | Cancer detection signature |
| `{sample}.WPS.panel.parquet` | Panel gene anchors | Targeted gene profiling |

This dual output provides both broad cancer signals and focused gene-level analysis.

### Minimal Panel Mode (No PON)

For quick analysis without a custom PON:

```bash
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ \
    --assay xs2 \
    --target-regions targets.bed
```

### Auto-Loading Assets with `--assay`

When you specify `--assay`, krewlyzer automatically loads bundled assets:

```bash
# Auto-loads PON and target regions for xs2
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ --assay xs2
```

| Assay | PON Model | Target Regions |
|-------|-----------|----------------|
| `xs1` | `xs1.all_unique.pon.parquet` | `xs1.targets.bed.gz` |
| `xs2` | `xs2.all_unique.pon.parquet` | `xs2.targets.bed.gz` |
| `wgs` | None | None (WGS mode) |

This auto-loading applies to **all tools** including:
- `extract`, `run-all`, `fsc`, `fsd`, `fsr`, `wps`, `ocf`, `region-entropy`, `motif`, `region-mds`, `build-pon`

### Forcing WGS Mode with `--skip-target-regions`

To force WGS-like behavior even when using a panel assay (e.g., for comparison):

```bash
# Use xs2 PON but disable panel mode (process as WGS)
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ \
    --assay xs2 \
    --skip-target-regions
```

This is useful when:
- Comparing panel samples to WGS baselines
- Running validation without on/off-target splitting
- Processing samples where target regions don't apply

> [!NOTE]
> `--skip-target-regions` only disables target region loading. The PON model is still loaded from `--assay` unless you also add `--skip-pon`.

### Flag Priority

Asset resolution follows this priority order:

1. **Explicit path** (`--target-regions path/to/file.bed`) - highest priority
2. **Skip flag** (`--skip-target-regions`) - forces WGS mode
3. **Bundled asset** (auto-loaded from `--assay`)
4. **None** - WGS mode (no targets)


## How It Works

### GC Correction

```
WGS Mode:     All fragments → GC model → Single correction curve
Panel Mode:   Off-target fragments → GC model → Unbiased correction curve
              On-target fragments → GC model → Capture-aware correction curve
```

The GC model is built from off-target fragments because:
- On-target fragments have probe-specific GC bias
- Off-target fragments represent natural cfDNA (similar to WGS)

### Dual Correction Factor Files

In panel mode, `krewlyzer extract` generates TWO correction factor files:

| File | Source | Used For |
|------|--------|----------|
| `{sample}.correction_factors.csv` | Off-target fragments | Primary biomarker analysis |
| `{sample}.correction_factors.ontarget.csv` | On-target fragments | Copy number, variant calling |

### Feature Splitting

In panel mode, each feature outputs two files:

| Feature | Off-Target File | On-Target File |
|---------|-----------------|----------------|
| FSC | `.FSC.tsv` | `.FSC.ontarget.tsv` |
| FSR | `.FSR.tsv` | `.FSR.ontarget.tsv` |
| FSD | `.FSD.tsv` | `.FSD.ontarget.tsv` |
| OCF | `.OCF.tsv` | `.OCF.ontarget.tsv` |
| TFBS | `.TFBS.tsv` | `.TFBS.ontarget.tsv` |
| ATAC | `.ATAC.tsv` | `.ATAC.ontarget.tsv` |

**Note on OCF ontarget**: OCF.ontarget uses **both** on-target fragments **AND** panel-filtered OCR regions. This dual-filter approach maximizes signal-to-noise for panel-specific tissue-of-origin detection. See [OCF Feature](../features/ocf.md#panel-mode---target-regions) for details.

**Note on TFBS/ATAC**: Off-target files use off-target reads across all ~800 TFs / 23 cancer types (WGS-like coverage). On-target files use only on-target reads. See [Region Entropy](../features/region-entropy.md#panel-mode-dual-output) for details.

**Off-target files** are used for:
- Fragment-based biomarkers
- GC-corrected coverage analysis
- Comparison with WGS baselines

**On-target files** are used for:
- Copy number analysis
- Integration with variant calling
- Panel-specific tissue signals (OCF)

## Target Regions File

The `--target-regions` BED file should contain the capture probe coordinates:

```
chr1    11166102    11166202    MTOR_exon1
chr1    27022522    27022622    ARID1A_exon1
...
```

- Use the **bait coordinates** (not target intervals)
- Standard BED format (0-based, half-open)
- Optional 4th column for region names


## Gene-Centric FSC (MSK-ACCESS)

For MSK-ACCESS panels (v1 and v2), krewlyzer provides **gene-level FSC aggregation**:

```bash
# FSC with gene-level output for MSK-ACCESS v2
krewlyzer fsc -i sample.bed.gz -o out/ --assay xs2
```

### Output Files

| File | Description |
|------|-------------|
| `{sample}.FSC.tsv` | Standard window-based FSC |
| `{sample}.FSC.gene.tsv` | Gene-level FSC (146 genes for xs2) |

### Gene FSC Output Format

```
gene    n_regions  total_bp  ultra_short  core_short  mono_nucl  di_nucl  long  total  ultra_short_ratio  ...
ATM     62         8432      1234         5678        9012       3456     789   20169  0.0612             ...
BRCA2   42         5689      ...
```

### Supported Assays

| Assay | Flag | Genes |
|-------|------|:-----:|
| MSK-ACCESS v1 | `--assay xs1` | 128 |
| MSK-ACCESS v2 | `--assay xs2` | 146 |

The gene groupings are bundled with krewlyzer in `data/genes/GRCh37/`.


## Panel WPS Anchors (MSK-ACCESS)

For MSK-ACCESS panels, WPS analysis uses **panel-specific anchors** filtered to genes in the panel:

```bash
# WPS with panel-specific anchors for MSK-ACCESS v2
krewlyzer wps -i sample.bed.gz -o out/ \
    --wps-anchors $(python -c "from krewlyzer.core.wps_anchor_filter import get_bundled_wps_anchors; print(get_bundled_wps_anchors('xs2', 'GRCh37'))")
```

### Bundled Panel Anchors

| Assay | File | Anchors |
|-------|------|:-------:|
| MSK-ACCESS v1 | `xs1.wps_anchors.bed.gz` | 1,611 |
| MSK-ACCESS v2 | `xs2.wps_anchors.bed.gz` | 1,820 |

### Anchor Types

- **TSS anchors**: Transcription start sites for panel genes
- **CTCF anchors**: CTCF binding sites within 100kb of panel TSS sites

> [!TIP] 
> Using panel-specific anchors reduces noise from irrelevant genome-wide signals and focuses WPS analysis on oncologically relevant regions.


## PON Compatibility

The PON model stores whether it was built in panel mode:

```python
from krewlyzer.pon.model import PonModel

pon = PonModel.load("msk-access.pon.parquet")
print(f"Panel mode: {pon.panel_mode}")
print(f"Target file: {pon.target_regions_file}")
```

For best results, use a PON built with the same `--target-regions` as sample processing.
