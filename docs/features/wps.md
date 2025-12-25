# Windowed Protection Score (WPS)

**Command**: `krewlyzer wps`

## Purpose
Computes unified nucleosome and transcription factor protection scores for each region in a transcript/region file.

## Biological Context
The WPS (Snyder et al., 2016) quantifies chromatin accessibility by comparing fragments that span a protection window to those ending within it. Krewlyzer calculates **both** Long and Short WPS in a single pass:

- **Long WPS** (120-180bp fragments): Measures **nucleosome positioning**
- **Short WPS** (35-80bp fragments): Measures **transcription factor binding footprints**
- **WPS Ratio** (Long/Short): Balance between nucleosome and TF protection

## Usage
```bash
krewlyzer wps sample.bed.gz -o output_dir/ --sample-name SAMPLE [options]
```

## Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | *required* | Output directory |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--tsv-input` | | PATH | | Transcript annotation TSV |
| `--pon-model` | `-P` | PATH | | PON model for hybrid correction |
| `--genome` | `-G` | TEXT | hg19 | Genome build (hg19/hg38) |
| `--reference` | `-r` | PATH | | Reference FASTA for GC correction |
| `--gc-correct` | | FLAG | | Apply GC bias correction |
| `--empty` | | FLAG | | Include regions with no coverage |
| `--verbose` | `-v` | FLAG | | Enable verbose logging |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |


## Output Format

Output file: `{sample}.WPS.tsv.gz`

| Column | Description |
|--------|-------------|
| `gene_id` | Transcript/region identifier |
| `chrom` | Chromosome |
| `pos` | Genomic position |
| `cov_long` | Coverage from long fragments (120-180bp) |
| `cov_short` | Coverage from short fragments (35-80bp) |
| `wps_long` | WPS for long fragments (nucleosome) |
| `wps_short` | WPS for short fragments (TF) |
| `wps_ratio` | Ratio of wps_long / wps_short |
| `wps_long_norm` | Normalized WPS long (per million fragments) |
| `wps_short_norm` | Normalized WPS short (per million fragments) |
| `wps_ratio_norm` | Normalized WPS ratio |

## Column Interpretation Guide

### Coverage Columns (`cov_long`, `cov_short`)

| Value | Interpretation |
|-------|----------------|
| **0** | No fragments cover this exact position |
| **Low (1-5)** | Sparse coverage, low confidence |
| **Medium (5-20)** | Adequate coverage for reliable WPS |
| **High (>20)** | Well-covered, high-confidence region |

> **Note**: Coverage counts fragments covering the **exact position**, while WPS considers fragments in the **protection window** around the position. Thus it's possible to have WPS ≠ 0 with coverage = 0.

### WPS Long (`wps_long`) - Nucleosome Positioning

| Value | Interpretation | Biological Meaning |
|-------|----------------|-------------------|
| **> 0** | Protected | Nucleosome present, stable chromatin |
| **≈ 0** | Neutral | Transitional region |
| **< 0** | Fragmented | Nucleosome-depleted, accessible chromatin |

### WPS Short (`wps_short`) - Transcription Factor Footprints

| Value | Interpretation | Biological Meaning |
|-------|----------------|-------------------|
| **> 0** | TF-protected | Transcription factor bound at position |
| **≈ 0** | No TF signal | No detectable TF footprint |
| **< 0** | TF-fragmented | Active TF turnover or open region |

### WPS Ratio (`wps_ratio`) - Nucleosome vs TF Balance

| Value | Interpretation | Biological Meaning |
|-------|----------------|-------------------|
| **> 1** | Nucleosome-dominated | Stable, closed chromatin |
| **≈ 1** | Balanced | Mixed nucleosome/TF regulation |
| **< 1** | TF-dominated | Actively regulated region |
| **negative** | Complex signal | Fragmented for both (very accessible) |

### When Coverage = 0 but WPS ≠ 0

This occurs when fragments overlap the **protection window** around a position but don't directly cover that position:

```
Position:        |------ pos 100 ------|
WPS Window:      [======40-160========]
Fragment:    [===10-50===]
                       ↑
                 Overlaps window → affects WPS
                 But doesn't reach pos 100 → cov = 0
```

## Cancer Applications

- **Tumor Detection**: Altered nucleosome patterns at promoters indicate aberrant gene regulation
- **Tissue of Origin**: Cell-type-specific nucleosome footprints reveal tumor origin
- **TF Activity**: Short WPS can detect cancer-specific TF binding (e.g., ER in breast cancer)

## Calculation Details

For each position $k$ in the region of interest:

$$ WPS(k) = N_{spanning}(k) - N_{ends}(k) $$

Where:
- $N_{spanning}(k)$ is the count of fragments that completely overlap the window $[k - P, k + P]$.
- $N_{ends}(k)$ is the count of fragments with an endpoint within the window $[k - P, k + P]$.
- **Protection Window ($P$)**:
    - **Long WPS**: $P=60$ bp (Total window ~120bp). Targets 120-180bp fragments.
    - **Short WPS**: $P=8$ bp (Total window ~16bp). Targets 35-80bp fragments.

## Depth Normalization

The normalized columns (`wps_long_norm`, `wps_short_norm`, `wps_ratio_norm`) enable comparison of WPS values **across samples** with different sequencing depths.

### Formula

```
wps_norm = wps_raw / (total_fragments / 1,000,000)
```

This gives **per-million fragment** values.

### Metadata Requirement

Normalization requires a metadata file from the `extract` command:

```
{sample}.metadata.json  ← Contains total_fragments
{sample}.bed.gz         ← Input to wps
```

> **Warning**: If metadata is missing, normalized columns default to raw values (not comparable across samples). Run `krewlyzer extract` first.

### When to Use Raw vs Normalized

| Use Case | Columns |
|----------|---------|
| **Within-sample** pattern analysis | `wps_long`, `wps_short`, `wps_ratio` |
| **Cross-sample** comparison | `wps_long_norm`, `wps_short_norm` |
| **Machine learning** across cohort | Normalized columns |

## References

> See [Citation & Scientific Background](../citation.md#wps) for detailed paper summary.
