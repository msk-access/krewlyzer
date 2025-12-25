# Orientation-aware Fragmentation (OCF)

**Command**: `krewlyzer ocf`

## Purpose
Computes orientation-aware cfDNA fragmentation (OCF) values in tissue-specific open chromatin regions.

## Biological Context
OCF ([Sun et al., 2019](../citation.md#ocf)) measures the phasing of upstream (U) and downstream (D) fragment ends in open chromatin, informing tissue-of-origin of cfDNA.

## Usage
```bash
krewlyzer ocf sample.bed.gz --output output_dir/ [options]
```
## Output
- `{sample}.OCF.tsv`: Summary of OCF calculations per tissue type.
- `{sample}.OCF.sync.tsv`: Detailed sync scores.

## Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | *required* | Output directory |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--ocr-input` | `-r` | PATH | | Open chromatin regions file |
| `--genome` | `-G` | TEXT | hg19 | Genome build (hg19/hg38) |
| `--gc-correct` | | FLAG | True | Apply GC bias correction |
| `--verbose` | `-v` | FLAG | | Enable verbose logging |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |

## Calculation Details

1.  **Alignment**: Fragments are mapped relative to the center of the Open Chromatin Region (OCR).
2.  **Counting**:
    - `Left` ends (Start) and `Right` ends (End) are counted in 10bp bins across a ±1000bp window.
    - Counts are normalized by total sequencing depth.
3.  **OCF Score**:
    $$ OCF = \sum_{Peak} P_{signal} - \sum_{Peak} P_{background} $$
    - **Signal**: Right ends at -60bp and Left ends at +60bp (Phased nucleosome boundaries).
    - **Background**: Left ends at -60bp and Right ends at +60bp (Unphased).

## Clinical Interpretation

### Healthy Plasma Baseline
In healthy individuals, cfDNA primarily originates from:

| Tissue | OCF Value |
|--------|-----------|
| **T-cells (hematopoietic)** | Highest |
| **Liver** | Second highest |
| Other tissues | Near zero |

### Detecting Tumor-Derived cfDNA
When comparing a sample to healthy plasma:

| Pattern | Interpretation |
|---------|----------------|
| ↑ Tissue-specific OCF | Tumor shedding from that tissue |
| ↓ T-cell OCF | Dilution of hematopoietic cfDNA by tumor DNA |
| OCF correlates with tumor fraction | Higher ctDNA → stronger tissue signal |

### Cancer-Specific Patterns
| Cancer Type | Expected OCF Change |
|-------------|---------------------|
| Hepatocellular carcinoma | ↑ Liver OCF |
| Colorectal cancer | ↑ Intestine OCF, ↓ T-cell OCF |
| Lung cancer | ↑ Lung OCF, ↓ T-cell OCF |

> **Reference:** See [Citation & Scientific Background](../citation.md#ocf) for detailed paper summary.
