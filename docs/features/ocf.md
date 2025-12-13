# Orientation-aware Fragmentation (OCF)

**Command**: `krewlyzer ocf`

## Purpose
Computes orientation-aware cfDNA fragmentation (OCF) values in tissue-specific open chromatin regions.

## Biological Context
OCF (Sun et al., Genome Res 2019) measures the phasing of upstream (U) and downstream (D) fragment ends in open chromatin, informing tissue-of-origin of cfDNA.

## Usage
```bash
krewlyzer ocf sample.bed.gz --output output_dir/ [options]
```
## Output
- `{sample}.OCF.tsv`: Summary of OCF calculations per tissue type.
- `{sample}.OCF.sync.tsv`: Detailed sync scores.

## Options
- `--threads`, `-t`: Number of processes

## Calculation Details

1.  **Alignment**: Fragments are mapped relative to the center of the Open Chromatin Region (OCR).
2.  **Counting**:
    - `Left` ends (Start) and `Right` ends (End) are counted in 10bp bins across a ±1000bp window.
    - Counts are normalized by total sequencing depth.
3.  **OCF Score**:
    $$ OCF = \sum_{Peak} P_{signal} - \sum_{Peak} P_{background} $$
    - **Signal**: Right ends at -60bp and Left ends at +60bp (Phased nucleosome boundaries).
    - **Background**: Left ends at -60bp and Right ends at +60bp (Unphased).

## Interpretation
- **High OCF**: Fragment ends align perfectly with nucleosome boundaries flanking the OCR, indicating the region is **open** in the tissue of origin.
- **Low OCF**: No phasing, indicating the region is **closed**.

