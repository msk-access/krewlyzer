# Getting Started

Get running with Krewlyzer in 5 minutes.

## Quick Install

=== "Docker (Recommended)"
    ```bash
    docker pull ghcr.io/msk-access/krewlyzer:X.Y.Z  # Replace X.Y.Z with latest release version
    ```

=== "Clone + Install"
    ```bash
    git clone https://github.com/msk-access/krewlyzer.git && cd krewlyzer
    git lfs pull && pip install -e .
    ```

=== "pip + Data Clone"
    ```bash
    pip install krewlyzer
    git clone --depth 1 https://github.com/msk-access/krewlyzer.git ~/.krewlyzer-data
    cd ~/.krewlyzer-data && git lfs pull
    export KREWLYZER_DATA_DIR=~/.krewlyzer-data/src/krewlyzer/data
    ```

!!! note "pip Users"
    See [Installation Guide](installation.md) for `KREWLYZER_DATA_DIR` setup.

## Your First Analysis

Run all fragmentomics features on a BAM file:

```bash
krewlyzer run-all sample.bam \
    --reference hg19.fa \
    --output results/
```

## Check Results

```bash
ls results/
# sample.bed.gz                    # Extracted fragments
# sample.metadata.json             # Run metadata and QC metrics
# sample.correction_factors.tsv    # GC correction factors
# sample.EndMotif.tsv              # End motif frequencies
# sample.EndMotif1mer.tsv          # 1-mer motifs + Jagged Index
# sample.BreakPointMotif.tsv       # Breakpoint motif frequencies
# sample.MDS.tsv                   # Motif Diversity Score
# sample.FSC.tsv                   # Fragment size coverage
# sample.FSR.tsv                   # Fragment size ratios
# sample.FSD.tsv                   # Size distribution by arm
# sample.WPS.parquet               # Windowed protection scores
# sample.WPS_background.parquet    # WPS background (Alu)
# sample.OCF.tsv                   # Orientation-aware fragmentation
# sample.TFBS.tsv                  # TFBS entropy (808 TFs)
# sample.ATAC.tsv                  # ATAC entropy (23 cancer types)
# sample.features.json             # Unified JSON for ML
```

## Common Workflows

### Targeted Panel (MSK-ACCESS)

For MSK-ACCESS v1 or v2 panels, use the `--assay` flag for panel-optimized analysis:

```bash
krewlyzer run-all sample.bam \
    --reference hg19.fa \
    --output results/ \
    --assay xs2 \
    --target-regions MSK-ACCESS-v2_targets.bed
```

This enables:
- **Gene-level FSC** aggregation (146 genes)
- **Dual WPS output** (genome-wide + panel-specific)
- **On/off-target splitting** for all features

For PON normalization:
```bash
krewlyzer run-all sample.bam \
    --reference hg19.fa \
    --output results/ \
    --assay xs2 \
    --target-regions targets.bed \
    --pon-model xs2.pon.parquet \
    --generate-json
```

**→ [Full MSK-ACCESS Quickstart](../guides/msk-access-quickstart.md)** for detailed workflows.

### With Variant Analysis

Add mutant fragment size analysis using a VCF/MAF:

```bash
krewlyzer run-all sample.bam \
    --reference hg19.fa \
    --output results/ \
    --variants variants.maf
```

### Individual Tools

Run specific features separately:

```bash
# Extract fragments first
krewlyzer extract -i sample.bam -r hg19.fa -o output/

# Then run feature tools on the BED
krewlyzer fsc -i output/sample.bed.gz -o output/
krewlyzer wps -i output/sample.bed.gz -o output/
```

## Next Steps

- [Installation Guide](installation.md) - Detailed setup instructions
- [Usage Guide](../cli/run-all.md) - Full CLI reference
- [Feature Documentation](../features/core/extract.md) - Per-feature details
- [Nextflow Pipeline](../nextflow/index.md) - Batch processing
