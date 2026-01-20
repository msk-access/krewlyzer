# Getting Started

Get running with Krewlyzer in 5 minutes.

## Quick Install

=== "pip"
    ```bash
    pip install krewlyzer
    ```

=== "Docker"
    ```bash
    docker pull ghcr.io/msk-access/krewlyzer:latest
    ```

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
# sample.bed.gz            # Extracted fragments
# sample.EndMotif.tsv      # End motif frequencies
# sample.FSC.tsv           # Fragment size coverage
# sample.FSR.tsv           # Fragment size ratios
# sample.FSD.tsv           # Size distribution by arm
# sample.WPS.tsv.gz        # Windowed protection scores
# sample.OCF.tsv           # Orientation-aware fragmentation
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

**→ [Full MSK-ACCESS Quickstart](advanced/msk-access-quickstart.md)** for detailed workflows.

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
- [Usage Guide](usage.md) - Full CLI reference
- [Feature Documentation](features/extract.md) - Per-feature details
- [Nextflow Pipeline](pipeline.md) - Batch processing
