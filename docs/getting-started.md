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

For targeted sequencing, provide a custom bins file:

```bash
krewlyzer run-all sample.bam \
    --reference hg19.fa \
    --output results/ \
    --bin-input MSK-ACCESS-targets.bed
```

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
