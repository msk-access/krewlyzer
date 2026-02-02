# Nextflow Pipeline

Run Krewlyzer at scale with the Nextflow pipeline.

## Quick Start

```bash
nextflow run msk-access/krewlyzer \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --outdir results/
```

## Workflow Architecture

The pipeline uses a Nextflow-native parallel pattern:

```mermaid
flowchart TB
    BAM["sample.bam"] --> EXTRACT["KREWLYZER_EXTRACT"]
    EXTRACT --> BED["sample.bed.gz"]
    
    BED --> MOTIF["KREWLYZER_MOTIF"]
    BED --> FSC["KREWLYZER_FSC"]
    BED --> FSD["KREWLYZER_FSD"]
    BED --> WPS["KREWLYZER_WPS"]
    BED --> OCF["KREWLYZER_OCF"]
    BED --> ENTROPY["KREWLYZER_REGION_ENTROPY"]
    BED --> RMDS["KREWLYZER_REGION_MDS"]
    
    FSC --> FSR["KREWLYZER_FSR"]
    
    subgraph "Parallel Paths"
        METH_BAM["meth.bam"] --> UXM["KREWLYZER_UXM"]
        BAM2["BAM + MAF"] --> MFSD["KREWLYZER_MFSD"]
    end
```

## Documentation

| Page | Description |
|------|-------------|
| [Samplesheet](samplesheet.md) | Input samplesheet format |
| [Parameters](parameters.md) | All pipeline parameters |
| [Outputs](outputs.md) | Output channels and files |
| [Examples](examples.md) | Workflow examples |

## Features

- **Parallel processing** - Process multiple samples simultaneously
- **Resume support** - Resume failed runs
- **Container support** - Docker/Singularity
- **Cloud ready** - AWS, Google Cloud, Azure

## See Also

- [CLI Reference](../cli/index.md) - Command-line usage
- [Panel Mode](../guides/panel-mode.md) - MSK-ACCESS workflows

