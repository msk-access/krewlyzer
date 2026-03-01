# CLI Reference

Krewlyzer provides a command-line interface for all feature extraction tools.

## Main Command

```bash
krewlyzer run-all -i sample.bam --reference hg19.fa --output results/
```

See [run-all](run-all.md) for the unified command that runs all features.

## Individual Commands

| Command | Description |
|---------|-------------|
| `extract` | Extract fragments from BAM |
| `fsc` | Fragment Size Coverage |
| `fsd` | Fragment Size Distribution |
| `fsr` | Fragment Size Ratio |
| `wps` | Windowed Protection Score |
| `motif` | End Motif extraction |
| `ocf` | Orientation-aware Fragmentation |
| `region-entropy` | TFBS/ATAC entropy |
| `region-mds` | Per-gene MDS |
| `mfsd` | Mutant Fragment Size Distribution |
| `uxm` | Fragment-level Methylation |
| `build-pon` | Build Panel of Normals |
| `build-gc-reference` | Build GC reference assets |

## See Also

- [Nextflow Pipeline](../nextflow/index.md) - Batch processing
- [Features](../features/index.md) - Feature documentation
