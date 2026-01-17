# Output Formats

Krewlyzer supports multiple output formats to suit different use cases: ML pipelines, human inspection, and data analysis.

## Format Summary

| Format | Extension | Use Case | Default For |
|--------|-----------|----------|-------------|
| **TSV** | `.tsv` | Human-readable, Excel | FSD, FSR, FSC, Motif, OCF |
| **Parquet** | `.parquet` | ML training, pandas | WPS (vector data), PON |
| **JSON** | `.json` | APIs, ML pipelines | Unified output |

## Command-Line Flags

### `run-all` Global Options

```bash
# Generate unified JSON with ALL features
krewlyzer run-all -i sample.bam -o out/ --generate-json

# Force all outputs to a specific format
krewlyzer run-all -i sample.bam -o out/ --output-format parquet
```

| Flag | Description |
|------|-------------|
| `--generate-json` | Create `sample.features.json` with all data |
| `--output-format` | Global format: `auto`, `tsv`, `parquet`, `json` |

### Individual Tool Override

Each tool supports `--format/-f` to override the default:

```bash
krewlyzer fsd -i sample.bed.gz -o out/ --format json
krewlyzer wps -i sample.bed.gz -o out/ --format tsv
```

## Unified JSON Output

The `--generate-json` flag creates a single JSON file with **all features** for ML pipelines:

```
out/sample.features.json
```

### Schema (v1.0)

```json
{
  "schema_version": "1.0",
  "sample_id": "sample_001",
  "krewlyzer_version": "0.3.2",
  "timestamp": "2026-01-17T14:52:33Z",
  
  "metadata": {
    "total_fragments": 8200000,
    "on_target_rate": 0.41,
    "gc_bias_corrected": true
  },
  
  "features": {
    "fsd": {
      "arms": ["1p", "1q", ...],
      "size_bins": ["65-69", "70-74", ...],
      "counts": [[120, 135, ...], ...]
    },
    "fsr": [{...}, ...],
    "wps": {"regions": [...], "wps_nuc": [...], ...},
    "motif": {"mds": 0.87, "edm": {...}, "bpm": {...}},
    "ocf": [{...}, ...]
  },
  
  "qc": {...}
}
```

### Features Included

| Feature | Data Included |
|---------|---------------|
| **FSD** | Full matrix (41 arms × size bins) |
| **FSR** | All regions with ratios |
| **FSC** | All windows with z-scores |
| **WPS** | Full vectors per region |
| **Motif** | All 256 k-mer frequencies + MDS |
| **OCF** | All regions with scores |
| **UXM** | If enabled (`--bisulfite-bam`) |
| **mFSD** | If enabled (`--variants`) |

## Smart Defaults

| Tool | Default Format | Reason |
|------|:--------------:|--------|
| fsd, fsr, fsc | TSV | Human-readable tables |
| wps | Parquet | Vector data (200+ values) |
| motif, ocf | TSV | Small tables |
| PON model | Parquet | Efficient storage |

## File Extensions

After implementation, output files follow this pattern:

```
out/
├── sample.FSD.tsv
├── sample.FSR.tsv
├── sample.FSC.tsv
├── sample.WPS.parquet
├── sample.EndMotif.tsv
├── sample.MDS.tsv
├── sample.OCF.tsv
├── sample.correction_factors.tsv  # Was .csv
├── sample.metadata.json
└── sample.features.json  # With --generate-json
```

## Migration Notes

### Breaking Change: `.csv` → `.tsv`

The `correction_factors.csv` file is now `correction_factors.tsv` for consistency. Update any downstream scripts:

```diff
- sample.correction_factors.csv
+ sample.correction_factors.tsv
```
