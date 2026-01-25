# Unified JSON Output

The `--generate-json` flag produces a single JSON file containing all features for ML integration.

## Enabling JSON Output

```bash
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ \
    --generate-json
```

This generates `{sample}.features.json` alongside the standard TSV/Parquet outputs.

## Output Structure

```json
{
  "sample_id": "sample_001",
  "metadata": {
    "genome": "hg19",
    "assay": "xs2",
    "panel_mode": true,
    "on_target_rate": 0.45,
    "timestamp": "2024-01-20T00:00:00"
  },
  "fsc": { ... },
  "fsc_gene": { ... },
  "fsc_region": { ... },
  "fsr": { ... },
  "fsd": { ... },
  "wps": { ... },
  "wps_panel": { ... },
  "wps_background": { ... },
  "motif": { ... },
  "ocf": { ... },
  "tfbs": { ... },
  "atac": { ... },
  "gc_factors": { ... }
}
```

---

## Feature Schemas

### FSC (Fragment Size Coverage)

Window-based fragment size counts with z-scores.

```json
"fsc": {
  "n_windows": 2534,
  "data": [
    {
      "region": "chr1:0-500000",
      "ultra_short": 123,
      "core_short": 4567,
      "mono_nucl": 8901,
      "di_nucl": 2345,
      "long": 678,
      "total": 16614,
      "log_ratio_core_short": -0.15,
      "log_ratio_mono_nucl": 0.02,
      "zscore_core_short": -1.23
    }
  ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `region` | string | Genomic window coordinates |
| `ultra_short` | int | Fragments 65-99bp |
| `core_short` | int | Fragments 100-149bp |
| `mono_nucl` | int | Fragments 150-259bp |
| `di_nucl` | int | Fragments 260-399bp |
| `long` | int | Fragments 400+bp |
| `log_ratio_*` | float | Log2(observed/expected) vs PON |
| `zscore_*` | float | Z-score vs PON (if PON provided) |

---

### FSC Gene (Panel Mode Only)

Gene-level fragment size aggregation. Only present with `--assay`.

```json
"fsc_gene": [
  {
    "gene": "ATM",
    "n_regions": 62,
    "total_bp": 8432,
    "ultra_short": 1234,
    "core_short": 5678,
    "mono_nucl": 9012,
    "core_short_ratio": 0.282,
    "normalized_depth": 1245.67,
    "z_core_short": -0.45
  }
]
```

---

### FSC Region (Panel Mode Only)

Per-exon/target fragment size data. More granular than gene-level.

```json
"fsc_region": [
  {
    "chrom": "1",
    "start": 11168235,
    "end": 11168345,
    "gene": "MTOR",
    "region_name": "MTOR_target_02",
    "region_bp": 110,
    "ultra_short": 8.0,
    "core_short": 229.0,
    "mono_nucl": 804.0,
    "di_nucl": 88.0,
    "total": 1129.0,
    "normalized_depth": 1272.71
  }
]
```

| Field | Type | Description |
|-------|------|-------------|
| `region_name` | string | Unique exon/target identifier |
| `normalized_depth` | float | RPKM-like depth: (count × 10⁹) / (bp × total_frags) |

### FSR (Fragment Size Ratios)

Biomarker ratios for tumor detection.

```json
"fsr": {
  "n_windows": 2534,
  "data": [
    {
      "region": "chr1:0-500000",
      "ultra_short_ratio": 0.0074,
      "core_short_ratio": 0.275,
      "mono_nucl_ratio": 0.536,
      "di_nucl_ratio": 0.141,
      "long_ratio": 0.041,
      "core_short_long_ratio": 6.73
    }
  ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `*_ratio` | float | Fraction of total fragments |
| `core_short_long_ratio` | float | Primary cancer biomarker (higher = more tumor) |

---

### FSD (Fragment Size Distribution)

Per-arm size distribution profiles.

```json
"fsd": {
  "arms": ["1p", "1q", "2p", ...],
  "size_bins": [65, 70, 75, ..., 395, 400],
  "data": {
    "1p": {
      "counts": [123, 456, 789, ...],
      "proportions": [0.001, 0.004, 0.007, ...]
    }
  }
}
```

---

### WPS (Windowed Protection Score)

Nucleosome positioning profiles around gene TSS/CTCF sites.

```json
"wps": {
  "n_anchors": 15234,
  "columns": ["region_id", "chrom", "start", "end", 
              "wps_nuc_mean", "wps_tf_mean", "prot_frac_nuc", "prot_frac_tf",
              "wps_nuc_z", "wps_tf_z", "ndr_depth"],
  "data": [
    {
      "region_id": "ENSG00000142611_TSS",
      "chrom": "chr1",
      "start": 11166102,
      "end": 11166502,
      "wps_nuc_mean": 24.5,
      "wps_tf_mean": -3.2,
      "prot_frac_nuc": 0.62,
      "prot_frac_tf": 0.41,
      "wps_nuc_z": 1.2,
      "ndr_depth": 15.3
    }
  ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `wps_nuc_mean` | float | Mean nucleosomal WPS (120-180bp fragments) |
| `wps_tf_mean` | float | Mean TF footprint WPS (35-80bp fragments) |
| `prot_frac_*` | float | Protected fraction (values > 0) |
| `wps_*_z` | float | Z-score vs PON |
| `ndr_depth` | float | Nucleosome-depleted region depth |

---

### WPS Panel (Panel Mode Only)

Same schema as `wps`, but filtered to panel gene anchors.

```json
"wps_panel": {
  "n_anchors": 1820,
  "data": [ ... ]
}
```

---

### WPS Background

Alu element stacking scores for global fragmentation.

```json
"wps_background": {
  "n_elements": 142567,
  "data": [
    {
      "region_id": "AluSx_chr1_12345",
      "stacking_score": 0.78,
      "coverage": 45.2
    }
  ]
}
```

---

### Motif

End motif (EDM) k-mer frequencies and diversity score.

```json
"motif": {
  "end_motif": {
    "AAAA": 0.0042,
    "AAAC": 0.0039,
    ...
  },
  "breakpoint_motif": {
    "AAAA": 0.0038,
    ...
  },
  "mds": 0.8234,
  "mds_z": -1.23
}
```

| Field | Type | Description |
|-------|------|-------------|
| `end_motif` | dict | 256 4-mer frequencies (fragment ends) |
| `breakpoint_motif` | dict | 256 4-mer frequencies (breakpoints) |
| `mds` | float | Motif Diversity Score |
| `mds_z` | float | MDS z-score vs PON (if PON with MDS baseline) |

---

### OCF (Orientation-aware cfDNA Fragmentation)

Open chromatin footprint scores by tissue type.

```json
"ocf": {
  "tissues": ["Liver", "Lung", "Colon", "Placenta", ...],
  "scores": {
    "Liver": 0.42,
    "Lung": 0.31,
    "Colon": 0.28
  }
}
```

---

### TFBS (Transcription Factor Binding Site Entropy)

Fragment size entropy at TFBS regions.

```json
"tfbs": {
  "off_target": [
    {
      "region": "CTCF_chr1_12345",
      "entropy": 3.45,
      "n_fragments": 234,
      "mean_size": 167.5
    }
  ],
  "on_target": [ ... ]
}
```

---

### ATAC (Chromatin Accessibility Regions)

Fragment size entropy at ATAC-seq accessible regions.

```json
"atac": {
  "off_target": [
    {
      "region": "peak_chr1_23456",
      "entropy": 3.21,
      "n_fragments": 189,
      "mean_size": 145.2
    }
  ],
  "on_target": [ ... ]
}
```

---

### GC Factors (Diagnostic)

GC bias correction factors used internally during processing.

> [!NOTE]
> **Not recommended for ML features.** These are intermediate diagnostic data, not predictive features.
> The GC correction is already applied to FSC/FSR/FSD counts. Use those corrected values instead.

```json
"gc_factors": {
  "off_target": [
    {
      "len_bin": 100,
      "gc_pct": 45,
      "correction_factor": 1.12
    }
  ],
  "on_target": [ ... ]
}
```

**When to use GC factors:**
- **QC/Diagnostics**: Visualize library prep bias, capture efficiency
- **Batch Effect Detection**: Compare correction factors across runs
- **Panel Development**: Characterize probe GC performance

**When NOT to use:**
- **ML models**: Skip these—use GC-corrected FSC/FSR/FSD instead

---

## ML Integration Example

```python
import json
import pandas as pd

# Load features
with open("sample.features.json") as f:
    features = json.load(f)

# Extract FSC gene-level for panel analysis
if "fsc_gene" in features:
    df_genes = pd.DataFrame(features["fsc_gene"])
    print(f"Gene FSC: {len(df_genes)} genes")

# Extract WPS for nucleosome signature
wps_data = pd.DataFrame(features["wps"]["data"])
print(f"WPS anchors: {len(wps_data)}")

# Use motif MDS z-score as feature
mds_z = features["motif"].get("mds_z", 0)
print(f"MDS z-score: {mds_z:.2f}")
```

---

## JSON Generation with Panel Mode

For MSK-ACCESS panels:

```bash
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ \
    --assay xs2 \
    --target-regions targets.bed \
    --pon-model xs2.pon.parquet \
    --generate-json
```

This produces JSON with all panel-specific features:
- `fsc_gene`: 146 genes
- `wps_panel`: 1,820 anchors
- `wps_background`: Alu stacking
- PON z-scores across all features
