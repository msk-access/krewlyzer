# Pipeline Integration

## Run All Features

**Command**: `krewlyzer run-all`

The unified pipeline executes all feature extraction modules for a single BAM file in one optimized pass.

---

## Architecture Flowchart

```mermaid
flowchart TB
    BAM["sample.bam"] --> EXTRACT["extract"]
    REF["Reference FASTA"] --> EXTRACT
    
    EXTRACT --> BED["sample.bed.gz"]
    EXTRACT --> MOTIF["EndMotif + MDS"]
    EXTRACT --> GC["GC Correction Factors"]
    
    BED --> PIPELINE["Unified Rust Pipeline"]
    GC --> PIPELINE
    
    PIPELINE --> FSC["FSC.tsv"]
    PIPELINE --> FSR["FSR.tsv"]
    PIPELINE --> FSD["FSD.tsv"]
    PIPELINE --> WPS["WPS.parquet"]
    PIPELINE --> OCF["OCF.tsv"]
    
    subgraph "With --variants"
        BAM --> MFSD["mFSD.tsv"]
        VCF["variants.vcf/maf"] --> MFSD
    end
    
    subgraph "With --assay (Panel Mode)"
        PIPELINE --> FSCG["FSC.gene.tsv"]
        FSCG --> FSCR["FSC.regions.tsv"]
        FSCR --> E1["FSC.regions.e1only.tsv"]
        BAM --> RMDS["Region-MDS"]
        RMDS --> MDSE["MDS.exon.tsv"]
        RMDS --> MDSG["MDS.gene.tsv"]
    end
    
    subgraph "With --target-regions"
        TARGETS["Target BED"] --> PIPELINE
        PIPELINE --> ON["*.ontarget.tsv files"]
    end
```

### Python/Rust Boundary

```mermaid
flowchart TB
    subgraph "Python Layer"
        CLI["run-all CLI (wrapper.py)"]
        UP["unified_processor.py"]
        PROC["*_processor.py (post-processing)"]
        PON["PonModel (normalization)"]
        ASSETS["AssetManager"]
    end
    
    subgraph "Rust Layer (_core)"
        EXTRACT["extract_fragments()"]
        UNIFIED["run_unified_pipeline()"]
        GC["GC correction"]
        FSC_R["FSC/FSR counting"]
        FSD_R["FSD per arm"]
        WPS_R["WPS profiling"]
        OCF_R["OCF calculation"]
    end
    
    CLI --> ASSETS
    CLI --> UP
    UP --> UNIFIED
    
    UNIFIED --> GC --> FSC_R & FSD_R & WPS_R & OCF_R
    
    FSC_R --> PROC
    FSD_R --> PROC
    WPS_R --> PROC
    OCF_R --> PROC
    
    PROC --> PON
```

---

## Usage

```bash
# Basic run-all
krewlyzer run-all -i sample.bam -r hg19.fa -o output/

# With variant analysis
krewlyzer run-all -i sample.bam -r hg19.fa -o output/ \
    --variants mutations.maf

# Panel data (MSK-ACCESS) with --assay for gene-centric output
krewlyzer run-all -i sample.bam -r hg19.fa -o output/ \
    --target-regions panel_targets.bed \
    --assay xs2

# Full options
krewlyzer run-all -i sample.bam -r hg19.fa -o output/ \
    --variants mutations.vcf \
    --target-regions targets.bed \
    --assay xs2 \
    --pon-model cohort.pon.parquet \
    --threads 8
```

---

## CLI Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | *required* | Input BAM file (sorted, indexed) |
| `--reference` | `-r` | PATH | *required* | Reference genome FASTA |
| `--output` | `-o` | PATH | *required* | Output directory |
| `--genome` | `-G` | TEXT | hg19 | Genome build (hg19/hg38) |
| `--mapq` | `-q` | INT | 20 | Minimum mapping quality |
| `--minlen` | | INT | 65 | Minimum fragment length |
| `--maxlen` | | INT | 1000 | Maximum fragment length (extended FSD range) |
| `--skip-duplicates` | | FLAG | True | Skip duplicate reads |
| `--require-proper-pair` | | FLAG | True | Require proper pairs |
| `--exclude-regions` | `-x` | PATH | | Exclude regions BED file |
| `--target-regions` | `-T` | PATH | | Target BED (panel mode: GC from off-target) |
| `--assay` | `-A` | TEXT | | Assay code (xs1, xs2) for gene-centric FSC |
| `--bisulfite-bam` | | PATH | | Bisulfite BAM for UXM (optional) |
| `--variants` | `-v` | PATH | | VCF/MAF for mFSD analysis |
| `--duplex` | `-D` | FLAG | False | Enable duplex weighting for mFSD |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--chromosomes` | | TEXT | | Comma-separated chromosomes to process |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |
| `--arms-file` | `-a` | PATH | | Custom arms file for FSD |
| `--bin-input` | `-b` | PATH | | Custom bins for FSC/FSR |
| `--ocr-file` | | PATH | | Custom OCR file for OCF |
| `--wps-anchors` | | PATH | | WPS anchors BED for dual-stream profiling |
| `--wps-background` | | PATH | | WPS background Alu BED |
| `--bait-padding` | | INT | 50 | Bait edge padding in bp |
| `--pon-model` | `-P` | PATH | | PON model for normalization |
| `--skip-pon` | | FLAG | False | Skip PON z-score normalization |
| `--skip-target-regions` | | FLAG | False | Force WGS mode (ignore bundled targets) |
| `--no-tfbs` | | FLAG | False | Skip TFBS region entropy analysis |
| `--no-atac` | | FLAG | False | Skip ATAC region entropy analysis |
| `--disable-e1-aggregation` | | FLAG | False | Skip E1-only FSC region filtering |
| `--region-mds-e1-only` | | FLAG | False | Run region-MDS on E1 (first exon) only |
| `--generate-json` | | FLAG | False | Generate unified sample.features.json |
| `--debug` | | FLAG | False | Enable debug logging |
| `--validate-assets` | | FLAG | False | Validate bundled assets before running |

---

## Output Files

| File | Module | Description |
|------|--------|-------------|
| `{sample}.bed.gz` | extract | Fragments with GC content |
| `{sample}.EndMotif.tsv` | extract | 4-mer end motif frequencies |
| `{sample}.MDS.tsv` | extract | Motif Diversity Score |
| `{sample}.FSC.tsv` | fsc | 5-channel coverage per bin |
| `{sample}.FSC.gene.tsv` | fsc | Gene-centric FSC (with --assay) |
| `{sample}.FSR.tsv` | fsr | Short/long fragment ratios |
| `{sample}.FSD.tsv` | fsd | Size distribution per arm |
| `{sample}.WPS.parquet` | wps | Nucleosome protection profiles |
| `{sample}.WPS_background.parquet` | wps | Alu element stacking |
| `{sample}.OCF.tsv` | ocf | Tissue-of-origin OCF |
| `{sample}.mFSD.tsv` | mfsd | Mutant vs WT sizes (with -v) |
| `{sample}.TFBS.tsv` | region_entropy | TFBS size entropy (808 TFs) |
| `{sample}.ATAC.tsv` | region_entropy | ATAC size entropy (23 types) |
| `{sample}.MDS.exon.tsv` | region_mds | Per-exon MDS (with --assay) |
| `{sample}.MDS.gene.tsv` | region_mds | Gene-level MDS (with --assay) |
| `{sample}.features.json` | run-all | Unified JSON (with --generate-json) |

---

## Panel Mode (--target-regions)

For targeted sequencing panels (MSK-ACCESS):

```bash
krewlyzer run-all -i sample.bam -r hg19.fa -o output/ \
    --target-regions MSK-ACCESS_targets.bed \
    --bin-input MSK-ACCESS_genes.bed
```

### What Happens in Panel Mode

```mermaid
flowchart TB
    subgraph "GC Correction"
        EXTRACT["extract"] -->|"Off-target only"| GC_MODEL["Unbiased GC model"]
    end
    
    subgraph "Feature Output"
        PIPELINE["Each tool"] --> OFF["*.tsv (off-target)"]
        PIPELINE --> ON["*.ontarget.tsv (on-target)"]
    end
    
    subgraph "OCF Warning"
        OCF_TOOL["OCF"] -->|"hg38 + no custom file"| SKIP["Skipped with warning"]
    end
```

| Effect | Description |
|--------|-------------|
| GC model | Built from **off-target** fragments only |
| FSC/FSR aggregation | **Disabled** (preserves gene resolution) |
| All outputs | Split into `.tsv` (off) and `.ontarget.tsv` (on) |
| OCF (hg38) | Skipped unless custom OCR file provided |

> [!IMPORTANT]
> **Off-target = unbiased** – preferred for fragmentomics biomarkers.

---

## Nextflow Pipeline

For batch processing, use the Nextflow pipeline:

```bash
nextflow run main.nf \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --asset_dir /path/to/krewlyzer/data/ \
    --outdir results/
```

### Workflow Architecture

The pipeline uses a Nextflow-native parallel pattern:

```
BAM → EXTRACT → BED.gz
              ↓
    ┌───┬───┬───┼───┬───┬───────────┬────────────┐
    ↓   ↓   ↓   ↓   ↓   ↓           ↓            ↓
  MOTIF FSC FSD WPS OCF REGION_ENTROPY      REGION_MDS
         ↓
        FSR

Meth BAM → UXM         (parallel path)
BAM + MAF → MFSD       (parallel path)
```

**Available Modules (14 total):**
| Module | Description |
|--------|-------------|
| `KREWLYZER_EXTRACT` | Fragment extraction from BAM |
| `KREWLYZER_FSC` | Fragment Size Coverage |
| `KREWLYZER_FSR` | Fragment Size Ratio |
| `KREWLYZER_FSD` | Fragment Size Distribution |
| `KREWLYZER_WPS` | Windowed Protection Score |
| `KREWLYZER_OCF` | Orientation cfDNA Fragmentation |
| `KREWLYZER_MOTIF` | End Motif & MDS |
| `KREWLYZER_REGION_ENTROPY` | TFBS/ATAC size entropy |
| `KREWLYZER_REGION_MDS` | Per-gene/exon Motif Diversity |
| `KREWLYZER_UXM` | Methylation deconvolution |
| `KREWLYZER_MFSD` | Mutant Fragment Size |
| `KREWLYZER_RUNALL` | Full pipeline (single process) |
| `KREWLYZER_BUILD_PON` | Build PON model |
| `FILTER_MAF` | MAF filtering |

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | *required* | CSV with sample information |
| `--ref` | *required* | Reference genome FASTA |
| `--outdir` | `./results` | Output directory |
| `--asset_dir` | | Base directory for PON/targets (enables assay resolution) |
| `--targets` | | Global target BED (fallback) |
| `--genome` | `hg19` | Genome build |
| `--pon_model` | | Global PON model (fallback) |
| `--bait_padding` | `50` | Bait edge padding for WPS |
| `--skip_pon` | `false` | Skip PON z-score normalization |
| `--no_tfbs` | `false` | Disable TFBS region entropy |
| `--no_atac` | `false` | Disable ATAC region entropy |
| `--maxlen` | `1000` | Maximum fragment length (extended range) |
| `--verbose` | `false` | Enable verbose logging |
| `--threads` | `8` | Threads per process |

### Samplesheet Format

```csv
sample,bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
```

| Column | Type | Description |
|--------|------|-------------|
| `sample` | TEXT | Sample identifier |
| `bam` | PATH | WGS/Panel BAM file |
| `meth_bam` | PATH | Bisulfite BAM for UXM |
| `vcf` | PATH | VCF for mFSD |
| `bed` | PATH | Pre-extracted .bed.gz |
| `maf` | PATH | MAF for mFSD |
| `single_sample_maf` | BOOL | Skip MAF filtering if `true` |
| `assay` | TEXT | Assay code: `XS1`, `XS2`, `WGS` |
| `pon` | PATH | Sample-specific PON (overrides assay) |
| `targets` | PATH | Sample-specific targets (overrides assay) |

### Assay Resolution (XS1/XS2)

When `assay` is set and `--asset_dir` is provided, the pipeline auto-resolves:

| Assay | PON File | Targets File | Genes | WPS Anchors |
|-------|----------|--------------|:-----:|:-----------:|
| `XS1` | `{asset_dir}/pon/GRCh37/xs1.pon.parquet` | `{asset_dir}/targets/GRCh37/xs1.targets.bed` | 128 | 1,611 |
| `XS2` | `{asset_dir}/pon/GRCh37/xs2.pon.parquet` | `{asset_dir}/targets/GRCh37/xs2.targets.bed` | 146 | 1,820 |
| `WGS` | `{asset_dir}/pon/GRCh37/wgs.pon.parquet` | None | - | 147,772 |

### Example Samplesheet

```csv
sample,bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
# MSK-ACCESS V1 samples (auto-resolve PON/targets)
ACCESS_001,/data/sample1.bam,,,,,false,XS1,,

# MSK-ACCESS V2 with MAF
ACCESS_002,/data/sample2.bam,,,,/data/cohort.maf,false,XS2,,

# WGS (no targets)
WGS_001,/data/wgs.bam,,/data/wgs.vcf,,,,WGS,,

# Custom PON/targets override
CUSTOM,/data/custom.bam,,,,,,,/data/custom.pon.parquet,/data/custom.bed
```

### Pipeline Logic

| Input | Triggered Workflow |
|-------|-------------------|
| `bam` only | Full run-all (extract → features) |
| `bam` + `vcf`/`maf` | run-all + mFSD |
| `meth_bam` only | UXM methylation |
| `bed` only | FSC, FSR, FSD, WPS, OCF (no extract) |

### Profiles
- `-profile docker` – Docker container
- `-profile slurm` – SLURM clusters (cmobic_cpu)

---

## See Also

- [Feature Documentation](features/extract.md) – Per-tool details
- [Architecture](advanced/architecture.md) – Rust/Python structure
- [PON Models](advanced/pon.md) – Normalization baselines
