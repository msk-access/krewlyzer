# Krewlyzer Scripts

Utility scripts for generating bundled asset files.

> **Note**: These scripts are for **developers/maintainers** who need to regenerate bundled assets. Regular users don't need to run these.

---

## Scripts

### `convert_alu.py`

Filters RepeatMasker data to extract full-length Alu elements for WPS background stacking.

**Purpose**: Generate `hg19.alu_consensus.bed.gz` for WPS background analysis.

**Input**: UCSC RepeatMasker BED file
**Output**: BED7 with subfamily classification

```bash
python scripts/convert_alu.py \
    --input hg19_rmsk.bed.gz \
    --output data/WpsBackground/GRCh37/hg19.alu_consensus.bed
```

**What it does**:
- Filters for Alu elements only (AluY, AluS, AluJ families)
- Keeps full-length elements (280-320bp)
- Adds subfamily column for hierarchical stacking
- Outputs BED7: `chrom, start, end, name, score, strand, subfamily`

---

### `convert_anchors.py`

Converts Ensembl BioMart TSS and CTCF files to WPS anchor BED format.

**Purpose**: Generate `hg19.wps_anchors.bed.gz` for WPS foreground analysis.

**Input**: 
- TSS file from Ensembl BioMart (TSV)
- Regulatory features file from Ensembl BioMart (TSV)

**Output**: Merged BED6 with TSS + CTCF anchors

```bash
python scripts/convert_anchors.py \
    --tss ensembl_tss.tsv.gz \
    --ctcf ensembl_regulatory.tsv.gz \
    --output data/WpsAnchors/GRCh37/hg19.wps_anchors.bed
```

**What it does**:
- Selects canonical TSS (longest transcript per gene)
- Extracts CTCF binding site midpoints
- Outputs BED6: `chrom, start, end, name, score, strand`
- Names: `TSS|GeneName|TranscriptID` or `CTCF|chrom:position`

---

## Data Sources

### RepeatMasker (for Alu)

Download from UCSC Table Browser:
- Assembly: hg19 or hg38
- Group: Repeats
- Track: RepeatMasker
- Table: rmsk
- Output format: BED

### TSS (for WPS Anchors)

Download from Ensembl BioMart:
- Dataset: Human genes (GRCh37)
- Attributes: Gene stable ID, Transcript stable ID, Gene name, Chromosome, TSS, Strand, Transcript length

### CTCF (for WPS Anchors)

Download from Ensembl BioMart:
- Dataset: Human Regulatory Features (GRCh37)
- Filter: Feature type = "CTCF Binding Site"
- Attributes: Chromosome, Start, End, Feature type

---

## When to Regenerate

You typically **don't need** to regenerate these assets. Only run these if:

1. Adding support for a new genome build
2. Updating to a new Ensembl release
3. Fixing issues with anchor/Alu definitions

The bundled assets in `src/krewlyzer/data/` are pre-generated and shipped with the package.
