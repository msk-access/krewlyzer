# Input File Formats

This page documents the expected formats for custom input files used as overrides. Krewlyzer validates these formats when you provide custom files.

## Quick Reference

| File Type | Columns | Used By | Example |
|-----------|---------|---------|---------|
| [Sample List](#sample-list) | paths | `build-pon` | `/path/to/sample.bam` |
| [BED3](#bed3) | chrom, start, end | `--bin-input`, `--target-regions` | `chr1\t0\t100000` |
| [Gene BED](#gene-bed) | chrom, start, end, gene, [name] | `--gene-bed` | `chr1\t100\t5000\tTP53\texon1` |
| [Arms BED](#arms-bed) | chrom, start, end, arm | `--arms-file` | `chr1\t0\t125000000\t1p` |
| [WPS Anchors](#wps-anchors) | BED6 format | `--wps-anchors`, `--wps-background` | `chr1\t1000\t2000\tGene_TSS\t0\t+` |
| [Region BED](#region-bed) | chrom, start, end, label | `--ocr-file`, `--tfbs-regions`, `--atac-regions` | `chr1\t500\t800\tLiver` |
| [GC Factors TSV](#gc-factors-tsv) | length_bin, gc_pct, factor | `--gc-factors` | `10\t45\t1.05` |

---

## Sample List {#sample-list}

Plain text file with one sample path per line for PON building.

### Format

```
/path/to/sample1.bam
/path/to/sample2.bam
/path/to/sample3.bed.gz
```

| Input Type | Description |
|------------|-------------|
| `.bam` / `.cram` | Full processing including MDS baseline |
| `.bed.gz` | Pre-extracted fragments (faster, no MDS) |

### Notes

- One path per line
- No header row
- Paths can be absolute or relative to working directory
- Mixing BAM and BED.gz inputs is allowed

### Used By

- `build-pon SAMPLE_LIST` - First positional argument

---

<a id="bins"></a>
## BED3 {#bed3}

Standard 3-column BED format for genomic intervals.

### Format

```
chrom    start    end
```

| Column | Type | Description |
|--------|------|-------------|
| chrom | string | Chromosome (e.g., `chr1`, `chrX`) |
| start | int | 0-based start position |
| end | int | 1-based end position (exclusive) |

### Example

```tsv
chr1	0	100000
chr1	100000	200000
chr2	0	100000
```

### Used By

- `--bin-input` / `-b` - Custom bins for FSC/FSR
- `--target-regions` / `-T` - Panel capture regions
- `--mark-input` / `-m` - UXM methylation markers

---

## Gene BED {#gene-bed}

Extended BED format for gene annotations with 4-5 columns.

### Format

```
chrom    start    end    gene    [name]
```

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| chrom | string | ✅ | Chromosome |
| start | int | ✅ | 0-based start |
| end | int | ✅ | 1-based end |
| gene | string | ✅ | Gene symbol (e.g., `TP53`) |
| name | string | Optional | Exon/region name |

### Example

```tsv
chr17	7676594	7676707	TP53	exon1
chr17	7676707	7676863	TP53	exon2
chr7	140719327	140724764	BRAF	exon15
```

### Used By

- Custom gene files for panel FSC

---

## Arms BED {#arms-bed}

Chromosome arm annotations for FSD analysis.

### Format

```
chrom    start    end    arm
```

| Column | Type | Description |
|--------|------|-------------|
| chrom | string | Chromosome (e.g., `chr1`) |
| start | int | 0-based start position |
| end | int | 1-based end position |
| arm | string | Arm identifier (must match pattern: `Np` or `Nq`) |

### Arm Pattern

The arm column **must** match the regex pattern: `^\d{1,2}[pq]$`

Valid examples: `1p`, `1q`, `22p`, `22q`  
Invalid: `Arm1`, `chr1p`, `p`

### Example

```tsv
chr1	0	125000000	1p
chr1	125000000	249250621	1q
chr2	0	93300000	2p
chr2	93300000	243199373	2q
```

### Used By

- `--arms-file` / `-a` - Custom chromosome arms for FSD

---

## WPS Anchors {#wps-anchors}

BED6 format for WPS anchor regions (TSS, CTCF sites, etc.).

### Format

```
chrom    start    end    name    score    strand
```

| Column | Type | Description |
|--------|------|-------------|
| chrom | string | Chromosome |
| start | int | 0-based start |
| end | int | 1-based end |
| name | string | Anchor name (e.g., `TP53_TSS`) |
| score | int | Score (typically 0) |
| strand | string | Strand: `+`, `-`, or `.` |

### Example

```tsv
chr1	11873	14409	DDX11L1_TSS	0	+
chr1	29553	31109	MIR1302-2_TSS	0	+
chr17	7676594	7676707	TP53_TSS	0	-
```

### Used By

- `--wps-anchors` - Custom WPS anchor regions
- `--wps-background` / `-B` - Background normalization regions

---

## Region BED {#region-bed}

Labeled genomic regions for OCF, TFBS, and ATAC analysis.

### Format

```
chrom    start    end    label
```

| Column | Type | Description |
|--------|------|-------------|
| chrom | string | Chromosome |
| start | int | 0-based start |
| end | int | 1-based end |
| label | string | Region label (tissue, TF name, cancer type) |

### Example

```tsv
chr1	100	500	Liver
chr1	600	900	Lung
chr2	1000	1500	Blood
```

### Used By

- `--ocr-file` / `-r` - Open chromatin regions for OCF
- `--tfbs-regions` - Transcription factor binding sites
- `--atac-regions` - ATAC-seq peaks

---

## GC Factors TSV {#gc-factors-tsv}

Tab-separated correction factors for GC bias normalization.

### Format

```
length_bin    gc_pct    factor
```

| Column | Type | Description |
|--------|------|-------------|
| length_bin | int | Fragment length bin: `(length - 60) // 5` |
| gc_pct | int | GC percentage (0-100) |
| factor | float | Correction factor (typically 0.5-2.0) |

### Example

```tsv
length_bin	gc_pct	factor
10	40	1.05
10	41	1.03
10	42	0.98
11	40	1.02
```

### Notes

- Length bin 10 corresponds to fragments 110-114bp
- GC percentage is rounded to the nearest integer
- Factors near 1.0 indicate minimal bias

### Used By

- `--gc-factors` / `-F` - Custom GC correction factors

---

## Compression Support

All BED files can be gzip-compressed (`.bed.gz`). Krewlyzer automatically detects and handles compression.

```bash
# Both work
krewlyzer validate --arms-bed my_arms.bed
krewlyzer validate --arms-bed my_arms.bed.gz
```

---

## Validation

Validate your files before analysis:

```bash
# Validate specific files
krewlyzer validate --gene-bed my_genes.bed
krewlyzer validate --arms-bed my_arms.bed --wps-anchors my_anchors.bed

# Validate bundled assets
krewlyzer validate --genome hg19
```

If validation fails, you'll see:
- Expected format and columns
- Line number of the first error
- Example of correct format

See [Troubleshooting > Asset Validation](../resources/troubleshooting.md#asset-validation) for common errors.
