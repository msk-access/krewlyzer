# Contributing Data Files

This guide explains how to add or modify data files in the Krewlyzer repository.

## Git LFS Setup

The `src/krewlyzer/data/` folder uses **Git LFS** (Large File Storage) for binary files.

### File Types Tracked

```gitattributes
src/krewlyzer/data/**/*.gz filter=lfs diff=lfs merge=lfs -text
src/krewlyzer/data/**/*.parquet filter=lfs diff=lfs merge=lfs -text
src/krewlyzer/data/**/*.bed filter=lfs diff=lfs merge=lfs -text
```

### Prerequisites

1. Install Git LFS:
   ```bash
   # macOS
   brew install git-lfs
   
   # Ubuntu
   apt-get install git-lfs
   ```

2. Initialize LFS in your clone:
   ```bash
   cd krewlyzer
   git lfs install
   ```

### Cloning with LFS

```bash
# Clone with LFS files
git clone https://github.com/msk-access/krewlyzer.git

# If LFS files weren't pulled, fetch them
git lfs pull
```

---

## Adding New Panel Assets

### 1. Target Regions

Place in `src/krewlyzer/data/targets/GRCh37/`:

```bash
# File: {assay}.targets.bed.gz
# Format: BED3 or BED4
chr1  11166102  11166202  MTOR_exon1
chr1  27022522  27022622  ARID1A_exon1
```

### 2. TFBS Regions (Pre-intersected)

Place in `src/krewlyzer/data/TFBS/GRCh37/`:

```bash
# Create panel-specific TFBS by intersecting with targets
bedtools intersect -u \
    -a Homo_sapiens_meta_clusters_hg19.bed.gz \
    -b ../targets/GRCh37/{assay}.targets.bed.gz \
    | bgzip > {assay}.meta_clusters_hg19.bed.gz
tabix -p bed {assay}.meta_clusters_hg19.bed.gz
```

### 3. ATAC Regions (Pre-intersected)

Place in `src/krewlyzer/data/ATAC/GRCh37/`:

```bash
# Create panel-specific ATAC by intersecting with targets
bedtools intersect -u \
    -a TCGA_ATAC_peak.hg19.bed.gz \
    -b ../targets/GRCh37/{assay}.targets.bed.gz \
    | bgzip > {assay}.TCGA_ATAC_peak.hg19.bed.gz
tabix -p bed {assay}.TCGA_ATAC_peak.hg19.bed.gz
```

### 4. WPS Anchors

Place in `src/krewlyzer/data/WpsAnchors/GRCh37/`:

```bash
# File: {assay}.wps_anchors.bed.gz
# Format: BED6 with TSS/CTCF annotations
chr1  10000  10500  TSS  100  +
chr1  15000  15500  CTCF  100  -
```

---

## File Format Requirements

| File Type | Format | Index | Notes |
|-----------|--------|-------|-------|
| Targets | BED3/4 | Optional | 0-based coordinates |
| TFBS | BED4 | `.tbi` | Col4 = TF name |
| ATAC | BED4 | `.tbi` | Col4 = cancer type |
| WPS Anchors | BED6 | `.tbi` | TSS/CTCF regions |
| Gene BEDs | BED6 | None | Per-gene exon regions |

### Compression

- Use **bgzip** (not gzip) for indexed files:
  ```bash
  bgzip file.bed
  tabix -p bed file.bed.gz
  ```

- BGZF format is required for tabix indexing

---

## Registering Assets

After adding files, update `src/krewlyzer/assets.py`:

```python
# In AssetManager.get_targets()
if assay == "new_assay":
    return self.targets_dir / "new_assay.targets.bed.gz"
```

---

## Committing LFS Files

```bash
# Stage new files
git add src/krewlyzer/data/targets/GRCh37/new_assay.targets.bed.gz

# Verify LFS tracking
git lfs status

# Commit
git commit -m "feat: add new_assay target regions"

# Push (uploads to LFS)
git push origin main
```

---

## Troubleshooting

### "LFS object not found"

```bash
git lfs fetch --all
git lfs checkout
```

### Large clone size

```bash
# Clone without LFS content (faster)
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/msk-access/krewlyzer.git

# Pull only needed LFS files
git lfs pull --include="src/krewlyzer/data/gc/**"
```
