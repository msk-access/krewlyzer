# Installation

## Quick Reference

| Method | Best For | Includes Data |
|--------|----------|---------------|
| **Docker** | Production & HPC | ✅ All bundled |
| **Clone + Install** | Development | ✅ Via Git LFS |
| **pip + Data Clone** | Custom environments | ⚠️ Requires env var |

---

## Option 1: Docker (Recommended)

The easiest way to run Krewlyzer with all dependencies and data:

```bash
docker pull ghcr.io/msk-access/krewlyzer:latest
```

### Running with Docker

```bash
docker run --rm -v $PWD:/data ghcr.io/msk-access/krewlyzer:latest \
    run-all -i /data/sample.bam \
    --reference /data/hg19.fa \
    --output /data/results/ \
    --assay xs2
```

!!! tip "Volume Mounting (Standalone Docker)"
    Use `-v $PWD:/data` to mount your current directory. All paths in the command should use the `/data/` prefix. **For Nextflow pipelines**, volume mounting is automatic—just use host paths in your samplesheet.

---

## Option 2: Clone Repository

Full installation with bundled data (for development or when Docker isn't available):

```bash
# Clone repository with LFS data
git clone https://github.com/msk-access/krewlyzer.git
cd krewlyzer
git lfs pull

# Install in development mode
pip install -e .

# Verify
krewlyzer --version
```

!!! info "Why This Works Without Configuration"
    With `pip install -e .` (editable mode), Python runs code directly from the source directory.
    Asset paths resolve to `src/krewlyzer/data/` where all LFS files exist.

---

## Option 3: pip Install + Data Clone

For environments where you want PyPI code with external data:

### Step 1: Install Package

```bash
pip install krewlyzer
```

### Step 2: Clone Data Repository

```bash
# Shallow clone (faster, code not needed)
git clone --depth 1 https://github.com/msk-access/krewlyzer.git ~/.krewlyzer-data
cd ~/.krewlyzer-data && git lfs pull
```

### Step 3: Configure Environment Variable

```bash
# Set for current session
export KREWLYZER_DATA_DIR=~/.krewlyzer-data/src/krewlyzer/data

# Add to shell profile for persistence
echo 'export KREWLYZER_DATA_DIR=~/.krewlyzer-data/src/krewlyzer/data' >> ~/.bashrc
```

!!! warning "Required for pip Install"
    The `KREWLYZER_DATA_DIR` environment variable is **required** when installing via pip.
    Without it, asset auto-loading will fail. You can still use explicit paths like `--pon-model`.

---

## Requirements

- **OS**: Linux or macOS (tested on Ubuntu 20.04+, macOS 12+)
- **Python**: 3.10+
- **RAM**: ≥16GB recommended for large BAM files
- **Reference**: Indexed FASTA file (hg19 or hg38)

---

## Reference Genome Setup

Download and index the reference genome:

=== "hg19 (GRCh37)"
    ```bash
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip hg19.fa.gz
    samtools faidx hg19.fa
    ```

=== "hg38 (GRCh38)"
    ```bash
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip hg38.fa.gz
    samtools faidx hg38.fa
    ```

---

## Bundled Data Files

Krewlyzer includes annotation files in `src/krewlyzer/data/`:

| Directory | Contents |
|-----------|----------|
| `ChromosomeBins/` | 100kb genome bins (hg19, hg38) |
| `ChromosomeArms/` | Chromosome arm definitions |
| `WpsAnchors/` | WPS anchor regions (TSS + CTCF) |
| `OpenChromatinRegion/` | Tissue-specific OCR for OCF |
| `MethMark/` | Methylation markers for UXM |
| `pon/` | Panel of Normals models |
| `TFBS/` | GTRD meta-clusters for region entropy |
| `ATAC/` | TCGA ATAC peaks for region entropy |

---

## Troubleshooting

### "Asset not found" or "PON not found"

If you installed via `pip install krewlyzer`, you need to set up the data directory:

```bash
git clone --depth 1 https://github.com/msk-access/krewlyzer.git ~/.krewlyzer-data
cd ~/.krewlyzer-data && git lfs pull
export KREWLYZER_DATA_DIR=~/.krewlyzer-data/src/krewlyzer/data
```

### "ModuleNotFoundError: krewlyzer._core"

The Rust extension failed to build. Ensure you have:

- Python 3.10+
- C compiler (`gcc` or `clang`)
- Rust toolchain (install via [rustup](https://rustup.rs/))

Reinstall with verbose output:

```bash
pip install krewlyzer -v
```

### "htslib not found"

Install htslib development files:

=== "Ubuntu/Debian"
    ```bash
    sudo apt-get install libhts-dev
    ```

=== "macOS"
    ```bash
    brew install htslib
    ```

### Memory Errors

For large BAM files, increase available memory or process chromosomes separately:

```bash
krewlyzer extract sample.bam -r hg19.fa -o output/ --chromosomes chr1,chr2
```
