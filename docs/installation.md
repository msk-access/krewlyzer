# Installation

## Requirements

- **OS**: Linux or macOS (tested on Ubuntu 20.04+, macOS 12+)
- **Python**: 3.10+
- **RAM**: ≥16GB recommended for large BAM files
- **Reference**: Indexed FASTA file (hg19 or hg38)

---

## Docker (Recommended)

The easiest way to run Krewlyzer with all dependencies:

```bash
docker pull ghcr.io/msk-access/krewlyzer:latest
```

### Running with Docker

```bash
docker run --rm -v $PWD:/data ghcr.io/msk-access/krewlyzer:latest \
    run-all -i /data/sample.bam \
    --reference /data/hg19.fa \
    --output /data/results/
```

!!! tip "Volume Mounting"
    Use `-v $PWD:/data` to mount your current directory. All paths in the command should then use `/data/` prefix.

---

## pip / uv

### Standard Installation

```bash
pip install krewlyzer
```

Or using [uv](https://github.com/astral-sh/uv) for faster installs:

```bash
uv pip install krewlyzer
```

### Verify Installation

```bash
krewlyzer --version
krewlyzer --help
```

---

## Development Setup

For contributing or development:

```bash
# Clone repository
git clone https://github.com/msk-access/krewlyzer.git
cd krewlyzer

# Create virtual environment
uv venv .venv
source .venv/bin/activate

# Install in development mode with Rust compilation
uv pip install -e ".[dev,test]"

# Verify Rust core is built
python -c "from krewlyzer import _core; print('Rust core loaded')"
```

### Rust Development

The Rust backend requires:

- Rust toolchain (install via [rustup](https://rustup.rs/))
- C compiler (clang recommended)
- htslib dependencies

```bash
# Build Rust extension
cd rust
maturin develop --release
```

---

## Reference Data

### Reference Genome

Download and index the reference genome:

```bash
# hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa

# hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

### Bundled Data Files

Krewlyzer includes default annotation files in `src/krewlyzer/data/`:

| Directory | Contents |
|-----------|----------|
| `ChromosomeBins/` | 100kb genome bins (hg19, hg38) |
| `ChromosomeArms/` | Chromosome arm definitions |
| `WpsAnchors/` | WPS anchor regions (TSS + CTCF) |
| `OpenChromatinRegion/` | Tissue-specific OCR for OCF |
| `MethMark/` | Methylation markers for UXM |
| `pon/` | Panel of Normals models |

---

## Troubleshooting

### "ModuleNotFoundError: krewlyzer._core"

The Rust extension failed to build. Ensure you have:

- Python 3.10+
- C compiler (`gcc` or `clang`)
- Rust toolchain

Reinstall with verbose output:

```bash
pip install krewlyzer -v
```

### "htslib not found"

Install htslib development files:

```bash
# Ubuntu/Debian
sudo apt-get install libhts-dev

# macOS
brew install htslib
```

### Memory Errors

For large BAM files, increase available memory or process chromosomes separately:

```bash
krewlyzer extract sample.bam -r hg19.fa -o output/ --chromosomes chr1,chr2
```
