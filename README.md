# Krewlyzer: Comprehensive cfDNA Feature Extraction Toolkit

<p align="center">
  <img src="src/krewlyzer/logo.svg" alt="Krewlyzer logo" width="120"/>
</p>

<p align="center">
  <a href="https://pypi.org/project/krewlyzer/"><img src="https://img.shields.io/pypi/v/krewlyzer.svg?color=blue" alt="PyPI version"></a>
  <a href="https://github.com/msk-access/krewlyzer/actions"><img src="https://github.com/msk-access/krewlyzer/workflows/Release/badge.svg" alt="GitHub Actions"></a>
  <a href="https://github.com/msk-access/krewlyzer/pkgs/container/krewlyzer"><img src="https://img.shields.io/badge/docker-ready-blue.svg" alt="Docker"></a>
</p>

**Krewlyzer** is a high-performance toolkit for extracting biological features from cell-free DNA (cfDNA) sequencing data. Designed for cancer genomics, liquid biopsy research, and clinical bioinformatics.

> [!TIP]
> **Full Documentation**: [msk-access.github.io/krewlyzer](https://msk-access.github.io/krewlyzer/)

---

## Why Krewlyzer?

Cancer cells leave molecular fingerprints in your blood. Krewlyzer finds them.

### The Fragmentomics Advantage

| Traditional Liquid Biopsy | Fragmentomics with Krewlyzer |
|---------------------------|------------------------------|
| Look for specific mutations | Analyze **how DNA is cut** |
| Need prior knowledge of tumor | Works without knowing mutations |
| Miss ~50% of early cancers | Detect more cancers, earlier |

**Key insight**: Tumor DNA fragments are **shorter** (~145bp) than healthy DNA (~166bp). Krewlyzer quantifies this difference and extracts ML-ready features.

### What You Get

| Feature | Clinical Use |
|---------|--------------|
| **Fragment size ratios** | Tumor burden estimation |
| **Cutting patterns** | Tissue of origin identification |
| **Nucleosome positioning** | Epigenetic profiling |
| **Mutation-specific sizes** | MRD monitoring |

> **New to cfDNA?** Read [What is Cell-Free DNA?](https://msk-access.github.io/krewlyzer/introduction/) for background.

---

## Quick Install

```bash
# Docker (recommended)
docker pull ghcr.io/msk-access/krewlyzer:latest

# pip
pip install krewlyzer
```

## Quick Start

```bash
# Run all fragmentomics features
krewlyzer run-all -i sample.bam --reference hg19.fa --output results/

# Individual tools
krewlyzer extract -i sample.bam -r hg19.fa -o output/
krewlyzer fsc -i output/sample.bed.gz -o output/
```

---

## Features

| Command | Description | Output |
|---------|-------------|--------|
| `extract` | Extract fragments from BAM | `.bed.gz` |
| `motif` | End motif & MDS scores | `.EndMotif.tsv`, `.MDS.tsv` |
| `fsc` | Fragment size coverage | `.FSC.tsv` |
| `fsr` | Fragment size ratios | `.FSR.tsv` |
| `fsd` | Size distribution by arm | `.FSD.tsv` |
| `wps` | Windowed protection score | `.WPS.tsv.gz` |
| `ocf` | Orientation-aware fragmentation | `.OCF.tsv` |
| `uxm` | Fragment-level methylation | `.UXM.tsv` |
| `mfsd` | Mutant vs wild-type sizes | `.mFSD.tsv` |
| `run-all` | All features in one pass | All outputs |

---

## Documentation

- [Getting Started](https://msk-access.github.io/krewlyzer/getting-started/) - 5-minute quickstart
- [Installation](https://msk-access.github.io/krewlyzer/installation/) - Docker, pip, development
- [Usage Guide](https://msk-access.github.io/krewlyzer/usage/) - CLI reference
- [Feature Details](https://msk-access.github.io/krewlyzer/features/extract/) - Per-feature documentation
- [Nextflow Pipeline](https://msk-access.github.io/krewlyzer/pipeline/) - Batch processing

---

## Citation

If you use Krewlyzer, please cite:

- **DELFI (FSR):** Cristiano S, et al. Nature 2019
- **WPS:** Snyder MW, et al. Cell 2016
- **OCF:** Sun K, et al. Genome Res 2019
- **UXM:** Loyfer N, et al. Nature 2022

See [Citation & Scientific Background](https://msk-access.github.io/krewlyzer/citation/) for full references.

---

## License

GNU Affero General Public License v3.0 (AGPL-3.0). See [LICENSE](./LICENSE).

---

*Developed by the MSK-ACCESS team at Memorial Sloan Kettering Cancer Center.*
