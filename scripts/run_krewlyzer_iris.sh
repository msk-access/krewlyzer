#!/bin/bash
#SBATCH --job-name=krewlyzer
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=krewlyzer_%j.log
#SBATCH --error=krewlyzer_%j.err

# ============================================================================
# Krewlyzer Nextflow Pipeline - IRIS SLURM Submission Script
# ============================================================================
# Usage:  sbatch run_krewlyzer_iris.sh
# Resume: sbatch run_krewlyzer_iris.sh --resume
#
# Resource math (14K+ samples):
#   - Head process: 4 CPUs + 32GB (Nextflow JVM tracking all tasks)
#   - Per-sample jobs: 8 CPUs + 32GB (submitted by Nextflow to SLURM)
#   - queue_size=200: up to 200 concurrent sample jobs
#   - Estimated: 14350 samples / 200 concurrent × ~1.5h each ≈ 3-4 days
# ============================================================================

set -euo pipefail

# Activate Nextflow environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate nf-env

# Optional: pass -resume if provided as argument
RESUME_FLAG=""
if [[ "${1:-}" == "--resume" ]]; then
    RESUME_FLAG="-resume"
    echo ">>> Resuming previous run..."
fi

echo ">>> Starting Krewlyzer pipeline at $(date)"
echo ">>> Working directory: $PWD"

nextflow run /usersoftware/shahr2/github/krewlyzer/nextflow/main.nf \
    -profile iris \
    --partition cmobic_cpu \
    --isolated true \
    --samplesheet $PWD/samplesheet.csv \
    --ref /data1/core006/access/production/resources/reference/versions/hg19/Homo_sapiens_assembly19.fasta \
    --outdir $PWD/results/ \
    --output_format both \
    --compress_tsv true \
    --generate_json true \
    --verbose true \
    --queue_size 200 \
    -qs 200 \
    ${RESUME_FLAG}

echo ">>> Pipeline completed at $(date)"
