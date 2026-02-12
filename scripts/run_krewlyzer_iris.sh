#!/bin/bash
#SBATCH --job-name=krewlyzer
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=krewlyzer_%j.log
#SBATCH --error=krewlyzer_%j.err

# ============================================================================
# Krewlyzer Nextflow Pipeline - IRIS SLURM Submission Script
# ============================================================================
# Usage:  sbatch run_krewlyzer_iris.sh
# Resume: sbatch run_krewlyzer_iris.sh --resume
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
    -qs 100 \
    ${RESUME_FLAG}

echo ">>> Pipeline completed at $(date)"
