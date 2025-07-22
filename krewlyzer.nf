#!/usr/bin/env nextflow

/*
 * Krewlyzer Nextflow pipeline: run-all features for cfDNA BAMs
 * Supports LSF and SLURM via Nextflow executor profiles
 * Usage:
 *   nextflow run krewlyzer.nf --samplesheet samplesheet.csv --ref /path/to/reference.fa
 *
 * For LSF:         nextflow run krewlyzer.nf -profile lsf ...
 * For SLURM:       nextflow run krewlyzer.nf -profile slurm ...
 * For Docker:      nextflow run krewlyzer.nf -profile docker ...
 * For Singularity: nextflow run krewlyzer.nf -profile singularity ...
 */

params.samplesheet = ''
params.ref = ''
params.krewlyzer = 'krewlyzer' // or path to CLI

// Parse the sample sheet CSV into a channel of [sample_id, bam]
Channel.fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample_id, file(row.bam)) }
    .set { SAMPLES }

process RUN_ALL {
    container 'ghcr.io/msk-access/krewlyzer:latest'
    tag "${sample_id}"
    cpus 8
    memory '32 GB'
    time '24h'

    input:
    val sample_id
    path bam_file

    script:
    """
    OUTDIR="$sample_id"
    mkdir -p "$OUTDIR"
    ${params.krewlyzer} run-all ${bam_file} --reference ${params.ref} --output "$OUTDIR" --threads 8
    """
    output:
    path "*", emit: out
}

workflow {
    RUN_ALL(SAMPLES)
}
