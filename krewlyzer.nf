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

// Parse the sample sheet CSV into a channel of [sample_id, bam, variant_file]
Channel.fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> 
        def bam = file(row.bam)
        def variant = row.vcf ? file(row.vcf) : (row.maf ? file(row.maf) : [])
        tuple(row.sample_id, bam, variant)
    }
    .set { SAMPLES }

process RUN_ALL {
    container 'ghcr.io/msk-access/krewlyzer:latest'
    tag "${sample_id}"
    cpus 8
    memory '32 GB'
    time '24h'

    input:
    tuple val(sample_id), path(bam_file), path(variant_file)

    script:
    def variant_arg = variant_file ? "--variant-input ${variant_file}" : ""
    """
    OUTDIR="$sample_id"
    mkdir -p "$OUTDIR"
    ${params.krewlyzer} run-all ${bam_file} --reference ${params.ref} --output "$OUTDIR" --threads 8 ${variant_arg}
    """
    output:
    path "*", emit: out
}

workflow {
    RUN_ALL(SAMPLES)
}
