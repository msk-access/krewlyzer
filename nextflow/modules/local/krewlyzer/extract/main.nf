/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_EXTRACT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Extract cfDNA fragments from BAM to BED.gz with GC correction factors.
    First step in krewlyzer pipeline - generates input for all feature tools.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_EXTRACT {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.9.1"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path targets            // Optional target regions for panel mode GC correction

    output:
    tuple val(meta), path("*.bed.gz")     , emit: bed
    // Metadata moved from .metadata.json to .metadata.tsv/.tsv.gz/.parquet in
    // 0.7.0. This still globbed "*.json", which the tool has not written
    // since, so the process failed on a missing required output and the
    // tool-level path could not run at all. Optional because the extension
    // now depends on --output-format/--compress.
    tuple val(meta), path("*.metadata.tsv")    , emit: metadata,         optional: true
    tuple val(meta), path("*.metadata.tsv.gz") , emit: metadata_gz,      optional: true
    tuple val(meta), path("*.metadata.parquet"), emit: metadata_parquet, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def assay_arg = params.assay ? "--assay ${params.assay}" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""
    def maxlen_arg = params.maxlen != 1000 ? "--maxlen ${params.maxlen}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def output_format_arg = params.output_format && params.output_format != 'tsv' ? "--output-format ${params.output_format}" : ""
    def compress_arg = params.compress_tsv ? "--compress" : ""

    """
    krewlyzer extract \\
        -i $bam \\
        -r $fasta \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $genome_arg \\
        $assay_arg \\
        $targets_arg \\
        $skip_targets_arg \\
        $maxlen_arg \\
        $verbose_arg \\
        $output_format_arg \\
        $compress_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi
    printf 'sample_id\\ttotal_fragments\\n${prefix}\\t1000\\n' > ${prefix}.metadata.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.8.3
    END_VERSIONS
    """
}
