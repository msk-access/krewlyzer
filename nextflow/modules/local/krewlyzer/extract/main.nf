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
    container "ghcr.io/msk-access/krewlyzer:0.3.2"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    tuple val(meta), path("*.json")  , emit: metadata
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def verbose_arg = params.verbose ? "--verbose" : ""

    """
    krewlyzer extract \\
        -i $bam \\
        -r $fasta \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
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
    echo '{"sample_id":"${prefix}","total_fragments":1000}' > ${prefix}.metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.3.2
    END_VERSIONS
    """
}
