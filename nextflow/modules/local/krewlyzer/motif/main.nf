/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_MOTIF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    End Motif and Motif Diversity Score (MDS) extraction.
    Captures cfDNA cleavage preferences for enzymatic signature analysis.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_MOTIF {
    tag "$meta.id"
    label 'process_high'
    container "ghcr.io/msk-access/krewlyzer:0.3.2"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.EndMotif.tsv"), emit: end_motif
    tuple val(meta), path("*.MDS.tsv")     , emit: mds
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def verbose_arg = params.verbose ? "--verbose" : ""

    """
    krewlyzer motif \\
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
    echo -e "motif\\tcount\\tfrequency" > ${prefix}.EndMotif.tsv
    echo -e "sample\\tmds" > ${prefix}.MDS.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.3.2
    END_VERSIONS
    """
}
