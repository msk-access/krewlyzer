process KREWLYZER_MOTIF {
    tag "$meta.id"
    label 'process_high'
    container "ghcr.io/msk-access/krewlyzer:0.2.0"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.EndMotif.tsv"), emit: end_motif
    tuple val(meta), path("*.MDS.tsv")     , emit: mds
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    krewlyzer motif \\
        $bam \\
        --reference $fasta \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
