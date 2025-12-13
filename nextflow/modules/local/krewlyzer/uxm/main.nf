process KREWLYZER_UXM {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.2.0"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.UXM.tsv"), emit: tsv
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    krewlyzer uxm \\
        $bam \\
        --reference $fasta \\
        --output ./ \\
        --sample-name $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
