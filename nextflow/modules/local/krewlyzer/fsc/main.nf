process KREWLYZER_FSC {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.2.0"

    input:
    tuple val(meta), path(bed)
    path targets   // Optional for custom bins

    output:
    tuple val(meta), path("*.FSC.tsv"), emit: tsv
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def targets_arg = targets ? "--bin-input ${targets}" : ""

    """
    krewlyzer fsc \\
        $bed \\
        --output ./ \\
        --sample-name $prefix \\
        $targets_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
