process KREWLYZER_OCF {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.3.2"

    input:
    tuple val(meta), path(bed)
    path fasta

    output:
    tuple val(meta), path("*.OCF.tsv"), emit: tsv
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""

    """
    krewlyzer ocf \\
        $bed \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $genome_arg \\
        $gc_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
