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
    def targets_arg = params.targets ? "--target-regions ${params.targets}" : ""
    def pon_arg = params.pon_model ? "--pon-model ${params.pon_model}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""

    """
    krewlyzer ocf \\
        -i $bed \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $genome_arg \\
        $gc_arg \\
        $targets_arg \\
        $pon_arg \\
        $verbose_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
