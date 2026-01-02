process KREWLYZER_WPS {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:latest"

    input:
    tuple val(meta), path(bed)
    path fasta
    path wps_anchors
    path wps_background

    output:
    tuple val(meta), path("*.WPS.parquet")           , emit: parquet
    tuple val(meta), path("*.WPS_background.parquet"), emit: background, optional: true
    path "versions.yml"                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def ref_arg = fasta ? "--reference ${fasta}" : ""
    def anchors_arg = wps_anchors ? "--wps-anchors ${wps_anchors}" : ""
    def background_arg = wps_background ? "--background ${wps_background}" : ""

    """
    krewlyzer wps \\
        $bed \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $ref_arg \\
        $genome_arg \\
        $gc_arg \\
        $anchors_arg \\
        $background_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
