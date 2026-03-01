/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_WPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Windowed Protection Score - nucleosome/TF profiling at TSS/CTCF anchors.
    Dual-stream: WPS-Nuc (150-180bp) and WPS-TF (35-80bp) fragments.
    Includes FFT periodicity extraction for Nucleosome Repeat Length (NRL).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_WPS {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.6.0"

    input:
    tuple val(meta), path(bed)
    path fasta
    path wps_anchors
    path wps_background

    output:
    tuple val(meta), path("*.WPS.parquet")           , emit: parquet
    tuple val(meta), path("*.WPS_background.parquet"), emit: background, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def ref_arg = fasta ? "-r ${fasta}" : ""
    def anchors_arg = wps_anchors ? "--wps-anchors ${wps_anchors}" : ""
    def background_arg = wps_background ? "--background ${wps_background}" : ""
    def targets_arg = params.targets ? "--target-regions ${params.targets}" : ""
    def bait_arg = params.bait_padding ? "--bait-padding ${params.bait_padding}" : ""
    def pon_arg = params.pon_model ? "--pon-model ${params.pon_model}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""

    """
    krewlyzer wps \\
        -i $bed \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $ref_arg \\
        $genome_arg \\
        $gc_arg \\
        $anchors_arg \\
        $background_arg \\
        $targets_arg \\
        $bait_arg \\
        $pon_arg \\
        $verbose_arg \\
        $skip_targets_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.WPS.parquet
    touch ${prefix}.WPS_background.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.6.0
    END_VERSIONS
    """
}
