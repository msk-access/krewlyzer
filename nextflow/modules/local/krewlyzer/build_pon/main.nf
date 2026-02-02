/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_BUILD_PON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Build Panel of Normals (PON) model from healthy plasma samples.
    Creates unified PON model with GC bias, FSD, WPS, and MDS baselines.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_BUILD_PON {
    tag "$assay"
    label 'process_high'
    container "ghcr.io/msk-access/krewlyzer:0.5.0"

    input:
    path sample_list          // Text file with BAM/CRAM/BED.gz paths (one per line)
    val assay                 // Assay name (e.g., msk-access-v2, xs1, xs2)
    path fasta                // Reference FASTA
    path targets              // Optional target regions BED for panel mode
    path wps_anchors          // Optional WPS anchors BED

    output:
    path "*.pon.parquet"      , emit: pon
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def anchors_arg = wps_anchors ? "--wps-anchors ${wps_anchors}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""

    """
    krewlyzer build-pon \\
        ${sample_list} \\
        --assay ${assay} \\
        --reference ${fasta} \\
        --output ${assay}.pon.parquet \\
        --threads $task.cpus \\
        $genome_arg \\
        $targets_arg \\
        $anchors_arg \\
        $verbose_arg \\
        $skip_targets_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer version: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${assay}.pon.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.5.0
    END_VERSIONS
    """
}
