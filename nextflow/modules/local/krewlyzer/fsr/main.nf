/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_FSR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fragment Size Ratio - short/long fragment ratios for cancer biomarkers.
    Key output: short_long_ratio (tumor fraction proxy)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_FSR {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.5.2"

    input:
    tuple val(meta), path(bed)
    path targets   // Optional

    output:
    tuple val(meta), path("*.FSR.tsv"), emit: tsv
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def pon_arg = params.pon_model ? "--pon-model ${params.pon_model}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""

    """
    krewlyzer fsr \\
        -i $bed \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $targets_arg \\
        $genome_arg \\
        $gc_arg \\
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
    echo -e "region\\tshort_count\\tlong_count\\tshort_long_ratio" > ${prefix}.FSR.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.5.2
    END_VERSIONS
    """
}
