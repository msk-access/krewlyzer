process KREWLYZER_FSC {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.3.2"

    input:
    tuple val(meta), path(bed)
    path targets   // Optional for custom bins

    output:
    tuple val(meta), path("*.FSC.tsv"), emit: tsv
    tuple val(meta), path("*.FSC.gene.tsv"), emit: gene_fsc, optional: true
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def pon_arg = params.pon_model ? "--pon-model ${params.pon_model}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""

    """
    krewlyzer fsc \\
        -i $bed \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $targets_arg \\
        $assay_arg \\
        $genome_arg \\
        $gc_arg \\
        $pon_arg \\
        $verbose_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
