/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_FSC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fragment Size Coverage - windowed depth per size channel.
    Channels: ultra_short, core_short, mono_nucl, di_nucl, long
    Panel mode: Also outputs gene-centric FSC with --assay flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_FSC {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.5.2"

    input:
    tuple val(meta), path(bed)
    path targets   // Optional for custom bins

    output:
    tuple val(meta), path("*.FSC.tsv"), emit: tsv
    tuple val(meta), path("*.FSC.gene.tsv"), emit: gene_fsc, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def pon_arg = params.pon_model ? "--pon-model ${params.pon_model}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""

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
    echo -e "chrom\\tstart\\tend\\tultra_short\\tcore_short\\tmono_nucl\\tdi_nucl\\tlong\\ttotal" > ${prefix}.FSC.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.5.2
    END_VERSIONS
    """
}
