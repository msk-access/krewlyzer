/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_REGION_ENTROPY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Region Entropy - Shannon entropy of fragment sizes at regulatory regions.
    Features: TFBS (808 TFs) and ATAC (23 cancer types) size entropy.
    Used for tissue-of-origin and epigenetic signal detection.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_REGION_ENTROPY {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.5.2"

    input:
    tuple val(meta), path(bed)
    path fasta
    path pon_model          // Optional PON model for z-score normalization
    path targets            // Optional target regions for panel mode

    output:
    tuple val(meta), path("*.TFBS.tsv")           , emit: tfbs, optional: true
    tuple val(meta), path("*.ATAC.tsv")           , emit: atac, optional: true
    tuple val(meta), path("*.TFBS.ontarget.tsv")  , emit: tfbs_ontarget, optional: true
    tuple val(meta), path("*.ATAC.ontarget.tsv")  , emit: atac_ontarget, optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def pon_arg = pon_model ? "--pon-model ${pon_model}" : ""
    def skip_pon_arg = params.skip_pon ? "--skip-pon" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def tfbs_arg = params.no_tfbs ? "--no-tfbs" : ""
    def atac_arg = params.no_atac ? "--no-atac" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""

    """
    krewlyzer region-entropy \\
        -i $bed \\
        --output ./ \\
        --sample-name $prefix \\
        $genome_arg \\
        $pon_arg \\
        $skip_pon_arg \\
        $targets_arg \\
        $tfbs_arg \\
        $atac_arg \\
        $verbose_arg \\
        $skip_targets_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "label\\tentropy\\tz_score" > ${prefix}.TFBS.tsv
    echo -e "label\\tentropy\\tz_score" > ${prefix}.ATAC.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.5.2
    END_VERSIONS
    """
}
