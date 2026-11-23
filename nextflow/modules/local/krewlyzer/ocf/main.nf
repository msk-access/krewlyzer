/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_OCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Orientation-aware cfDNA Fragmentation - strand asymmetry at Open Chromatin Regions.
    Used for tissue-of-origin deconvolution (Sun et al. 2019).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_OCF {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.9.2"

    // Into the sample's folder, as run-all publishes. `{sample}.metadata.*` is never among a
    // tool's outputs: it is downstream's completion marker.
    publishDir path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    input:
    // The sample's PON and targets, resolved once by INPUT_CHECK as run-all's are ([] if none).
    // gc_factors: extract's tables, staged beside the BED ([] if none).
    tuple val(meta), path(bed), path(gc_factors), path(pon), path(targets)
    path fasta

    output:
    // In params.output_format. On- and off-target tables in panel mode only, and each only when
    // a fragment of that kind falls in an OCR (the synthetic test sample has no on-target one).
    tuple val(meta), path("*.OCF.{parquet,tsv,tsv.gz}")               , emit: table
    tuple val(meta), path("*.OCF.sync.{parquet,tsv,tsv.gz}")          , emit: sync
    tuple val(meta), path("*.OCF.ontarget.{parquet,tsv,tsv.gz}")      , emit: ontarget, optional: true
    tuple val(meta), path("*.OCF.ontarget.sync.{parquet,tsv,tsv.gz}") , emit: ontarget_sync, optional: true
    tuple val(meta), path("*.OCF.offtarget.{parquet,tsv,tsv.gz}")     , emit: offtarget, optional: true
    tuple val(meta), path("*.OCF.offtarget.sync.{parquet,tsv,tsv.gz}"), emit: offtarget_sync, optional: true
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def pon_arg = pon ? "--pon-model ${pon}" : ""
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    def skip_pon_arg = params.skip_pon ? "--skip-pon" : ""
    def pon_variant_arg = params.pon_variant != 'all_unique' ? "--pon-variant ${params.pon_variant}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def output_format_arg = params.output_format && params.output_format != 'tsv' ? "--output-format ${params.output_format}" : ""
    def compress_arg = params.compress_tsv ? "--compress" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""

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
        $assay_arg \\
        $skip_pon_arg \\
        $pon_variant_arg \\
        $verbose_arg \\
        $output_format_arg \\
        $compress_arg \\
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
    echo -e "tissue\\tocf_score\\tstrand_bias" > ${prefix}.OCF.tsv
    touch ${prefix}.OCF.sync.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.9.2
    END_VERSIONS
    """
}
