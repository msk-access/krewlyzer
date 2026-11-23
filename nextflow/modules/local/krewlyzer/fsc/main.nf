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
    container "ghcr.io/msk-access/krewlyzer:0.9.2"

    // Into the sample's folder, as run-all publishes. `{sample}.metadata.*` is never among a
    // tool's outputs: it is downstream's completion marker.
    publishDir path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    input:
    // The sample's PON and targets, resolved once by INPUT_CHECK as run-all's are ([] if none).
    // gc_factors: extract's tables, staged beside the BED ([] if none).
    tuple val(meta), path(bed), path(gc_factors), path(pon), path(targets)

    output:
    // In params.output_format; on-target and panel tables only in panel mode.
    tuple val(meta), path("*.FSC.{parquet,tsv,tsv.gz}")                , emit: table
    tuple val(meta), path("*.FSC.ontarget.{parquet,tsv,tsv.gz}")       , emit: ontarget, optional: true
    tuple val(meta), path("*.FSC.gene.{parquet,tsv,tsv.gz}")           , emit: gene_fsc, optional: true
    tuple val(meta), path("*.FSC.regions.{parquet,tsv,tsv.gz}")        , emit: regions, optional: true
    tuple val(meta), path("*.FSC.regions.e1only.{parquet,tsv,tsv.gz}") , emit: regions_e1only, optional: true
    tuple val(meta), path("*.fsc_counts.{parquet,tsv,tsv.gz}")         , emit: counts
    tuple val(meta), path("*.fsc_counts.ontarget.{parquet,tsv,tsv.gz}"), emit: counts_ontarget, optional: true
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def pon_arg = pon ? "--pon-model ${pon}" : ""
    def skip_pon_arg = params.skip_pon ? "--skip-pon" : ""
    def pon_variant_arg = params.pon_variant != 'all_unique' ? "--pon-variant ${params.pon_variant}" : ""
    def e1_disable_arg = params.disable_e1_aggregation ? "--disable-e1-aggregation" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def output_format_arg = params.output_format && params.output_format != 'tsv' ? "--output-format ${params.output_format}" : ""
    def compress_arg = params.compress_tsv ? "--compress" : ""
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
        $skip_pon_arg \\
        $pon_variant_arg \\
        $e1_disable_arg \\
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
    echo -e "chrom\\tstart\\tend\\tultra_short\\tcore_short\\tmono_nucl\\tdi_nucl\\tlong\\ttotal" > ${prefix}.FSC.tsv
    touch ${prefix}.fsc_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.9.2
    END_VERSIONS
    """
}
