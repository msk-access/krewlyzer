/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_MOTIF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    End Motif and Motif Diversity Score (MDS) extraction.
    Captures cfDNA cleavage preferences for enzymatic signature analysis.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_MOTIF {
    tag "$meta.id"
    label 'process_high'
    container "ghcr.io/msk-access/krewlyzer:0.9.2"

    // Into the sample's folder, as run-all publishes. `{sample}.metadata.*` is never among a
    // tool's outputs: it is downstream's completion marker.
    publishDir path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    input:
    // The BAM: motifs are read from it, as run-all reads them during extraction. The sample's PON
    // and targets, resolved once by INPUT_CHECK as run-all's are ([] if none).
    tuple val(meta), path(bam), path(bai), path(pon), path(targets)
    path fasta

    output:
    // In params.output_format; on-target tables only in panel mode.
    tuple val(meta), path("*.EndMotif.{parquet,tsv,tsv.gz}")                , emit: end_motif
    tuple val(meta), path("*.EndMotif1mer.{parquet,tsv,tsv.gz}")            , emit: end_motif_1mer
    tuple val(meta), path("*.BreakPointMotif.{parquet,tsv,tsv.gz}")         , emit: breakpoint_motif
    tuple val(meta), path("*.MDS.{parquet,tsv,tsv.gz}")                     , emit: mds
    tuple val(meta), path("*.EndMotif.ontarget.{parquet,tsv,tsv.gz}")       , emit: end_motif_ontarget, optional: true
    tuple val(meta), path("*.BreakPointMotif.ontarget.{parquet,tsv,tsv.gz}"), emit: breakpoint_motif_ontarget, optional: true
    tuple val(meta), path("*.MDS.ontarget.{parquet,tsv,tsv.gz}")            , emit: mds_ontarget, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""
    def pon_arg = pon ? "--pon-model ${pon}" : ""
    def skip_pon_arg = params.skip_pon ? "--skip-pon" : ""
    def pon_variant_arg = params.pon_variant != 'all_unique' ? "--pon-variant ${params.pon_variant}" : ""
    // run-all's read filters; the tool applies the genome's bundled blacklist, as run-all does.
    def mapq_arg = params.mapq != 20 ? "--mapq ${params.mapq}" : ""
    def minlen_arg = params.minlen != 65 ? "--minlen ${params.minlen}" : ""
    def maxlen_arg = params.maxlen != 1000 ? "--maxlen ${params.maxlen}" : ""
    def skip_dup_arg = params.skip_duplicates == false ? "--no-skip-duplicates" : ""
    def proper_pair_arg = params.require_proper_pair == false ? "--no-require-proper-pair" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def output_format_arg = params.output_format && params.output_format != 'tsv' ? "--output-format ${params.output_format}" : ""
    def compress_arg = params.compress_tsv ? "--compress" : ""

    """
    krewlyzer motif \\
        -i $bam \\
        -r $fasta \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $assay_arg \\
        $genome_arg \\
        $targets_arg \\
        $skip_targets_arg \\
        $pon_arg \\
        $skip_pon_arg \\
        $pon_variant_arg \\
        $mapq_arg \\
        $minlen_arg \\
        $maxlen_arg \\
        $skip_dup_arg \\
        $proper_pair_arg \\
        $verbose_arg \\
        $output_format_arg \\
        $compress_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "motif\\tcount\\tfrequency" > ${prefix}.EndMotif.tsv
    echo -e "sample\\tmds" > ${prefix}.MDS.tsv
    touch ${prefix}.EndMotif1mer.tsv
    touch ${prefix}.BreakPointMotif.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.9.2
    END_VERSIONS
    """
}
