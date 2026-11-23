/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_EXTRACT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Extract cfDNA fragments from BAM to BED.gz with GC correction factors.
    First step in krewlyzer pipeline - generates input for all feature tools.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_EXTRACT {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.9.2"

    // Into the sample's folder, as run-all publishes, but not the record. `{sample}.metadata.*`
    // is downstream's completion marker, and extraction's record has neither
    // `krewlyzer_version` nor `pon_applied` and would appear before any feature table.
    publishDir path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' || filename.contains('.metadata.') ? null : filename }

    input:
    tuple val(meta), path(bam), path(bai), path(targets)  // targets: the sample's ([] if none)
    path fasta

    output:
    tuple val(meta), path("*.bed.gz")     , emit: bed
    tuple val(meta), path("*.bed.gz.tbi") , emit: bed_index,        optional: true
    // Metadata moved from .metadata.json to .metadata.tsv/.tsv.gz/.parquet in
    // 0.7.0. This still globbed "*.json", which the tool has not written
    // since, so the process failed on a missing required output and the
    // tool-level path could not run at all. Optional because the extension
    // now depends on --output-format/--compress.
    tuple val(meta), path("*.metadata.tsv")    , emit: metadata,         optional: true
    tuple val(meta), path("*.metadata.tsv.gz") , emit: metadata_gz,      optional: true
    tuple val(meta), path("*.metadata.parquet"), emit: metadata_parquet, optional: true
    // GC factor tables (.parquet / .tsv.gz / .tsv, off- and on-target). The tool-level
    // subworkflow stages them beside the BED so each tool uses extract's factors, as run-all does.
    tuple val(meta), path("*.correction_factors*"), emit: gc_factors, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""
    // run-all's read filters, so the BED every tool reads is the one run-all extracts.
    def mapq_arg = params.mapq != 20 ? "--mapq ${params.mapq}" : ""
    def minlen_arg = params.minlen != 65 ? "--minlen ${params.minlen}" : ""
    def maxlen_arg = params.maxlen != 1000 ? "--maxlen ${params.maxlen}" : ""
    def skip_dup_arg = params.skip_duplicates == false ? "--no-skip-duplicates" : ""
    def proper_pair_arg = params.require_proper_pair == false ? "--no-require-proper-pair" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def output_format_arg = params.output_format && params.output_format != 'tsv' ? "--output-format ${params.output_format}" : ""
    def compress_arg = params.compress_tsv ? "--compress" : ""

    """
    krewlyzer extract \\
        -i $bam \\
        -r $fasta \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        $genome_arg \\
        $assay_arg \\
        $targets_arg \\
        $skip_targets_arg \\
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
    touch ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi
    printf 'sample_id\\ttotal_fragments\\n${prefix}\\t1000\\n' > ${prefix}.metadata.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.9.2
    END_VERSIONS
    """
}
