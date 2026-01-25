/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_RUNALL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run all krewlyzer fragmentomics feature extraction tools on a single BAM.
    
    Outputs: FSC, FSR, FSD, WPS, OCF, Motif (EndMotif, BreakpointMotif, MDS)
    Optional: mFSD (with --variants), UXM (with bisulfite BAM)
    
    Panel mode: Use meta.assay (xs1/xs2) for gene-centric FSC and dual WPS output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_RUNALL {
    tag "$meta.id"
    label 'process_high'

    container "ghcr.io/msk-access/krewlyzer:0.3.2"

    input:
    tuple val(meta), path(bam), path(bai), path(variants), path(pon), path(targets)
    path fasta

    output:
    tuple val(meta), path("*.{txt,tsv,bed.gz,tsv.gz,parquet}"), emit: results
    tuple val(meta), path("*.FSC.gene.tsv"), emit: gene_fsc, optional: true
    tuple val(meta), path("*.metadata.json")    , emit: metadata, optional: true
    tuple val(meta), path("*.features.json")    , emit: unified_json, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def variant_arg = variants ? "--variants ${variants}" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def pon_arg = pon ? "--pon-model ${pon}" : ""
    def skip_pon_arg = params.skip_pon ? "--skip-pon" : ""
    def no_tfbs_arg = params.no_tfbs ? "--no-tfbs" : ""
    def no_atac_arg = params.no_atac ? "--no-atac" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def json_arg = params.generate_json ? "--generate-json" : ""
    def format_arg = params.output_format != 'auto' ? "--output-format ${params.output_format}" : ""
    def maxlen_arg = params.maxlen != 1000 ? "--maxlen ${params.maxlen}" : ""
    
    // Construct CLI command
    """
    krewlyzer run-all \\
        -i $bam \\
        -r $fasta \\
        --output ./ \\
        --threads $task.cpus \\
        --sample-name $prefix \\
        $variant_arg \\
        $targets_arg \\
        $assay_arg \\
        $genome_arg \\
        $gc_arg \\
        $pon_arg \\
        $skip_pon_arg \\
        $no_tfbs_arg \\
        $no_atac_arg \\
        $verbose_arg \\
        $json_arg \\
        $format_arg \\
        $maxlen_arg \\
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
    touch ${prefix}.FSC.tsv
    touch ${prefix}.FSR.tsv
    touch ${prefix}.FSD.tsv
    touch ${prefix}.WPS.parquet
    touch ${prefix}.WPS_background.parquet
    touch ${prefix}.OCF.tsv
    touch ${prefix}.EndMotif.tsv
    touch ${prefix}.BreakPointMotif.tsv
    touch ${prefix}.MDS.tsv
    touch ${prefix}.correction_factors.tsv
    echo '{"sample_id":"${prefix}","total_fragments":0}' > ${prefix}.metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.3.2
    END_VERSIONS
    """
}
