/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_REGION_MDS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-Region Motif Diversity Score (MDS) at exon/gene level.
    Calculates MDS for each gene's exons, with E1 (first exon) tracking.
    Based on Helzer et al. (2025) methodology.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_REGION_MDS {
    tag "$meta.id"
    label 'process_high'
    container "ghcr.io/msk-access/krewlyzer:0.5.1"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path pon_model   // Optional PON model for z-score normalization
    val assay        // Assay code (xs1, xs2, wgs) for bundled gene BED

    output:
    tuple val(meta), path("*.MDS.exon.tsv"), emit: mds_exon
    tuple val(meta), path("*.MDS.gene.tsv"), emit: mds_gene
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_arg = params.genome ? "--genome ${params.genome}" : "--genome hg19"
    def assay_arg = assay ? "--assay ${assay}" : ""
    def pon_arg = pon_model ? "--pon-model ${pon_model}" : ""
    def skip_pon_arg = params.skip_pon ? "--skip-pon" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def silent_arg = params.silent ? "--silent" : ""

    """
    krewlyzer region-mds \\
        $bam \\
        $fasta \\
        ./ \\
        $genome_arg \\
        $assay_arg \\
        $pon_arg \\
        $skip_pon_arg \\
        $verbose_arg \\
        $silent_arg \\
        $args

    # Rename outputs to include sample prefix if needed
    if [ -f "*.MDS.exon.tsv" ] && [ ! -f "${prefix}.MDS.exon.tsv" ]; then
        for f in *.MDS.exon.tsv; do
            mv "\$f" "${prefix}.MDS.exon.tsv"
        done
    fi
    if [ -f "*.MDS.gene.tsv" ] && [ ! -f "${prefix}.MDS.gene.tsv" ]; then
        for f in *.MDS.gene.tsv; do
            mv "\$f" "${prefix}.MDS.gene.tsv"
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "gene\\tname\\tchrom\\tstart\\tend\\tstrand\\tn_fragments\\tmds" > ${prefix}.MDS.exon.tsv
    echo -e "gene\\tn_exons\\tn_fragments\\tmds_mean\\tmds_e1\\tmds_std" > ${prefix}.MDS.gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.5.1
    END_VERSIONS
    """
}
