/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_MFSD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mutant Fragment Size Distribution - variant-level fragment profiles.
    Compares fragment sizes near somatic variants vs. background.
    
    Supports GC correction via --reference or --correction-factors.
    Supports duplex weighting for consensus BAMs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_MFSD {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.6.0"

    input:
    tuple val(meta), path(bam), path(bai), path(variants)
    path fasta                    // Optional: reference FASTA for GC correction
    path correction_factors       // Optional: pre-computed GC correction factors

    output:
    tuple val(meta), path("*.mFSD.tsv")              , emit: mfsd
    tuple val(meta), path("*.mFSD.distributions.tsv"), emit: distributions, optional: true
    tuple val(meta), path("*.filtered.maf")          , emit: filtered_maf, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // GC correction: prefer pre-computed factors, fallback to FASTA
    def gc_factors_arg = correction_factors ? "--correction-factors ${correction_factors}" : ""
    def ref_arg = (!correction_factors && fasta) ? "--reference ${fasta}" : ""
    
    // Filter and algorithm settings (use pipeline params with defaults)
    def mapq_arg = params.mapq != 20 ? "--mapq ${params.mapq}" : ""
    def minlen_arg = params.minlen != 65 ? "--minlen ${params.minlen}" : ""
    def maxlen_arg = params.maxlen != 1000 ? "--maxlen ${params.maxlen}" : ""
    def min_baseq_arg = params.min_baseq != 20 ? "--min-baseq ${params.min_baseq}" : ""
    def skip_dup_arg = params.skip_duplicates == false ? "--no-skip-duplicates" : ""
    def proper_pair_arg = params.require_proper_pair == false ? "--no-require-proper-pair" : ""
    def duplex_arg = params.duplex ? "--duplex" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""

    """
    krewlyzer mfsd \\
        -i $bam \\
        -V $variants \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        --output-distributions \\
        $gc_factors_arg \\
        $ref_arg \\
        $mapq_arg \\
        $minlen_arg \\
        $maxlen_arg \\
        $min_baseq_arg \\
        $skip_dup_arg \\
        $proper_pair_arg \\
        $duplex_arg \\
        $verbose_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "variant_id\\tmean_size\\tmedian_size" > ${prefix}.mFSD.tsv
    echo -e "variant_id\\tsize\\tcount" > ${prefix}.mFSD.distributions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.6.0
    END_VERSIONS
    """
}
