/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_MFSD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mutant Fragment Size Distribution - variant-level fragment profiles.
    Compares fragment sizes near somatic variants vs. background.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_MFSD {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.3.2"

    input:
    tuple val(meta), path(bam), path(bai), path(variants)

    output:
    tuple val(meta), path("*.mFSD.tsv")              , emit: mfsd
    tuple val(meta), path("*.distributions.tsv")     , emit: distributions
    tuple val(meta), path("*.filtered.maf")          , emit: filtered_maf, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def verbose_arg = params.verbose ? "--verbose" : ""
    
    """
    krewlyzer mfsd \\
        -i $bam \\
        -V $variants \\
        --output ./ \\
        --sample-name $prefix \\
        --threads $task.cpus \\
        --output-distributions \\
        $verbose_arg \\
        $args

    # Copy filtered MAF to output if it exists (keeps the filtered subset with results)
    if [[ "$variants" == *.filtered.maf ]]; then
        cp $variants ${prefix}.filtered.maf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "variant_id\\tmean_size\\tmedian_size" > ${prefix}.mFSD.tsv
    echo -e "variant_id\\tsize\\tcount" > ${prefix}.distributions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.3.2
    END_VERSIONS
    """
}
