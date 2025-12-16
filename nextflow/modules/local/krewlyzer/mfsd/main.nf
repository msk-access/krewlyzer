process KREWLYZER_MFSD {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.2.3"

    input:
    tuple val(meta), path(bam), path(bai), path(variants)

    output:
    tuple val(meta), path("*.mFSD.tsv")     , emit: mfsd
    tuple val(meta), path("*.filtered.maf") , emit: filtered_maf, optional: true
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    krewlyzer mfsd \\
        $bam \\
        --input-file $variants \\
        --output ./ \\
        --sample-name $prefix \\
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
}
