process KREWLYZER_RUNALL {
    tag "$meta.id"
    label 'process_high'

    container "ghcr.io/msk-access/krewlyzer:0.2.0"

    input:
    tuple val(meta), path(bam), path(bai), path(variants)
    path fasta
    path targets   // Optional: bin_input/arms_file

    output:
    tuple val(meta), path("*.{txt,tsv,csv,bed.gz,tsv.gz}"), emit: results
    tuple val(meta), path("*.json")             , emit: metadata, optional: true
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def variant_arg = variants ? "--variants ${variants}" : ""
    def targets_arg = targets ? "--bin-input ${targets}" : ""
    
    // Construct CLI command
    """
    krewlyzer run-all \\
        $bam \\
        --reference $fasta \\
        --output ./ \\
        --threads $task.cpus \\
        --sample-name $prefix \\
        $variant_arg \\
        $targets_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
