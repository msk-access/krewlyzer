process KREWLYZER_RUNALL {
    tag "$meta.id"
    label 'process_high'

    container "ghcr.io/msk-access/krewlyzer:0.3.2"

    input:
    tuple val(meta), path(bam), path(bai), path(variants), path(pon), path(targets)
    path fasta

    output:
    tuple val(meta), path("*.{txt,tsv,bed.gz,tsv.gz,parquet}"), emit: results
    tuple val(meta), path("*.metadata.json")    , emit: metadata, optional: true
    tuple val(meta), path("*.features.json")    , emit: unified_json, optional: true
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def variant_arg = variants ? "--variants ${variants}" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def gc_arg = params.gc_correct == false ? "--no-gc-correct" : ""
    def pon_arg = pon ? "--pon-model ${pon}" : ""
    def verbose_arg = params.verbose ? "--verbose" : ""
    def json_arg = params.generate_json ? "--generate-json" : ""
    def format_arg = params.output_format != 'auto' ? "--output-format ${params.output_format}" : ""
    
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
        $genome_arg \\
        $gc_arg \\
        $pon_arg \\
        $verbose_arg \\
        $json_arg \\
        $format_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
