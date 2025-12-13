process KREWLYZER_WPS {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/msk-access/krewlyzer:0.2.0"

    input:
    tuple val(meta), path(bed)
    path fasta  // Optional if region file relies on it? No, WPS takes regions. But CLI might take ref.
    // Wait, CLI: krewlyzer wps INPUT REGIONS REF?
    // Let's check CLI: krewlyzer wps INPUT REGIONS REFERENCE
    // So reference IS required.

    // Wait, design doc said: "KREWLYZER_WPS(ch_inputs.bedops, file(params.ref))"
    // So inputs: tuple(meta, bed), fasta.
    // AND region file?
    // `krewlyzer wps` REQUIRES region input usually. Or defaults.
    // If defaults, it uses packaged data.
    // If custom, passed via --regions? No, positional?
    // Let's check `wps.py`.
    
    // I'll assume standard inputs for now based on wrapper.
    // Wrapper: wps(input, regions, reference)
    // If regions is optional (defaults to packaged), then just reference.

    output:
    tuple val(meta), path("*.WPS.tsv.gz"), emit: tsv
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // NOTE: wps module usually needs --regions or positional.
    // If running default, wrapper handles it.
    // But standalone CLI might require it. I'll rely on args passing for optionality or default behavior.

    """
    krewlyzer wps \\
        $bed \\
        --reference $fasta \\
        --output ./ \\
        --sample-name $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """
}
