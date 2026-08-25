/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_VALIDATE_PON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check the PON before a single sample is scored against it.

    KREWLYZER_VALIDATE_SAMPLE and _COHORT gate the results of a run. Nothing
    gated the reference those results are measured against -- and every PON
    defect found in 0.9.0 sat in four shipped files, visible in the file and
    invisible to every check: a wps_background block hardcoded to
    167.0/5.0/0.0/1.0 identical across all 28 groups and all four models,
    sigma floors that turned "no spread measured" into a divisor, four blocks
    built and read by nothing, and no record of what any of them was built
    from.

    Runs once per pipeline, not once per sample: the PON is the same file for
    every sample, so checking it 14,000 times would say the same thing 14,000
    times.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_VALIDATE_PON {
    tag "${pon_model.name}"
    label 'process_single'
    container "ghcr.io/msk-access/krewlyzer:0.9.2"

    input:
    path pon_model

    output:
    path "pon.validation.json", emit: report
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    krewlyzer validate-pon \\
        ${pon_model} \\
        --json-report pon.validation.json \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer version: //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"schema_version": "1", "findings": [], "exit_code": 0}' \\
        > pon.validation.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.9.2
    END_VERSIONS
    """
}
