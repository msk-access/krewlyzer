/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_VALIDATE_COHORT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Gather half of the output-contract gate: the checks that need more than
    one sample.

    A single sample cannot distinguish "this metric is a constant" from "this
    is its value here", which is exactly how four of the five WPS_background
    columns reached production carrying no information at all. This step
    compares fingerprints across the cohort and fails when a metric never
    varies.

    Reads the ~20 KB fingerprints, never the outputs, so a 14k-sample cohort
    reduces from megabytes rather than terabytes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_VALIDATE_COHORT {
    tag "cohort"
    label 'process_single'
    container "ghcr.io/msk-access/krewlyzer:0.9.1"

    input:
    path fingerprints, stageAs: "fingerprints/*"

    output:
    path "cohort.validation.json", emit: report
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def min_samples = params.validate_min_samples ?: 3
    """
    krewlyzer validate-cohort \\
        fingerprints/ \\
        --min-samples ${min_samples} \\
        --json-report cohort.validation.json \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer version: //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"schema_version": "1", "findings": [], "exit_code": 0}' \\
        > cohort.validation.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.8.3
    END_VERSIONS
    """
}
