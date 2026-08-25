/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_VALIDATE_SAMPLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Scatter half of the output-contract gate. Everything decidable from a
    single sample: that every table a consumer reads is present and correctly
    shaped, and that the domain invariants hold.

    Also emits a small fingerprint (a hash and a couple of counts per column,
    ~20 KB) which KREWLYZER_VALIDATE_COHORT reduces. Degeneracy -- a metric
    that never varies -- cannot be judged one sample at a time, so it is
    deliberately not attempted here.

    Exit codes: 0 satisfied, 1 contract violation, 2 structural (unreadable or
    incomplete input). 2 is retryable; 1 is not.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_VALIDATE_SAMPLE {
    tag "$meta.id"
    label 'process_single'
    container "ghcr.io/msk-access/krewlyzer:0.9.2"

    input:
    tuple val(meta), path(sample_dir)

    output:
    tuple val(meta), path("*.fingerprint.json"), emit: fingerprint
    tuple val(meta), path("*.validation.json") , emit: report
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # A contract violation must not abort the scatter: the cohort step needs
    # every fingerprint to judge degeneracy, and one bad sample is a finding,
    # not a reason to lose the other 13,999. Exit 2 (structural) does fail,
    # because an unreadable sample cannot produce a usable fingerprint.
    set +e
    krewlyzer validate-output \\
        $sample_dir \\
        --fingerprint-out ${prefix}.fingerprint.json \\
        --json-report ${prefix}.validation.json \\
        $args
    status=\$?
    set -e

    if [ "\$status" -eq 2 ]; then
        echo "krewlyzer validate-output: structural failure on ${prefix}" >&2
        exit 2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"fingerprint_version": "1", "sample": "${prefix}", "observations": {}}' \\
        > ${prefix}.fingerprint.json
    echo '{"schema_version": "1", "findings": [], "exit_code": 0}' \\
        > ${prefix}.validation.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.9.2
    END_VERSIONS
    """
}
