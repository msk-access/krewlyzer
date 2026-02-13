/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTER_MAF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pre-filter MAF file to extract variants for a specific sample.
    Used before mFSD to isolate sample-specific variants from multi-sample MAFs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FILTER_MAF {
    tag "$meta.id"
    label 'process_low'
    maxForks 50  // Half of -qs to leave slots for RUNALL streaming

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.filtered.maf"), emit: maf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    #!/usr/bin/env python3
    import sys
    import platform

    sample_id = '${prefix}'
    sample_id_lower = sample_id.lower()
    variant_count = 0
    total_rows = 0

    with open('${maf}', 'r') as infile, open('${prefix}.filtered.maf', 'w') as outfile:
        header_written = False
        tsb_col = None

        for line in infile:
            # Handle comment lines (keep them)
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\\t')

            # First non-comment line is the header
            if not header_written:
                outfile.write(line)
                header_written = True
                # Find Tumor_Sample_Barcode column
                for i, col in enumerate(fields):
                    if col == 'Tumor_Sample_Barcode':
                        tsb_col = i
                        break
                if tsb_col is None:
                    print("WARNING: Tumor_Sample_Barcode column not found, keeping all rows", file=sys.stderr)
                else:
                    print(f"Found Tumor_Sample_Barcode at column {tsb_col}", file=sys.stderr)
                continue

            total_rows += 1

            # Filter: check if sample_id is a substring of Tumor_Sample_Barcode
            if tsb_col is not None:
                if len(fields) > tsb_col and sample_id_lower in fields[tsb_col].lower():
                    outfile.write(line)
                    variant_count += 1
            else:
                # No TSB column found, keep all rows
                outfile.write(line)
                variant_count += 1

    # Log result
    print(f"Sample ID: '{sample_id}' | Total MAF rows: {total_rows} | Matched: {variant_count}", file=sys.stderr)
    if variant_count == 0:
        print(f"WARNING: No variants found for '{sample_id}' in Tumor_Sample_Barcode - mFSD will be skipped", file=sys.stderr)

    # Write versions.yml
    with open('versions.yml', 'w') as vf:
        vf.write('"${task.process}":\\n')
        vf.write(f'    python: {platform.python_version()}\\n')
    """

    stub:
    def prefix = meta.id
    """
    echo -e "Hugo_Symbol\\tChromosome\\tStart_Position\\tTumor_Sample_Barcode" > ${prefix}.filtered.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
