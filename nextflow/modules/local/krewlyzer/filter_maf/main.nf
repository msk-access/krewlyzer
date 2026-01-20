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
    import re
    
    sample_pattern = re.compile(r'.*${prefix}.*', re.IGNORECASE)
    variant_count = 0
    
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
                continue
            
            # Filter data rows by sample pattern
            if tsb_col is not None:
                if len(fields) > tsb_col and sample_pattern.match(fields[tsb_col]):
                    outfile.write(line)
                    variant_count += 1
            else:
                # No TSB column found, keep all rows
                outfile.write(line)
                variant_count += 1
    
    # Log result
    if variant_count == 0:
        print(f"WARNING: No variants found for sample pattern '.*${prefix}.*' - MFSD will be skipped", file=sys.stderr)
    else:
        print(f"Filtered {variant_count} variants for sample pattern: .*${prefix}.*", file=sys.stderr)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
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
