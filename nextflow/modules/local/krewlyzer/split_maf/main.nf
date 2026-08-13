/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_MAF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Split ONE cohort MAF into per-sample MAFs in a single pass.

    FILTER_MAF does the same work per sample: a full linear scan of the shared
    MAF to keep the rows matching one Tumor_Sample_Barcode. On a cohort where
    every row of the samplesheet points at the same MAF -- which is the normal
    case -- that is N full scans of one file, each as its own SLURM job whose
    scheduling overhead dwarfs the grep it performs. At 16,552 samples it also
    held `queue_size.intdiv(2)` = 250 of 500 slots, upstream of RUNALL, so the
    opening wave of the run was nothing but MAF filtering and no real work
    could start.

    One scan, one job, N outputs.

    Batching is keyed on (maf, single_sample) by the caller, so only samples
    that would have produced identical work are grouped. A sheet with per-sample
    MAFs yields one task per distinct MAF and degrades to the old behaviour
    rather than doing anything surprising.

    EVERY requested sample gets a file, including those with no matching rows.
    A sample with zero variants that received no file would be dropped by the
    downstream re-pairing, and RUNALL would never run for it -- a silent loss of
    exactly the kind `validate-output --expect` exists to catch. The header-only
    output preserves FILTER_MAF's behaviour for that case.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SPLIT_MAF {
    tag "${metas.size()} sample(s)"
    label 'process_low'

    input:
    tuple val(metas), val(single_sample), path(maf)

    output:
    tuple val(metas), path("split/*.filtered.maf"), emit: maf
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The metas ride through so the caller can re-pair by id in a map rather
    // than a join -- this workflow is deliberately join-free so RUNALL streams.
    def ids_arg = metas.collect { it.id }.join('\n')
    def single_flag = single_sample ? 'true' : 'false'
    """
    #!/usr/bin/env python3
    import os
    import platform
    import sys

    os.makedirs('split', exist_ok=True)

    sample_ids = [s for s in '''${ids_arg}'''.split('\\n') if s]
    single = '${single_flag}' == 'true'

    # Case-insensitive, matching FILTER_MAF. Keep the first spelling seen for
    # each id so output filenames match the samplesheet exactly.
    by_lower = {}
    for sid in sample_ids:
        by_lower.setdefault(sid.lower(), sid)

    handles = {}
    header = None
    tsb_col = None
    total_rows = 0
    matched = {sid: 0 for sid in sample_ids}

    def handle_for(sid):
        h = handles.get(sid)
        if h is None:
            h = open(os.path.join('split', sid + '.filtered.maf'), 'w')
            if header is not None:
                h.write(header)
            handles[sid] = h
        return h

    with open('${maf}', 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue

            if header is None:
                header = line
                fields = line.rstrip('\\n').split('\\t')
                for i, col in enumerate(fields):
                    if col == 'Tumor_Sample_Barcode':
                        tsb_col = i
                        break
                if tsb_col is None:
                    print('WARNING: Tumor_Sample_Barcode column not found, '
                          'keeping all rows for every sample', file=sys.stderr)
                else:
                    print('Found Tumor_Sample_Barcode at column %d' % tsb_col,
                          file=sys.stderr)
                continue

            total_rows += 1

            # single_sample_maf: every row belongs to every sample in the group.
            # Grouped on this flag by the caller, so the whole group agrees.
            if single or tsb_col is None:
                for sid in sample_ids:
                    handle_for(sid).write(line)
                    matched[sid] += 1
                continue

            fields = line.rstrip('\\n').split('\\t')
            if len(fields) <= tsb_col:
                continue
            sid = by_lower.get(fields[tsb_col].lower())
            if sid is not None:
                handle_for(sid).write(line)
                matched[sid] += 1

    for h in handles.values():
        h.close()

    # A sample with no matching rows still gets a header-only file. Without it
    # the re-pairing downstream drops the sample and RUNALL never runs for it.
    if header is None:
        header = 'Hugo_Symbol\\tChromosome\\tStart_Position\\tTumor_Sample_Barcode\\n'
    empty = 0
    for sid in sample_ids:
        path = os.path.join('split', sid + '.filtered.maf')
        if not os.path.exists(path):
            with open(path, 'w') as fh:
                fh.write(header)
            empty += 1

    print('Split %d MAF rows across %d sample(s); %d had no variants'
          % (total_rows, len(sample_ids), empty), file=sys.stderr)

    # Refuse to emit a short set. The downstream map pairs by sample id, and a
    # missing file there is indistinguishable from a sample nobody asked for.
    produced = len(os.listdir('split'))
    if produced != len(sample_ids):
        sys.exit('SPLIT_MAF produced %d file(s) for %d sample(s)'
                 % (produced, len(sample_ids)))

    with open('versions.yml', 'w') as vf:
        vf.write('"${task.process}":\\n')
        vf.write('    python: %s\\n' % platform.python_version())
    """

    stub:
    // `metas`, not `sample_ids` -- the input was renamed so the metas could
    // ride through and the caller could re-pair with a map instead of a join.
    // The script block was updated and this was not, and no unit test reaches
    // a stub block, so it took a `-stub-run` to surface it.
    """
    mkdir -p split
    for s in ${metas.collect { it.id }.join(' ')}; do
        echo -e "Hugo_Symbol\\tChromosome\\tStart_Position\\tTumor_Sample_Barcode" > split/\$s.filtered.maf
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
