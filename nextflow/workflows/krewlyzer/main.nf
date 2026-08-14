/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow for krewlyzer fragmentomics pipeline.
    
    Modes:
    - use_runall = true:  Unified run-all command per sample
    - use_runall = false: Individual tool-level modules (legacy)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../../subworkflows/local/input_check/main'
include { TOOL_LEVEL } from '../../subworkflows/local/tool_level/main'
include { KREWLYZER_RUNALL } from '../../modules/local/krewlyzer/runall/main'
include { FILTER_MAF } from '../../modules/local/krewlyzer/filter_maf/main'
include { SPLIT_MAF }  from '../../modules/local/krewlyzer/split_maf/main'
include { KREWLYZER_VALIDATE_PON    } from '../../modules/local/krewlyzer/validate_pon/main'
include { KREWLYZER_VALIDATE_COHORT } from '../../modules/local/krewlyzer/validate_cohort/main'

workflow KREWLYZER {
    take:
    samplesheet
    fasta

    main:
    ch_versions = Channel.empty()

    // Parquet is the default as of 0.9.0, so this only fires when someone has
    // deliberately excluded it. Both consequences are silent downstream --
    // every reader there swallows exceptions and yields an empty feature dict
    // -- so they are worth saying once, where they can still be changed.
    if (!params.output_format || params.output_format == 'tsv') {
        log.warn "output_format is '${params.output_format}': this cohort will " +
                 "contain no Parquet, which is the only format the downstream " +
                 "consumer reads, and per-sample output validation will be " +
                 "skipped. Use 'parquet' or 'both' for data destined for modelling."
    }

    // =====================================================
    // 1. PARSE AND VALIDATE SAMPLESHEET
    // =====================================================
    INPUT_CHECK(samplesheet)

    // =====================================================
    // 2. MODE SELECTION
    // =====================================================
    
    if (params.use_runall) {
        // =====================================================
        // RUN-ALL MODE: FILTER_MAF → map → RUNALL (join-free)
        // File paths ride through FILTER_MAF as meta fields,
        // then .map() reconstructs the RUNALL tuple.
        // Guaranteed streaming — no join operator blocking.
        // =====================================================

        // --- MAF preparation -------------------------------------------
        //
        // Samples that share a MAF are split in ONE pass instead of one job
        // each. FILTER_MAF re-scanned the whole cohort MAF per sample, so a
        // 16,552-sample run did 16,552 full scans of one file, each a separate
        // SLURM job, while holding half the queue upstream of RUNALL.
        //
        // Grouped on (maf, single_sample) so only samples that would have done
        // identical work are batched. Per-sample MAFs yield one task each and
        // behave exactly as before. Samples with no MAF keep FILTER_MAF, whose
        // Mode 3 already writes the header-only placeholder.
        ch_maf_split = INPUT_CHECK.out.maf_for_filter
            .branch { meta, maf ->
                has_maf: maf
                no_maf : true
            }

        // groupTuple is a barrier, and deliberately so: this channel comes
        // straight from the samplesheet, so it completes in seconds, and
        // nothing downstream could start before the MAFs existed anyway.
        // The sample ids go to SPLIT_MAF as a FILE, never interpolated into its
        // script. Injecting them cost a 16,552-sample run: Nextflow strips a
        // script block's common indentation *after* interpolation, so 16,552
        // zero-indent id lines left no common prefix, the Python kept its block
        // indent, and the task died on line 2 with IndentationError.
        //
        // The filename is derived from the group key, so it is stable across
        // runs and `-resume` still hits. Written here rather than in a process
        // because it is one small text file per distinct MAF, not work.
        ch_grouped = ch_maf_split.has_maf
            .map { meta, maf -> [ [maf, meta.single_sample ?: false], meta ] }
            .groupTuple()
            .map { key, metas ->
                def dir = file("${workflow.workDir}/split-maf-ids")
                dir.mkdirs()
                def ids = file("${dir}/${Math.abs(key.toString().hashCode())}.ids.txt")
                ids.text = metas.collect { it.id }.join('\n') + '\n'
                [ metas, key[1], key[0], ids ]
            }

        SPLIT_MAF(ch_grouped)
        ch_versions = ch_versions.mix(SPLIT_MAF.out.versions.first())

        // The metas ride through SPLIT_MAF, so re-pairing is a map and not a
        // join -- this workflow is deliberately join-free so RUNALL streams.
        // SPLIT_MAF guarantees one file per requested sample, so a null here is
        // a bug in that contract, not an expected absence: fail loudly rather
        // than drop the sample and shorten the cohort in silence.
        ch_split_paired = SPLIT_MAF.out.maf
            .flatMap { metas, files ->
                def as_list = files instanceof List ? files : [files]
                def by_id = as_list.collectEntries { f ->
                    [ (f.name - '.filtered.maf'): f ]
                }
                metas.collect { meta ->
                    def f = by_id[meta.id]
                    if (f == null) {
                        error "SPLIT_MAF emitted no file for '${meta.id}'. One " +
                              "file per requested sample is its contract; a " +
                              "missing file would silently drop the sample."
                    }
                    [meta, f]
                }
            }

        // --- FILTER_MAF: only samples with no MAF (header-only placeholder) ---
        FILTER_MAF(ch_maf_split.no_maf)
        ch_versions = ch_versions.mix(FILTER_MAF.out.versions.first())

        // --- Reconstruct RUNALL tuple from meta fields + filtered MAF ---
        ch_runall = ch_split_paired.mix(FILTER_MAF.out.maf)
            .map { meta, filtered_maf ->
                [meta, meta._bam, meta._bai,
                 meta._mfsd_bam ?: [], meta._mfsd_bai ?: [],
                 meta._bisulfite ?: [],
                 filtered_maf,
                 meta._pon ?: [], meta._targets ?: [],
                 meta._wps_anchors ?: [], meta._wps_bg ?: []]
            }

        KREWLYZER_RUNALL(ch_runall, fasta)
        ch_versions = ch_versions.mix(KREWLYZER_RUNALL.out.versions.first())
        
        // Map outputs
        ch_fsc = KREWLYZER_RUNALL.out.fsc
        ch_fsr = KREWLYZER_RUNALL.out.fsr
        ch_fsd = KREWLYZER_RUNALL.out.fsd
        ch_wps = KREWLYZER_RUNALL.out.wps
        ch_ocf = KREWLYZER_RUNALL.out.ocf
        ch_mds = KREWLYZER_RUNALL.out.mds
        ch_features_json = KREWLYZER_RUNALL.out.features_json
        ch_fingerprints  = KREWLYZER_RUNALL.out.fingerprint.map { _meta, fp -> fp }
        
    } else {
        // =====================================================
        // TOOL-LEVEL MODE: Individual modules
        // =====================================================
        
        // Filter multi-sample MAFs (only needed for tool-level mode)
        FILTER_MAF(INPUT_CHECK.out.maf_multi.map { meta, bam, bai, maf -> [meta, maf] })
        ch_versions = ch_versions.mix(FILTER_MAF.out.versions.first())
        
        // Filter out empty MAFs (header-only files < 100 bytes)
        ch_filtered_valid = FILTER_MAF.out.maf
            .filter { meta, filtered_maf -> filtered_maf.size() > 100 }
        
        ch_mfsd_filtered = INPUT_CHECK.out.maf_multi
            .map { meta, bam, bai, maf -> [meta.id, meta, bam, bai] }
            .join(ch_filtered_valid.map { meta, filtered_maf -> [meta.id, filtered_maf] })
            .map { id, meta, bam, bai, filtered_maf -> [meta, bam, bai, filtered_maf] }
        
        ch_mfsd_all = ch_mfsd_filtered.mix(INPUT_CHECK.out.maf_single)
        
        TOOL_LEVEL(
            INPUT_CHECK.out.extract,
            INPUT_CHECK.out.beds,
            INPUT_CHECK.out.methyl,
            ch_mfsd_all,
            fasta
        )
        ch_versions = ch_versions.mix(TOOL_LEVEL.out.versions)
        
        // Map outputs
        ch_fsc = TOOL_LEVEL.out.fsc
        ch_fsr = TOOL_LEVEL.out.fsr
        ch_fsd = TOOL_LEVEL.out.fsd
        ch_wps = TOOL_LEVEL.out.wps
        ch_ocf = TOOL_LEVEL.out.ocf
        ch_mds = TOOL_LEVEL.out.mds
        ch_features_json = Channel.empty()
        // Tool-level mode has no run-all step to emit fingerprints. Cohort
        // validation is skipped rather than run on a partial set, which would
        // report degeneracy findings that are artefacts of the missing samples.
        ch_fingerprints = Channel.empty()
    }

    // =====================================================
    // PON VALIDATION
    // =====================================================
    // Check the reference before anything is measured against it. Runs once
    // per pipeline rather than once per sample -- the PON is the same file for
    // all of them.
    //
    // Wired here rather than added as a module and left to be called later:
    // an unwired module is exactly the defect #34 fixed, and it looks like
    // coverage while providing none.
    ch_pon_report = Channel.empty()
    if (params.pon_model) {
        KREWLYZER_VALIDATE_PON(Channel.fromPath(params.pon_model, checkIfExists: true))
        ch_versions   = ch_versions.mix(KREWLYZER_VALIDATE_PON.out.versions)
        ch_pon_report = KREWLYZER_VALIDATE_PON.out.report
    }

    // =====================================================
    // COHORT VALIDATION (gather)
    // =====================================================
    // Degeneracy is inherently cross-sample: every sample can pass on its own
    // while a metric is constant across all of them. This is the gather half
    // of that check, and it was never wired up -- the module existed and
    // nothing called it.
    KREWLYZER_VALIDATE_COHORT(ch_fingerprints.collect(sort: true))
    ch_versions = ch_versions.mix(KREWLYZER_VALIDATE_COHORT.out.versions)

    emit:
    fsc           = ch_fsc
    fsr           = ch_fsr
    fsd           = ch_fsd
    wps           = ch_wps
    ocf           = ch_ocf
    mds           = ch_mds
    features_json = ch_features_json
    fingerprints  = ch_fingerprints
    cohort_report = KREWLYZER_VALIDATE_COHORT.out.report
    pon_report    = ch_pon_report
    versions      = ch_versions
}
