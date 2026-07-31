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

        // --- FILTER_MAF: all samples (multi-MAF filters, single-MAF passes through, no-MAF creates placeholder) ---
        FILTER_MAF(INPUT_CHECK.out.maf_for_filter)
        ch_versions = ch_versions.mix(FILTER_MAF.out.versions.first())

        // --- Reconstruct RUNALL tuple from meta fields + filtered MAF ---
        ch_runall = FILTER_MAF.out.maf
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
    versions      = ch_versions
}
