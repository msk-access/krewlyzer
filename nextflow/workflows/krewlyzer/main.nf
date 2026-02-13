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

workflow KREWLYZER {
    take:
    samplesheet
    fasta

    main:
    ch_versions = Channel.empty()

    // =====================================================
    // 1. PARSE AND VALIDATE SAMPLESHEET
    // =====================================================
    INPUT_CHECK(samplesheet)

    // =====================================================
    // 2. MODE SELECTION
    // =====================================================
    
    if (params.use_runall) {
        // =====================================================
        // RUN-ALL MODE: Direct pipe — each sample is independent
        // MAF filtering is inlined in the RUNALL process script
        // =====================================================
        
        ch_runall = INPUT_CHECK.out.runall
            .map { meta, bam, bai, mfsd_bam, mfsd_bai, bisulfite_bam, variants, pon, targets, wps_anchors, wps_bg ->
                [meta, bam, bai, mfsd_bam ?: [], mfsd_bai ?: [], bisulfite_bam ?: [], variants ?: [], pon ?: [], targets ?: [], wps_anchors ?: [], wps_bg ?: []]
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
    }

    emit:
    fsc           = ch_fsc
    fsr           = ch_fsr
    fsd           = ch_fsd
    wps           = ch_wps
    ocf           = ch_ocf
    mds           = ch_mds
    features_json = ch_features_json
    versions      = ch_versions
}
