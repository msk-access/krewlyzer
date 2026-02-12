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

// Helper: Check if file has data (not just header)
def hasData = { file ->
    def dataLineCount = 0
    file.withReader { reader ->
        String line
        while ((line = reader.readLine()) != null && dataLineCount < 2) {
            if (!line.startsWith('#') && line.trim()) {
                dataLineCount++
            }
        }
    }
    return dataLineCount > 1
}

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
    // 2. MAF FILTERING (Shared for both modes)
    // =====================================================
    
    // Filter multi-sample MAFs
    FILTER_MAF(INPUT_CHECK.out.maf_multi.map { meta, bam, bai, maf -> [meta, maf] })
    
    // Validate filtered MAF has data
    ch_filtered_valid = FILTER_MAF.out.maf
        .filter { meta, filtered_maf ->
            if (!hasData(filtered_maf)) {
                log.warn "Sample ${meta.id}: No variants after MAF filtering - skipping mFSD"
                return false
            }
            return true
        }
    
    // Join filtered MAF back with BAM info
    ch_mfsd_filtered = INPUT_CHECK.out.maf_multi
        .map { meta, bam, bai, maf -> [meta.id, meta, bam, bai] }
        .join(ch_filtered_valid.map { meta, filtered_maf -> [meta.id, filtered_maf] })
        .map { id, meta, bam, bai, filtered_maf -> [meta, bam, bai, filtered_maf] }
    
    // Combine with single-sample MAFs (bypass filtering)
    ch_mfsd_all = ch_mfsd_filtered.mix(INPUT_CHECK.out.maf_single)

    // =====================================================
    // 3. MODE SELECTION
    // =====================================================
    
    if (params.use_runall) {
        // =====================================================
        // RUN-ALL MODE: Unified command per sample
        // =====================================================
        
        // INPUT_CHECK.out.runall already has: [meta, bam, bai, mfsd_bam, mfsd_bai, bisulfite_bam, variants, pon, targets, wps_anchors, wps_background]
        // For samples with multi-sample MAF, we need to replace the variants with the filtered MAF
        ch_runall_base = INPUT_CHECK.out.runall
            .map { meta, bam, bai, mfsd_bam, mfsd_bai, bisulfite_bam, variants, pon, targets, wps_anchors, wps_bg ->
                [meta.id, meta, bam, bai, mfsd_bam, mfsd_bai, bisulfite_bam, pon, targets, wps_anchors, wps_bg]
            }
        
        // Get filtered MAFs by sample ID
        ch_filtered_mafs = ch_mfsd_all.map { meta, bam, bai, maf -> [meta.id, maf] }
        
        // Join and reconstruct the full tuple for RUNALL
        ch_runall = ch_runall_base
            .join(ch_filtered_mafs)
            .map { id, meta, bam, bai, mfsd_bam, mfsd_bai, bisulfite_bam, pon, targets, wps_anchors, wps_bg, maf ->
                [meta, bam, bai, mfsd_bam ?: [], mfsd_bai ?: [], bisulfite_bam ?: [], maf ?: [], pon ?: [], targets ?: [], wps_anchors ?: [], wps_bg ?: []]
            }
        
        KREWLYZER_RUNALL(ch_runall, fasta)
        ch_versions = ch_versions.mix(KREWLYZER_RUNALL.out.versions)
        
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
        ch_features_json = Channel.empty()  // Not available in tool-level mode
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
