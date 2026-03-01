/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TOOL_LEVEL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Individual tool-level subworkflow for krewlyzer pipeline.
    
    Runs each fragmentomics tool as a separate process.
    Supports pre-extracted BED files (skips extraction).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { KREWLYZER_EXTRACT } from '../../../modules/local/krewlyzer/extract/main'
include { KREWLYZER_FSC } from '../../../modules/local/krewlyzer/fsc/main'
include { KREWLYZER_FSR } from '../../../modules/local/krewlyzer/fsr/main'
include { KREWLYZER_FSD } from '../../../modules/local/krewlyzer/fsd/main'
include { KREWLYZER_WPS } from '../../../modules/local/krewlyzer/wps/main'
include { KREWLYZER_OCF } from '../../../modules/local/krewlyzer/ocf/main'
include { KREWLYZER_MOTIF } from '../../../modules/local/krewlyzer/motif/main'
include { KREWLYZER_MFSD } from '../../../modules/local/krewlyzer/mfsd/main'
include { KREWLYZER_REGION_ENTROPY } from '../../../modules/local/krewlyzer/region_entropy/main'
include { KREWLYZER_REGION_MDS } from '../../../modules/local/krewlyzer/region_mds/main'
include { KREWLYZER_UXM } from '../../../modules/local/krewlyzer/uxm/main'

workflow TOOL_LEVEL {
    take:
    ch_extract    // [meta, bam, bai, pon, targets]
    ch_beds       // [meta, bed, pon, targets] - pre-extracted
    ch_methyl     // [meta, mbam, mbai]
    ch_mfsd       // [meta, bam, bai, maf] - already filtered
    fasta

    main:
    ch_versions = Channel.empty()

    // =====================================================
    // 1. EXTRACTION (BAM → BED.gz)
    // =====================================================
    ch_extract_input = ch_extract.map { meta, bam, bai, pon, targets -> 
        [meta, bam, bai, targets ?: []] 
    }
    KREWLYZER_EXTRACT(
        ch_extract_input.map { meta, bam, bai, targets -> [meta, bam, bai] },
        fasta,
        ch_extract_input.map { it[3] }  // targets
    )
    ch_versions = ch_versions.mix(KREWLYZER_EXTRACT.out.versions)
    
    // Combine extracted BEDs with PON/targets for downstream
    ch_extracted_with_meta = KREWLYZER_EXTRACT.out.bed
        .join(ch_extract.map { meta, bam, bai, pon, targets -> [meta, pon, targets] })
        .map { meta, bed, pon, targets -> [meta, bed, pon ?: [], targets ?: []] }
    
    // Merge extracted BEDs with pre-extracted BEDs
    ch_all_beds = ch_extracted_with_meta.mix(ch_beds)

    // =====================================================
    // 2. PARALLEL FEATURE EXTRACTION
    // =====================================================
    
    // Convenience channels
    ch_bed_only = ch_all_beds.map { meta, bed, pon, targets -> [meta, bed] }
    ch_bed_targets = ch_all_beds.map { meta, bed, pon, targets -> [meta, bed, targets ?: []] }
    
    // FSC - Fragment Size Coverage
    KREWLYZER_FSC(
        ch_bed_targets.map { meta, bed, targets -> [meta, bed] }, 
        ch_bed_targets.map { it[2] }
    )
    ch_versions = ch_versions.mix(KREWLYZER_FSC.out.versions)
    
    // FSR - Fragment Size Ratio
    KREWLYZER_FSR(
        ch_bed_targets.map { meta, bed, targets -> [meta, bed] }, 
        ch_bed_targets.map { it[2] }
    )
    ch_versions = ch_versions.mix(KREWLYZER_FSR.out.versions)
    
    // FSD - Fragment Size Distribution
    KREWLYZER_FSD(ch_bed_only)
    ch_versions = ch_versions.mix(KREWLYZER_FSD.out.versions)
    
    // WPS - Windowed Protection Score
    KREWLYZER_WPS(
        ch_bed_only,
        fasta,
        params.wps_anchors ? file(params.wps_anchors) : [],
        params.wps_background ? file(params.wps_background) : []
    )
    ch_versions = ch_versions.mix(KREWLYZER_WPS.out.versions)
    
    // OCF - Orientation-aware cfDNA Fragmentation
    KREWLYZER_OCF(ch_bed_only, fasta)
    ch_versions = ch_versions.mix(KREWLYZER_OCF.out.versions)
    
    // MOTIF - End motif and MDS
    KREWLYZER_MOTIF(ch_bed_only, fasta)
    ch_versions = ch_versions.mix(KREWLYZER_MOTIF.out.versions)
    
    // REGION_ENTROPY - TFBS/ATAC size entropy
    ch_bed_pon_targets = ch_all_beds.map { meta, bed, pon, targets -> 
        [meta, bed, pon ?: [], targets ?: []] 
    }
    KREWLYZER_REGION_ENTROPY(
        ch_bed_pon_targets.map { meta, bed, pon, targets -> [meta, bed] },
        fasta,
        ch_bed_pon_targets.map { it[2] },  // pon
        ch_bed_pon_targets.map { it[3] }   // targets
    )
    ch_versions = ch_versions.mix(KREWLYZER_REGION_ENTROPY.out.versions)

    // REGION_MDS - Per-gene/exon Motif Diversity Score (BAM-based, needs assay)
    ch_bam_for_mds = ch_extract.map { meta, bam, bai, pon, targets -> 
        [meta, bam, bai, pon ?: [], meta.assay ?: ''] 
    }
    KREWLYZER_REGION_MDS(
        ch_bam_for_mds.map { meta, bam, bai, pon, assay -> [meta, bam, bai] },
        fasta,
        ch_bam_for_mds.map { it[3] },  // pon
        ch_bam_for_mds.map { it[4] }   // assay
    )
    ch_versions = ch_versions.mix(KREWLYZER_REGION_MDS.out.versions)

    // =====================================================
    // 3. METHYLATION (UXM)
    // =====================================================
    KREWLYZER_UXM(ch_methyl, fasta)
    ch_versions = ch_versions.mix(KREWLYZER_UXM.out.versions)

    // =====================================================
    // 4. mFSD (Mutant Fragment Size Distribution)
    // =====================================================
    // mFSD: pass FASTA for runtime GC correction; pre-computed factors are not
    // available in tool-level mode (run-all generates them during extraction)
    KREWLYZER_MFSD(ch_mfsd, fasta, [])
    ch_versions = ch_versions.mix(KREWLYZER_MFSD.out.versions)

    emit:
    fsc      = KREWLYZER_FSC.out.tsv
    fsr      = KREWLYZER_FSR.out.tsv
    fsd      = KREWLYZER_FSD.out.tsv
    wps      = KREWLYZER_WPS.out.parquet
    ocf      = KREWLYZER_OCF.out.tsv
    mds      = KREWLYZER_MOTIF.out.mds
    mfsd     = KREWLYZER_MFSD.out.mfsd
    uxm      = KREWLYZER_UXM.out.tsv
    versions = ch_versions
}
