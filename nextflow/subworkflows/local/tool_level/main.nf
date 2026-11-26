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
    ch_extract    // [meta, bam, bai, pon, targets, wps_anchors]
    ch_beds       // [meta, bed, pon, targets, wps_anchors] - pre-extracted
    ch_methyl     // [meta, mbam, mbai]
    ch_mfsd       // [meta, bam, bai, maf] - already filtered; the mFSD BAM (meta.mfsd_duplex)
    fasta

    main:
    ch_versions = Channel.empty()

    // Every tool is given the sample's PON, targets and WPS anchors, as INPUT_CHECK resolved them
    // for run-all. They came from the pipeline-wide params for five tools, so a sample whose row
    // or assay named another PON or BED was computed against the wrong one.

    // =====================================================
    // 1. EXTRACTION (BAM → BED.gz)
    // =====================================================
    KREWLYZER_EXTRACT(
        ch_extract.map { meta, bam, bai, pon, targets, anchors -> [meta, bam, bai, targets] },
        fasta
    )
    ch_versions = ch_versions.mix(KREWLYZER_EXTRACT.out.versions)

    // Combine extracted BEDs with extract's GC factor tables and the sample's settings.
    // The tables are staged beside the BED in every tool's work directory, so each tool resolves
    // them as run-all does. Without them every tool runs without GC correction: only extraction
    // fits a table.
    ch_extracted = KREWLYZER_EXTRACT.out.bed
        .join(KREWLYZER_EXTRACT.out.gc_factors, remainder: true)
        .join(ch_extract.map { meta, bam, bai, pon, targets, anchors -> [meta, pon, targets, anchors] })
        .map { meta, bed, gc, pon, targets, anchors -> [meta, bed, gc ?: [], pon, targets, anchors] }

    // Merge with pre-extracted BEDs, which come without factor tables
    ch_all_beds = ch_extracted.mix(
        ch_beds.map { meta, bed, pon, targets, anchors -> [meta, bed, [], pon, targets, anchors] }
    )

    // =====================================================
    // 2. PARALLEL FEATURE EXTRACTION
    // =====================================================

    // [meta, bed, gc_factors, pon, targets]: what every BED-reading tool but WPS takes
    ch_bed_settings = ch_all_beds.map { meta, bed, gc, pon, targets, anchors -> [meta, bed, gc, pon, targets] }

    // FSC - Fragment Size Coverage
    KREWLYZER_FSC(ch_bed_settings)
    ch_versions = ch_versions.mix(KREWLYZER_FSC.out.versions)

    // FSR - Fragment Size Ratio
    KREWLYZER_FSR(ch_bed_settings)
    ch_versions = ch_versions.mix(KREWLYZER_FSR.out.versions)

    // FSD - Fragment Size Distribution
    KREWLYZER_FSD(ch_bed_settings)
    ch_versions = ch_versions.mix(KREWLYZER_FSD.out.versions)

    // WPS - Windowed Protection Score (also takes the sample's anchors)
    KREWLYZER_WPS(
        ch_all_beds,
        fasta,
        params.wps_background ? file(params.wps_background) : []
    )
    ch_versions = ch_versions.mix(KREWLYZER_WPS.out.versions)

    // OCF - Orientation-aware cfDNA Fragmentation
    KREWLYZER_OCF(ch_bed_settings, fasta)
    ch_versions = ch_versions.mix(KREWLYZER_OCF.out.versions)

    // REGION_ENTROPY - TFBS/ATAC size entropy
    KREWLYZER_REGION_ENTROPY(ch_bed_settings, fasta)
    ch_versions = ch_versions.mix(KREWLYZER_REGION_ENTROPY.out.versions)

    // MOTIF - End motif and MDS, read from the BAM as run-all reads them during extraction.
    // It was handed the BED, in a tuple one item short of its input, so it could not run.
    ch_bam_settings = ch_extract.map { meta, bam, bai, pon, targets, anchors -> [meta, bam, bai, pon, targets] }
    KREWLYZER_MOTIF(ch_bam_settings, fasta)
    ch_versions = ch_versions.mix(KREWLYZER_MOTIF.out.versions)

    // REGION_MDS - Per-gene/exon Motif Diversity Score (BAM-based; the assay comes from meta).
    // Only for a sample with an assay, as run-all: its gene BED is the assay's. Without one the
    // tool falls back to the WGS gene BED and wrote tables run-all never writes.
    KREWLYZER_REGION_MDS(
        ch_extract
            .filter { meta, bam, bai, pon, targets, anchors -> meta.assay && !(meta.assay in ['UNKNOWN', 'WGS']) }
            .map { meta, bam, bai, pon, targets, anchors -> [meta, bam, bai, pon] },
        fasta
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
    // mFSD: no factor tables in tool-level mode yet, where run-all passes extract's (on-target
    // first on a panel). How a BAM-only tool chooses its table is decided in the fragment-model
    // Step 4 (mutant and wild-type fragments GC-matched by locus); docs/nextflow/outputs.md
    // lists the difference.
    KREWLYZER_MFSD(ch_mfsd, fasta, [])
    ch_versions = ch_versions.mix(KREWLYZER_MFSD.out.versions)

    emit:
    fsc      = KREWLYZER_FSC.out.table
    fsr      = KREWLYZER_FSR.out.table
    fsd      = KREWLYZER_FSD.out.table
    wps      = KREWLYZER_WPS.out.parquet
    ocf      = KREWLYZER_OCF.out.table
    mds      = KREWLYZER_MOTIF.out.mds
    mfsd     = KREWLYZER_MFSD.out.mfsd
    uxm      = KREWLYZER_UXM.out.table
    versions = ch_versions
}
