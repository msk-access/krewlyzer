#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import Modules
include { KREWLYZER_RUNALL } from './modules/local/krewlyzer/runall/main'
include { KREWLYZER_UXM } from './modules/local/krewlyzer/uxm/main'
include { KREWLYZER_EXTRACT } from './modules/local/krewlyzer/extract/main'
include { KREWLYZER_FSC } from './modules/local/krewlyzer/fsc/main'
include { KREWLYZER_FSR } from './modules/local/krewlyzer/fsr/main'
include { KREWLYZER_WPS } from './modules/local/krewlyzer/wps/main'
include { KREWLYZER_OCF } from './modules/local/krewlyzer/ocf/main'
include { KREWLYZER_FSD } from './modules/local/krewlyzer/fsd/main'
include { KREWLYZER_MOTIF } from './modules/local/krewlyzer/motif/main'
include { KREWLYZER_MFSD } from './modules/local/krewlyzer/mfsd/main'
include { FILTER_MAF } from './modules/local/krewlyzer/filter_maf/main'

// Function to auto-detect index
def get_index = { bam ->
    def bai = file("${bam}.bai")
    if (!bai.exists()) bai = file("${bam.toString().replace('.bam', '.bai')}")
    return bai
}

workflow {
    
    // Help Message
    if (params.help) {
        log.info """
        ================================================================
         K R E W L Y Z E R   P I P E L I N E  (v${workflow.manifest.version})
        ================================================================
         Usage:
         nextflow run main.nf --samplesheet samples.csv --ref hg19.fa [options]

         Input:
         --samplesheet     CSV with columns: sample, bam, meth_bam, vcf, bed, maf, single_sample_maf
         --ref             Reference genome FASTA

         Options:
         --outdir          Output directory (default: ./results)
         --targets         Target regions BED (optional)
         --mapq            Min MAPQ (default: 20)
         --threads         Threads per process (default: 8)
         --minlen          Min fragment length (default: 65)
         --maxlen          Max fragment length (default: 400)
         --skip_duplicates Skip duplicates (default: true)
         
         Output Format:
         --generate_json   Generate unified sample.features.json for ML (default: false)
         --output_format   Output format: auto, tsv, parquet, json (default: auto)
         
         Profiles:
         -profile docker   Run with Docker
         -profile slurm    Run on Slurm cluster (cmobic_cpu)
        ================================================================
        """.stripIndent()
        exit 0
    }

    // =====================================================
    // ASSAY RESOLUTION FUNCTIONS (XS1, XS2, WGS)
    // =====================================================
    
    // Assay code to PON filename mapping
    def assayToPon = [
        'XS1': 'msk-access-v1.pon.parquet',
        'XS2': 'msk-access-v2.pon.parquet',
        'WGS': 'wgs.pon.parquet'
    ]
    
    // Assay code to targets filename mapping
    def assayToTargets = [
        'XS1': 'XS1_targets.bed',
        'XS2': 'XS2_targets.bed'
        // WGS has no targets
    ]

    // Resolve PON model: explicit > assay-based > global param
    def resolvePon = { row ->
        if (row.pon) return file(row.pon)
        if (row.assay && params.asset_dir && assayToPon[row.assay]) {
            def assay_pon = file("${params.asset_dir}/pon/${assayToPon[row.assay]}")
            if (assay_pon.exists()) return assay_pon
            log.warn "PON not found for assay ${row.assay}: ${assay_pon}"
        }
        if (params.pon_model) return file(params.pon_model)
        return []
    }

    // Resolve targets: explicit > assay-based > global param
    def resolveTargets = { row ->
        if (row.targets) return file(row.targets)
        if (row.assay && params.asset_dir && assayToTargets[row.assay]) {
            def assay_targets = file("${params.asset_dir}/targets/${assayToTargets[row.assay]}")
            if (assay_targets.exists()) return assay_targets
            log.warn "Targets not found for assay ${row.assay}: ${assay_targets}"
        }
        if (params.targets) return file(params.targets)
        return []
    }

    // 1. Parse Sample Sheet
    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row ->
            def meta = [id: row.sample, assay: row.assay ?: 'UNKNOWN']
            
            // Debug logging
            if (params.verbose) {
                log.info "Processing ${row.sample}: ASSAY=${row.assay ?: 'N/A'}, BAM=${row.bam}, BED=${row.bed}"
            }

            // Auto-detect Index for WGS BAM
            def bam = row.bam ? file(row.bam) : null
            def bai = (bam && get_index(bam).exists()) ? get_index(bam) : []
            
            // Auto-detect Index for Meth BAM
            def mbam = row.meth_bam ? file(row.meth_bam) : null
            def mbai = (mbam && get_index(mbam).exists()) ? get_index(mbam) : []
            
            def bed = row.bed ? file(row.bed) : null
            def vcf = row.vcf ? file(row.vcf) : []
            def maf = row.maf ? file(row.maf) : null
            // single_sample_maf: if true, skip filtering and use MAF directly
            def single_sample = row.single_sample_maf?.toLowerCase() in ['true', 'yes', '1']
            
            // Resolve PON and targets based on assay
            def pon = resolvePon(row)
            def targets = resolveTargets(row)

            [ meta, bam, bai, vcf, mbam, mbai, bed, maf, single_sample, pon, targets ]
        }
        .multiMap { meta, bam, bai, vcf, mbam, mbai, bed, maf, single_sample, pon, targets ->
            runall: bam  ? [ meta, bam, bai, vcf, pon, targets ] : null
            methyl: mbam ? [ meta, mbam, mbai ]    : null
            bedops: bed  ? [ meta, bed, pon, targets ]           : null
            // Multi-sample MAF: needs filtering
            maf_multi:  (bam && maf && !single_sample) ? [ meta, bam, bai, maf ] : null
            // Single-sample MAF: skip filtering, pass directly
            maf_single: (bam && maf && single_sample)  ? [ meta, bam, bai, maf ] : null
        }
        .set { ch_inputs }

    // Filter Channels (Fix multiMap nulls)
    def ch_runall     = ch_inputs.runall.filter { it }
    def ch_methyl     = ch_inputs.methyl.filter { it }
    def ch_bedops     = ch_inputs.bedops.filter { it }
    def ch_maf_multi  = ch_inputs.maf_multi.filter { it }
    def ch_maf_single = ch_inputs.maf_single.filter { it }

    // 2. Run-All (Optimal path for WGS BAMs)
    KREWLYZER_RUNALL(
        ch_runall,
        file(params.ref)
    )

    // 3. Methylation (UXM path)
    KREWLYZER_UXM(
        ch_methyl,
        file(params.ref)
    )

    // 4. Bed Operations (Lightweight path for pre-extracted BEDs)
    // Extract bed-only channel (meta, bed) for simpler tools
    ch_bed_only = ch_bedops.map { meta, bed, pon, targets -> [ meta, bed ] }
    
    // FSC/FSR get optional targets from channel
    ch_bed_targets = ch_bedops.map { meta, bed, pon, targets -> [ meta, bed, targets ?: [] ] }
    KREWLYZER_FSC(ch_bed_targets.map { meta, bed, targets -> [ meta, bed ] }, ch_bed_targets.map { it[2] })
    KREWLYZER_FSR(ch_bed_targets.map { meta, bed, targets -> [ meta, bed ] }, ch_bed_targets.map { it[2] })
    
    // WPS with optional anchors
    KREWLYZER_WPS(
        ch_bed_only,
        file(params.ref),
        params.wps_anchors ? file(params.wps_anchors) : [],
        params.wps_background ? file(params.wps_background) : []
    )
    
    // OCF
    KREWLYZER_OCF(ch_bed_only, file(params.ref))
    
    // FSD - only needs bed
    KREWLYZER_FSD(ch_bed_only)

    // 5. MFSD (Mutant Fragment Size Distribution)
    
    // 5a. Multi-sample MAF: Filter by Tumor_Sample_Barcode, then run mfsd
    ch_maf_to_filter = ch_maf_multi.map { meta, bam, bai, maf -> [ meta, maf ] }
    FILTER_MAF(ch_maf_to_filter)
    
    // Filter out MAFs with zero variants (only header line = no data)
    ch_filtered_valid = FILTER_MAF.out.maf
        .filter { meta, filtered_maf ->
            // Efficiently check for data rows without loading entire file
            // Stream through file, stop as soon as we find 2+ non-comment lines
            def dataLineCount = 0
            filtered_maf.withReader { reader ->
                String line
                while ((line = reader.readLine()) != null && dataLineCount < 2) {
                    if (!line.startsWith('#') && line.trim()) {
                        dataLineCount++
                    }
                }
            }
            if (dataLineCount <= 1) {  // Only header, no data rows
                log.warn "Sample ${meta.id}: No variants after MAF filtering - skipping MFSD"
                return false
            }
            return true
        }
    
    // Join filtered MAF back with BAM info
    ch_mfsd_filtered = ch_maf_multi
        .map { meta, bam, bai, maf -> [ meta.id, meta, bam, bai ] }
        .join(ch_filtered_valid.map { meta, filtered_maf -> [ meta.id, filtered_maf ] })
        .map { id, meta, bam, bai, filtered_maf -> [ meta, bam, bai, filtered_maf ] }
    
    // 5b. Single-sample MAF: Pass directly to mfsd (no filtering)
    ch_mfsd_direct = ch_maf_single
    
    // Combine and run MFSD
    ch_mfsd_all = ch_mfsd_filtered.mix(ch_mfsd_direct)
    KREWLYZER_MFSD(ch_mfsd_all)

}
