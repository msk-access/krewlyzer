#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import Modules
include { KREWLYZER_UXM } from './modules/local/krewlyzer/uxm/main'
include { KREWLYZER_EXTRACT } from './modules/local/krewlyzer/extract/main'
include { KREWLYZER_FSC } from './modules/local/krewlyzer/fsc/main'
include { KREWLYZER_FSR } from './modules/local/krewlyzer/fsr/main'
include { KREWLYZER_WPS } from './modules/local/krewlyzer/wps/main'
include { KREWLYZER_OCF } from './modules/local/krewlyzer/ocf/main'
include { KREWLYZER_FSD } from './modules/local/krewlyzer/fsd/main'
include { KREWLYZER_MOTIF } from './modules/local/krewlyzer/motif/main'
include { KREWLYZER_MFSD } from './modules/local/krewlyzer/mfsd/main'
include { KREWLYZER_REGION_ENTROPY } from './modules/local/krewlyzer/region_entropy/main'
include { KREWLYZER_REGION_MDS } from './modules/local/krewlyzer/region_mds/main'
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
         --samplesheet     CSV with columns: sample, bam, meth_bam, vcf, bed, maf, single_sample_maf, assay, pon, targets
         --ref             Reference genome FASTA

         Options:
         --outdir          Output directory (default: ./results)
         --targets         Target regions BED (optional)
         --mapq            Min MAPQ (default: 20)
         --threads         Threads per process (default: 8)
         --minlen          Min fragment length (default: 65)
         --maxlen          Max fragment length (default: 1000, extended FSD range)
         --skip_duplicates Skip duplicates (default: true)
         
         Region Entropy:
         --no_tfbs         Disable TFBS entropy analysis
         --no_atac         Disable ATAC entropy analysis
         --skip_pon        Skip PON z-score normalization
         
         Region MDS (Per-Gene Motif Diversity):
         Enabled automatically when 'assay' column in samplesheet is set (XS1, XS2, WGS).
         Uses bundled gene BED files for each assay.
         
         Assay Codes (samplesheet 'assay' column):
         XS1               MSK-ACCESS v1 (128 genes)
         XS2               MSK-ACCESS v2 (146 genes)
         WGS               Whole Genome Sequencing
         
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
    
    // Assay code to PON filename mapping (matches data/pon/GRCh37/)
    def assayToPon = [
        'XS1': 'GRCh37/xs1.pon.parquet',
        'XS2': 'GRCh37/xs2.pon.parquet',
        'WGS': 'GRCh37/wgs.pon.parquet'
    ]
    
    // Assay code to targets filename mapping (matches data/targets/GRCh37/)
    def assayToTargets = [
        'XS1': 'GRCh37/xs1.targets.bed',
        'XS2': 'GRCh37/xs2.targets.bed'
        // WGS has no targets
    ]

    // Assay code to WPS anchors filename mapping (matches data/WpsAnchors/GRCh37/)
    def assayToWpsAnchors = [
        'XS1': 'GRCh37/xs1.wps_anchors.bed.gz',
        'XS2': 'GRCh37/xs2.wps_anchors.bed.gz'
        // WGS uses genome-wide anchors (default)
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

    // Resolve WPS anchors: assay-based > global param > default
    def resolveWpsAnchors = { row ->
        if (params.wps_anchors) return file(params.wps_anchors)
        if (row.assay && params.asset_dir && assayToWpsAnchors[row.assay]) {
            def assay_anchors = file("${params.asset_dir}/WpsAnchors/${assayToWpsAnchors[row.assay]}")
            if (assay_anchors.exists()) return assay_anchors
        }
        return []
    }

    // =====================================================
    // 1. PARSE SAMPLE SHEET
    // =====================================================
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
            // BAM path: needs extraction first
            extract: bam ? [ meta, bam, bai, pon, targets ] : null
            // Pre-extracted BED path: skip extraction
            bedops: bed ? [ meta, bed, pon, targets ] : null
            // Methylation path
            methyl: mbam ? [ meta, mbam, mbai ] : null
            // MFSD paths
            maf_multi:  (bam && maf && !single_sample) ? [ meta, bam, bai, maf ] : null
            maf_single: (bam && maf && single_sample)  ? [ meta, bam, bai, maf ] : null
        }
        .set { ch_inputs }

    // Filter Channels (Fix multiMap nulls)
    ch_extract    = ch_inputs.extract.filter { it }
    ch_bedops     = ch_inputs.bedops.filter { it }
    ch_methyl     = ch_inputs.methyl.filter { it }
    ch_maf_multi  = ch_inputs.maf_multi.filter { it }
    ch_maf_single = ch_inputs.maf_single.filter { it }

    // =====================================================
    // 2. EXTRACTION (BAM → BED.gz)
    // =====================================================
    
    // Extract fragments from BAM files
    ch_extract_input = ch_extract.map { meta, bam, bai, pon, targets -> 
        [ meta, bam, bai, targets ?: [] ] 
    }
    KREWLYZER_EXTRACT(
        ch_extract_input.map { meta, bam, bai, targets -> [ meta, bam, bai ] },
        file(params.ref),
        ch_extract_input.map { it[3] }  // targets
    )
    
    // Combine extracted BEDs with PON/targets for downstream
    ch_extracted_with_meta = KREWLYZER_EXTRACT.out.bed
        .join(ch_extract.map { meta, bam, bai, pon, targets -> [ meta, pon, targets ] })
        .map { meta, bed, pon, targets -> [ meta, bed, pon ?: [], targets ?: [] ] }
    
    // Merge extracted BEDs with pre-extracted BEDs
    ch_all_beds = ch_extracted_with_meta.mix(ch_bedops)

    // =====================================================
    // 3. PARALLEL FEATURE EXTRACTION (All run in parallel!)
    // =====================================================
    
    // Convenience channels
    ch_bed_only = ch_all_beds.map { meta, bed, pon, targets -> [ meta, bed ] }
    ch_bed_targets = ch_all_beds.map { meta, bed, pon, targets -> [ meta, bed, targets ?: [] ] }
    
    // FSC - Fragment Size Coverage
    KREWLYZER_FSC(
        ch_bed_targets.map { meta, bed, targets -> [ meta, bed ] }, 
        ch_bed_targets.map { it[2] }
    )
    
    // FSR - Fragment Size Ratio (depends on FSC conceptually, but runs in parallel)
    KREWLYZER_FSR(
        ch_bed_targets.map { meta, bed, targets -> [ meta, bed ] }, 
        ch_bed_targets.map { it[2] }
    )
    
    // FSD - Fragment Size Distribution
    KREWLYZER_FSD(ch_bed_only)
    
    // WPS - Windowed Protection Score
    KREWLYZER_WPS(
        ch_bed_only,
        file(params.ref),
        params.wps_anchors ? file(params.wps_anchors) : [],
        params.wps_background ? file(params.wps_background) : []
    )
    
    // OCF - Orientation-aware cfDNA Fragmentation
    KREWLYZER_OCF(ch_bed_only, file(params.ref))
    
    // MOTIF - End motif and MDS
    KREWLYZER_MOTIF(
        ch_bed_only,
        file(params.ref)
    )
    
    // REGION_ENTROPY - TFBS/ATAC size entropy
    ch_bed_pon_targets = ch_all_beds.map { meta, bed, pon, targets -> 
        [ meta, bed, pon ?: [], targets ?: [] ] 
    }
    KREWLYZER_REGION_ENTROPY(
        ch_bed_pon_targets.map { meta, bed, pon, targets -> [ meta, bed ] },
        file(params.ref),
        ch_bed_pon_targets.map { it[2] },  // pon
        ch_bed_pon_targets.map { it[3] }   // targets
    )

    // REGION_MDS - Per-gene/exon Motif Diversity Score (BAM-based, needs assay)
    ch_bam_for_mds = ch_extract.map { meta, bam, bai, pon, targets -> 
        [ meta, bam, bai, pon ?: [], meta.assay ?: '' ] 
    }
    KREWLYZER_REGION_MDS(
        ch_bam_for_mds.map { meta, bam, bai, pon, assay -> [ meta, bam, bai ] },
        file(params.ref),
        ch_bam_for_mds.map { it[3] },  // pon
        ch_bam_for_mds.map { it[4] }   // assay
    )

    // =====================================================
    // 4. METHYLATION (UXM) - Parallel path
    // =====================================================
    KREWLYZER_UXM(
        ch_methyl,
        file(params.ref)
    )

    // =====================================================
    // 5. MFSD (Mutant Fragment Size Distribution) - Parallel path
    // =====================================================
    
    // 5a. Multi-sample MAF: Filter by Tumor_Sample_Barcode, then run mfsd
    ch_maf_to_filter = ch_maf_multi.map { meta, bam, bai, maf -> [ meta, maf ] }
    FILTER_MAF(ch_maf_to_filter)
    
    // Filter out MAFs with zero variants (only header line = no data)
    ch_filtered_valid = FILTER_MAF.out.maf
        .filter { meta, filtered_maf ->
            // Efficiently check for data rows without loading entire file
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
