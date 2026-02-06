/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT_CHECK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Samplesheet validation and asset resolution for krewlyzer pipeline.
    
    Features:
    - Parse samplesheet CSV
    - Auto-detect BAM indexes
    - Resolve PON, targets, WPS anchors from assay codes
    - Handle single_sample_maf bypass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Auto-detect BAM index
def get_index = { bam ->
    def bai = file("${bam}.bai")
    if (!bai.exists()) bai = file("${bam.toString().replace('.bam', '.bai')}")
    return bai
}

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    // =====================================================
    // ASSAY RESOLUTION MAPS
    // =====================================================
    def assayToPon = [
        'XS1': 'GRCh37/xs1.pon.parquet',
        'XS2': 'GRCh37/xs2.pon.parquet',
        'WGS': 'GRCh37/wgs.pon.parquet'
    ]
    
    def assayToTargets = [
        'XS1': 'GRCh37/xs1.targets.bed',
        'XS2': 'GRCh37/xs2.targets.bed'
    ]

    def assayToWpsAnchors = [
        'XS1': 'GRCh37/xs1.wps_anchors.bed.gz',
        'XS2': 'GRCh37/xs2.wps_anchors.bed.gz'
    ]

    // =====================================================
    // RESOLUTION FUNCTIONS
    // =====================================================
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

    def resolveWpsAnchors = { row ->
        if (params.wps_anchors) return file(params.wps_anchors)
        if (row.assay && params.asset_dir && assayToWpsAnchors[row.assay]) {
            def assay_anchors = file("${params.asset_dir}/WpsAnchors/${assayToWpsAnchors[row.assay]}")
            if (assay_anchors.exists()) return assay_anchors
        }
        return []
    }

    // =====================================================
    // PARSE SAMPLESHEET
    // =====================================================
    samplesheet
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, assay: row.assay ?: 'UNKNOWN']
            
            // Debug logging
            if (params.verbose) {
                log.info "Processing ${row.sample}: ASSAY=${row.assay ?: 'N/A'}, BAM=${row.bam}, BED=${row.bed}"
            }

            // Auto-detect BAM indexes
            def bam = row.bam ? file(row.bam) : null
            def bai = (bam && get_index(bam).exists()) ? get_index(bam) : []
            
            // mFSD-specific BAM (duplex consensus for mFSD, falls back to main BAM)
            def mfsd_bam = row.mfsd_bam ? file(row.mfsd_bam) : null
            def mfsd_bai = (mfsd_bam && get_index(mfsd_bam).exists()) ? get_index(mfsd_bam) : []
            
            def mbam = row.meth_bam ? file(row.meth_bam) : null
            def mbai = (mbam && get_index(mbam).exists()) ? get_index(mbam) : []
            
            def bed = row.bed ? file(row.bed) : null
            def vcf = row.vcf ? file(row.vcf) : []
            def maf = row.maf ? file(row.maf) : null
            def single_sample = row.single_sample_maf?.toLowerCase() in ['true', 'yes', '1']
            
            // Resolve assets
            def pon = resolvePon(row)
            def targets = resolveTargets(row)
            def wps_anchors = resolveWpsAnchors(row)
            def wps_background = params.wps_background ? file(params.wps_background) : []

            [meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single_sample, pon, targets, wps_anchors, wps_background]
        }
        .set { ch_parsed }

    // =====================================================
    // BRANCH INTO CHANNELS
    // =====================================================
    ch_parsed
        .branch {
            // Samples with BAM for run-all or tool_level
            bam_samples: it[1] != null
            // Samples with pre-extracted BED only (tool_level only)
            bed_samples: it[5] != null
        }
        .set { ch_branched }

    // Format for run-all: [meta, bam, bai, mfsd_bam, mfsd_bai, bisulfite_bam, variants, pon, targets, wps_anchors, wps_background]
    ch_runall = ch_branched.bam_samples.map { 
        meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
        // Use MAF if present, else VCF, else empty
        def variants = maf ?: (vcf ?: [])
        [meta, bam, bai, mfsd_bam ?: [], mfsd_bai ?: [], mbam ?: [], variants, pon ?: [], targets ?: [], wps_anchors ?: [], wps_bg ?: []]
    }

    // For tool_level: extract path
    ch_extract = ch_branched.bam_samples.map {
        meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
        [meta, bam, bai, pon ?: [], targets ?: []]
    }

    // Pre-extracted BEDs (tool_level only)
    ch_beds = ch_branched.bed_samples.map {
        meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
        [meta, bed, pon ?: [], targets ?: []]
    }

    // Methylation samples
    ch_methyl = ch_parsed
        .filter { it[5] != null }  // mbam exists (index shifted by 2 for mfsd_bam)
        .map { 
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
            [meta, mbam, mbai]
        }

    // MAF samples for filtering (multi-sample, not single_sample_maf)
    ch_maf_multi = ch_parsed
        .filter { it[9] != null && !it[10] }  // maf exists AND not single_sample (indices shifted)
        .map {
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
            [meta, bam, bai, maf]
        }

    // MAF samples bypass (single_sample_maf = true)
    ch_maf_single = ch_parsed
        .filter { it[9] != null && it[10] }  // maf exists AND single_sample (indices shifted)
        .map {
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
            [meta, bam, bai, maf]
        }

    emit:
    runall     = ch_runall      // For run-all module
    extract    = ch_extract     // For tool_level extract
    beds       = ch_beds        // Pre-extracted BEDs (tool_level only)
    methyl     = ch_methyl      // Methylation samples
    maf_multi  = ch_maf_multi   // MAFs needing filtering
    maf_single = ch_maf_single  // MAFs bypassing filter
}
