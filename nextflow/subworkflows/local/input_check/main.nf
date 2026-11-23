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

// Auto-detect BAM index.
//
// A function, not a closure assigned to a variable. Nextflow 26 rejects the
// latter at the top level of a script -- "Statements cannot be mixed with
// script declarations" -- and the file will not compile. A function
// declaration IS a script declaration, so it is allowed in both 25 and 26.
//
// Safe to change because every call site invokes it as `get_index(bam)`; it is
// never passed around as a value, which is the one usage a function could not
// serve.
def get_index(bam) {
    def bai = file("${bam}.bai")
    if (!bai.exists()) bai = file("${bam.toString().replace('.bam', '.bai')}")
    return bai
}

// Asset resolution, as top-level FUNCTIONS.
//
// These were closures assigned inside the workflow body. Nextflow 26 parses
// that but cannot resolve the names from inside the nested `.map { }` that
// calls them -- "`resolvePon` is not defined". Placement was not the only
// thing 26 tightened; closure scoping was too.
//
// Each lookup table moves into the function that uses it. Each was used by
// exactly one, so nothing is duplicated -- and leaving them at the top level
// as `def x = [...]` would just be the statement-vs-declaration error again.
def resolvePon(row) {
    def assayToPon = [
        'XS1': 'GRCh37/xs1.pon.parquet',
        'XS2': 'GRCh37/xs2.pon.parquet',
        'WGS': 'GRCh37/wgs.pon.parquet'
    ]
    if (row.pon) return file(row.pon)
    if (row.assay && params.asset_dir && assayToPon[row.assay]) {
        def assay_pon = file("${params.asset_dir}/pon/${assayToPon[row.assay]}")
        if (assay_pon.exists()) return assay_pon
        log.warn "PON not found for assay ${row.assay}: ${assay_pon}"
    }
    if (params.pon_model) return file(params.pon_model)
    return []
}

def resolveTargets(row) {
    def assayToTargets = [
        'XS1': 'GRCh37/xs1.targets.bed',
        'XS2': 'GRCh37/xs2.targets.bed'
    ]
    if (row.targets) return file(row.targets)
    if (row.assay && params.asset_dir && assayToTargets[row.assay]) {
        def assay_targets = file("${params.asset_dir}/targets/${assayToTargets[row.assay]}")
        if (assay_targets.exists()) return assay_targets
        log.warn "Targets not found for assay ${row.assay}: ${assay_targets}"
    }
    if (params.targets) return file(params.targets)
    return []
}

def resolveWpsAnchors(row) {
    def assayToWpsAnchors = [
        'XS1': 'GRCh37/xs1.wps_anchors.bed.gz',
        'XS2': 'GRCh37/xs2.wps_anchors.bed.gz'
    ]
    if (params.wps_anchors) return file(params.wps_anchors)
    if (row.assay && params.asset_dir && assayToWpsAnchors[row.assay]) {
        def assay_anchors = file("${params.asset_dir}/WpsAnchors/${assayToWpsAnchors[row.assay]}")
        if (assay_anchors.exists()) return assay_anchors
    }
    return []
}


workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    // =====================================================
    // ASSAY RESOLUTION MAPS
    // =====================================================
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

            // mFSD reads one variant file. A MAF wins over a VCF, and since the
            // VCF is then ignored, say so once rather than leave no trace.
            if (maf && vcf) {
                log.warn "${row.sample}: both maf and vcf are set; mFSD uses the MAF and ignores the VCF"
            }
            
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

    // =====================================================
    // MULTIMAP: Send every bam_sample to ALL output channels
    // (Using separate .map() calls on a queue channel would
    //  round-robin items, splitting them across consumers)
    // =====================================================
    ch_branched.bam_samples
        .multiMap {
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->

            // tool-level: extract path
            extract:
                [meta, bam, bai, pon ?: [], targets ?: []]

            // run-all: meta carries ALL paths so join is unnecessary. The VCF
            // rides in meta too: the main workflow routes a VCF-only sample
            // around FILTER_MAF/SPLIT_MAF, which only read MAF.
            maf_for_filter:
                def new_meta = meta + [
                    single_sample: single ?: false,
                    _bam: bam, _bai: bai,
                    _mfsd_bam: mfsd_bam, _mfsd_bai: mfsd_bai,
                    _bisulfite: mbam,
                    _pon: pon, _targets: targets,
                    _wps_anchors: wps_anchors, _wps_bg: wps_bg
                ] + (vcf ? [_vcf: vcf] : [:])  // only when set: meta is hashed, so a key on every sample would void -resume for all of them
                [new_meta, maf ?: []]
        }
        .set { ch_bam_outputs }

    // Pre-extracted BEDs (tool_level only — uses bed_samples, not bam_samples)
    ch_beds = ch_branched.bed_samples.map {
        meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
        [meta, bed, pon ?: [], targets ?: []]
    }

    // Methylation samples (from ch_parsed, not affected by multiMap)
    ch_methyl = ch_parsed
        .filter { it[5] != null }  // mbam exists (index shifted by 2 for mfsd_bam)
        .map { 
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
            [meta, mbam, mbai]
        }

    // MAF samples for filtering (multi-sample, not single_sample_maf) — tool-level mode
    ch_maf_multi = ch_parsed
        .filter { it[9] != null && !it[10] }  // maf exists AND not single_sample
        .map {
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
            [meta, bam, bai, maf]
        }

    // MAF samples bypass (single_sample_maf = true) — tool-level mode
    ch_maf_single = ch_parsed
        .filter { it[9] != null && it[10] }  // maf exists AND single_sample
        .map {
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
            [meta, bam, bai, maf]
        }

    // VCF samples (no MAF) — tool-level mode. A VCF is not filtered by
    // sample: mFSD reads every record, as run-all does. A MAF wins when
    // both are set, hence `it[9] == null`. The BAM is picked as the MAF
    // channels pick theirs, so mFSD's input does not depend on the format.
    ch_vcf_only = ch_parsed
        .filter { it[8] && it[9] == null }  // vcf exists AND no maf
        .map {
            meta, bam, bai, mfsd_bam, mfsd_bai, mbam, mbai, bed, vcf, maf, single, pon, targets, wps_anchors, wps_bg ->
            [meta, bam, bai, vcf]
        }

    ch_extract        = ch_bam_outputs.extract
    ch_maf_for_filter = ch_bam_outputs.maf_for_filter

    emit:
    extract        = ch_extract         // For tool_level extract
    beds           = ch_beds            // Pre-extracted BEDs (tool_level only)
    methyl         = ch_methyl          // Methylation samples
    maf_multi      = ch_maf_multi       // MAFs needing filtering (tool-level mode)
    maf_single     = ch_maf_single      // MAFs bypassing filter (tool-level mode)
    vcf_only       = ch_vcf_only        // VCFs, no MAF (tool-level mode)
    maf_for_filter = ch_maf_for_filter  // ALL samples for run-all (MAF or VCF in meta)
}
