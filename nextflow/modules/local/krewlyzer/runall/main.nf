/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER_RUNALL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run all krewlyzer fragmentomics feature extraction tools on a single BAM.
    
    Outputs: Extract, Motif, FSC, FSR, FSD, WPS, OCF, Region Entropy, Region MDS
    Optional: mFSD (with --variants), UXM (with --bisulfite-bam)
    
    Panel mode: Use meta.assay (xs1/xs2) for gene-centric FSC and dual WPS output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KREWLYZER_RUNALL {
    tag "$meta.id"
    label 'process_high'

    container "ghcr.io/msk-access/krewlyzer:0.5.3"

    input:
    tuple val(meta), path(bam), path(bai), path(mfsd_bam), path(mfsd_bai), path(bisulfite_bam), path(variants), path(pon), path(targets), path(wps_anchors), path(wps_background)
    path fasta

    output:
    // Core outputs
    tuple val(meta), path("*.bed.gz"),               emit: bed
    tuple val(meta), path("*.bed.gz.tbi"),           emit: bed_index, optional: true
    tuple val(meta), path("*.metadata.json"),        emit: metadata, optional: true
    tuple val(meta), path("*.features.json"),        emit: features_json, optional: true
    
    // FSC outputs
    tuple val(meta), path("*.FSC.tsv"),              emit: fsc, optional: true
    tuple val(meta), path("*.FSC.ontarget.tsv"),     emit: fsc_ontarget, optional: true
    tuple val(meta), path("*.FSC.gene.tsv"),         emit: fsc_gene, optional: true
    tuple val(meta), path("*.FSC.regions.tsv"),      emit: fsc_regions, optional: true
    tuple val(meta), path("*.FSC.regions.e1only.tsv"), emit: fsc_e1, optional: true
    
    // FSR outputs
    tuple val(meta), path("*.FSR.tsv"),              emit: fsr, optional: true
    tuple val(meta), path("*.FSR.ontarget.tsv"),     emit: fsr_ontarget, optional: true
    
    // FSD outputs
    tuple val(meta), path("*.FSD.tsv"),              emit: fsd, optional: true
    tuple val(meta), path("*.FSD.ontarget.tsv"),     emit: fsd_ontarget, optional: true
    
    // WPS outputs
    tuple val(meta), path("*.WPS.parquet"),          emit: wps, optional: true
    tuple val(meta), path("*.WPS.panel.parquet"),    emit: wps_panel, optional: true
    tuple val(meta), path("*.WPS_background.parquet"), emit: wps_background, optional: true
    
    // OCF outputs
    tuple val(meta), path("*.OCF.tsv"),              emit: ocf, optional: true
    tuple val(meta), path("*.OCF.sync.tsv"),         emit: ocf_sync, optional: true
    tuple val(meta), path("*.OCF.ontarget.tsv"),     emit: ocf_ontarget, optional: true
    tuple val(meta), path("*.OCF.ontarget.sync.tsv"), emit: ocf_ontarget_sync, optional: true
    tuple val(meta), path("*.OCF.offtarget.tsv"),    emit: ocf_offtarget, optional: true
    tuple val(meta), path("*.OCF.offtarget.sync.tsv"), emit: ocf_offtarget_sync, optional: true
    
    // Motif outputs
    tuple val(meta), path("*.EndMotif.tsv"),         emit: end_motif, optional: true
    tuple val(meta), path("*.EndMotif.ontarget.tsv"), emit: end_motif_ontarget, optional: true
    tuple val(meta), path("*.EndMotif1mer.tsv"),     emit: end_motif_1mer, optional: true
    tuple val(meta), path("*.BreakPointMotif.tsv"),  emit: bp_motif, optional: true
    tuple val(meta), path("*.BreakPointMotif.ontarget.tsv"), emit: bp_motif_ontarget, optional: true
    tuple val(meta), path("*.MDS.tsv"),              emit: mds, optional: true
    tuple val(meta), path("*.MDS.ontarget.tsv"),     emit: mds_ontarget, optional: true
    
    // TFBS/ATAC outputs
    tuple val(meta), path("*.TFBS.tsv"),             emit: tfbs, optional: true
    tuple val(meta), path("*.TFBS.sync.tsv"),        emit: tfbs_sync, optional: true
    tuple val(meta), path("*.TFBS.ontarget.tsv"),    emit: tfbs_ontarget, optional: true
    tuple val(meta), path("*.TFBS.ontarget.sync.tsv"), emit: tfbs_ontarget_sync, optional: true
    tuple val(meta), path("*.ATAC.tsv"),             emit: atac, optional: true
    tuple val(meta), path("*.ATAC.sync.tsv"),        emit: atac_sync, optional: true
    tuple val(meta), path("*.ATAC.ontarget.tsv"),    emit: atac_ontarget, optional: true
    tuple val(meta), path("*.ATAC.ontarget.sync.tsv"), emit: atac_ontarget_sync, optional: true
    
    // Region MDS outputs
    tuple val(meta), path("*.MDS.exon.tsv"),         emit: mds_exon, optional: true
    tuple val(meta), path("*.MDS.gene.tsv"),         emit: mds_gene, optional: true
    
    // mFSD outputs
    tuple val(meta), path("*.mFSD.tsv"),             emit: mfsd, optional: true
    tuple val(meta), path("*.mFSD.distributions.tsv"), emit: mfsd_dist, optional: true
    
    // UXM outputs
    tuple val(meta), path("*.UXM.tsv"),              emit: uxm, optional: true
    
    // GC correction outputs
    tuple val(meta), path("*.correction_factors.tsv"), emit: gc_factors, optional: true
    tuple val(meta), path("*.correction_factors.ontarget.tsv"), emit: gc_factors_ontarget, optional: true
    
    // Versions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Input file arguments
    def mfsd_bam_arg = mfsd_bam ? "--mfsd-bam ${mfsd_bam}" : ""
    def bisulfite_arg = bisulfite_bam ? "--bisulfite-bam ${bisulfite_bam}" : ""
    def variant_arg = variants ? "--variants ${variants}" : ""
    def targets_arg = targets ? "--target-regions ${targets}" : ""
    def pon_arg = pon ? "--pon-model ${pon}" : ""
    def wps_anchors_arg = wps_anchors ? "--wps-anchors ${wps_anchors}" : ""
    def wps_bg_arg = wps_background ? "--wps-background ${wps_background}" : ""
    
    // Meta-based arguments
    def assay_arg = meta.assay && meta.assay != 'UNKNOWN' && meta.assay != 'WGS' ? "--assay ${meta.assay.toLowerCase()}" : ""
    
    // Params-based arguments
    def genome_arg = params.genome ? "--genome ${params.genome}" : ""
    def mapq_arg = params.mapq != 20 ? "--mapq ${params.mapq}" : ""
    def minlen_arg = params.minlen != 65 ? "--minlen ${params.minlen}" : ""
    def maxlen_arg = params.maxlen != 1000 ? "--maxlen ${params.maxlen}" : ""
    def skip_dup_arg = params.skip_duplicates == false ? "--no-skip-duplicates" : ""
    def proper_pair_arg = params.require_proper_pair == false ? "--no-require-proper-pair" : ""
    def duplex_arg = params.duplex ? "--duplex" : ""
    def bait_padding_arg = params.bait_padding != 50 ? "--bait-padding ${params.bait_padding}" : ""
    
    // Skip/disable flags
    def skip_pon_arg = params.skip_pon ? "--skip-pon" : ""
    def pon_variant_arg = params.pon_variant != 'all_unique' ? "--pon-variant ${params.pon_variant}" : ""
    def skip_targets_arg = params.skip_target_regions ? "--skip-target-regions" : ""
    def no_tfbs_arg = params.no_tfbs ? "--no-tfbs" : ""
    def no_atac_arg = params.no_atac ? "--no-atac" : ""
    def e1_disable_arg = params.disable_e1_aggregation ? "--disable-e1-aggregation" : ""
    def e1_mds_arg = params.region_mds_e1_only ? "--region-mds-e1-only" : ""
    
    // Output flags
    def json_arg = params.generate_json ? "--generate-json" : ""
    def debug_arg = params.verbose ? "--debug" : ""

    """
    krewlyzer run-all \\
        -i $bam \\
        -r $fasta \\
        --output ./ \\
        --threads $task.cpus \\
        --sample-name $prefix \\
        $mfsd_bam_arg \\
        $bisulfite_arg \\
        $variant_arg \\
        $targets_arg \\
        $pon_arg \\
        $wps_anchors_arg \\
        $wps_bg_arg \\
        $assay_arg \\
        $genome_arg \\
        $mapq_arg \\
        $minlen_arg \\
        $maxlen_arg \\
        $skip_dup_arg \\
        $proper_pair_arg \\
        $duplex_arg \\
        $bait_padding_arg \\
        $skip_pon_arg \\
        $pon_variant_arg \\
        $skip_targets_arg \\
        $no_tfbs_arg \\
        $no_atac_arg \\
        $e1_disable_arg \\
        $e1_mds_arg \\
        $json_arg \\
        $debug_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: \$(krewlyzer --version | sed 's/krewlyzer //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi
    touch ${prefix}.FSC.tsv
    touch ${prefix}.FSR.tsv
    touch ${prefix}.FSD.tsv
    touch ${prefix}.WPS.parquet
    touch ${prefix}.WPS_background.parquet
    touch ${prefix}.OCF.tsv
    touch ${prefix}.EndMotif.tsv
    touch ${prefix}.BreakpointMotif.tsv
    touch ${prefix}.MDS.tsv
    touch ${prefix}.correction_factors.tsv
    echo '{"sample_id":"${prefix}","total_fragments":0}' > ${prefix}.metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krewlyzer: 0.5.3
    END_VERSIONS
    """
}
