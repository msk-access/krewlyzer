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
         --samplesheet     CSV with columns: sample, bam, meth_bam, vcf, bed
         --ref             Reference genome FASTA

         Options:
         --outdir          Output directory (default: ./results)
         --targets         Target regions BED (optional)
         --mapq            Min MAPQ (default: 20)
         --threads         Threads per process (default: 8)
         --minlen          Min fragment length (default: 65)
         --maxlen          Max fragment length (default: 400)
         --skip_duplicates Skip duplicates (default: true)
         
         Profiles:
         -profile docker   Run with Docker
         -profile slurm    Run on Slurm cluster (cmobic_cpu)
        ================================================================
        """.stripIndent()
        exit 0
    }

    // 1. Parse Sample Sheet
    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row ->
            def meta = [id: row.sample]
            
            // Debug logging
            // log.info "Processing ${row.sample}: BAM=${row.bam}, METH=${row.meth_bam}, BED=${row.bed}, VCF=${row.vcf}"

            // Auto-detect Index for WGS BAM
            def bam = row.bam ? file(row.bam) : null
            def bai = (bam && get_index(bam).exists()) ? get_index(bam) : []
            
            // Auto-detect Index for Meth BAM
            def mbam = row.meth_bam ? file(row.meth_bam) : null
            def mbai = (mbam && get_index(mbam).exists()) ? get_index(mbam) : []
            
            def bed = row.bed ? file(row.bed) : null
            def vcf = row.vcf ? file(row.vcf) : []

            [ meta, bam, bai, vcf, mbam, mbai, bed ]
        }
        .multiMap { meta, bam, bai, vcf, mbam, mbai, bed ->
            runall: bam  ? [ meta, bam, bai, vcf ] : null
            methyl: mbam ? [ meta, mbam, mbai ]    : null
            bedops: bed  ? [ meta, bed ]           : null
        }
        .set { ch_inputs }

    // Filter Channels (Fix multiMap nulls)
    def ch_runall = ch_inputs.runall.filter { it }
    def ch_methyl = ch_inputs.methyl.filter { it }
    def ch_bedops = ch_inputs.bedops.filter { it }

    // 2. Run-All (Optimal path for WGS BAMs)
    KREWLYZER_RUNALL(
        ch_runall,
        file(params.ref),
        params.targets ? file(params.targets) : []
    )

    // 3. Methylation (UXM path)
    KREWLYZER_UXM(
        ch_methyl,
        file(params.ref)
    )

    // 4. Bed Operations (Lightweight path for pre-extracted BEDs)
    // Run compatible tools in parallel
    
    // Tools that support optional targets
    KREWLYZER_FSC(ch_bedops, params.targets ? file(params.targets) : [])
    KREWLYZER_FSR(ch_bedops, params.targets ? file(params.targets) : [])
    
    // Tools that require Reference
    KREWLYZER_WPS(ch_bedops, file(params.ref))
    KREWLYZER_OCF(ch_bedops, file(params.ref))
    
    // Tools that need only BED
    KREWLYZER_FSD(ch_bedops)

}
