#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KREWLYZER PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fragmentomics feature extraction for liquid biopsy cfDNA analysis.
    
    Modes:
    - use_runall = true:  Unified run-all command per sample (default, recommended)
    - use_runall = false: Individual tool-level modules (legacy)
    
    https://github.com/msk-access/krewlyzer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { KREWLYZER } from './workflows/krewlyzer/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
    log.info """
    ================================================================
     K R E W L Y Z E R   P I P E L I N E  (v${workflow.manifest.version ?: '0.5.0'})
    ================================================================
     Usage:
     nextflow run main.nf --samplesheet samples.csv --ref hg19.fa [options]

     Required:
     --samplesheet     CSV with columns: sample, bam, meth_bam, vcf, bed, maf, single_sample_maf, assay, pon, targets
     --ref             Reference genome FASTA

     Mode Selection:
     --use_runall      Run unified run-all command (default: true)
                       Set to false for tool-level individual modules

     Options:
     --outdir          Output directory (default: ./results)
     --genome          Genome build: hg19, hg38 (default: hg19)
     --mapq            Min MAPQ (default: 20)
     --minlen          Min fragment length (default: 65)
     --maxlen          Max fragment length (default: 1000)
     --generate_json   Generate unified sample.features.json (default: true)
     
     Panel Mode:
     --skip_pon        Skip PON z-score normalization
     --skip_target_regions  Force WGS mode (ignore bundled targets)
     --bait_padding    Bait edge padding in bp (default: 50)
     
     Feature Control:
     --no_tfbs         Disable TFBS entropy analysis
     --no_atac         Disable ATAC entropy analysis
     --disable_e1_aggregation  Skip E1-only FSC aggregation
     --region_mds_e1_only      Run region-MDS on E1 only

     Assay Codes (samplesheet 'assay' column):
     XS1               MSK-ACCESS v1 (128 genes)
     XS2               MSK-ACCESS v2 (146 genes)
     WGS               Whole Genome Sequencing
     
     Profiles:
     -profile docker   Run with Docker
     -profile slurm    Run on Slurm cluster
    ================================================================
    """.stripIndent()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_KREWLYZER {
    main:
    
    // Check required parameters
    if (!params.samplesheet) {
        log.error "ERROR: --samplesheet is required"
        helpMessage()
        System.exit(1)
    }
    if (!params.ref) {
        log.error "ERROR: --ref is required"
        helpMessage()
        System.exit(1)
    }

    ch_samplesheet = Channel.value(file(params.samplesheet, checkIfExists: true))
    ch_fasta = file(params.ref, checkIfExists: true)

    KREWLYZER(
        ch_samplesheet,
        ch_fasta
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENTRY POINT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Show help if requested
    if (params.help) {
        helpMessage()
        System.exit(0)
    }
    
    NFCORE_KREWLYZER()
}
