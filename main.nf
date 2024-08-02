#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/raredisease
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/raredisease
    Website: https://nf-co.re/raredisease
    Slack  : https://nfcore.slack.com/channels/raredisease
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta                          = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fai                            = WorkflowMain.getGenomeAttribute(params, 'fai')
params.bwa                            = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2                        = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.call_interval                  = WorkflowMain.getGenomeAttribute(params, 'call_interval')
params.cadd_resources                 = WorkflowMain.getGenomeAttribute(params, 'cadd_resources')
params.gcnvcaller_model               = WorkflowMain.getGenomeAttribute(params, 'gcnvcaller_model')
params.gens_interval_list             = WorkflowMain.getGenomeAttribute(params, 'gens_interval_list')
params.gens_pon                       = WorkflowMain.getGenomeAttribute(params, 'gens_pon')
params.gens_gnomad_pos                = WorkflowMain.getGenomeAttribute(params, 'gens_gnomad_pos')
params.gnomad_af                      = WorkflowMain.getGenomeAttribute(params, 'gnomad_af')
params.gnomad_af_idx                  = WorkflowMain.getGenomeAttribute(params, 'gnomad_af_idx')
params.intervals_wgs                  = WorkflowMain.getGenomeAttribute(params, 'intervals_wgs')
params.intervals_y                    = WorkflowMain.getGenomeAttribute(params, 'intervals_y')
params.dbsnp                          = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi                      = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.known_indels                   = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi               = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')
params.known_indels_2                 = WorkflowMain.getGenomeAttribute(params, 'known_indels_2')
params.known_indels_2_tbi             = WorkflowMain.getGenomeAttribute(params, 'known_indels_2_tbi')
params.ml_model                       = WorkflowMain.getGenomeAttribute(params, 'ml_model')
params.mt_fasta                       = WorkflowMain.getGenomeAttribute(params, 'mt_fasta')
params.ploidy_model                   = WorkflowMain.getGenomeAttribute(params, 'ploidy_model')
params.reduced_penetrance             = WorkflowMain.getGenomeAttribute(params, 'reduced_penetrance')
params.readcount_intervals            = WorkflowMain.getGenomeAttribute(params, 'readcount_intervals')
params.sequence_dictionary            = WorkflowMain.getGenomeAttribute(params, 'sequence_dictionary')
params.score_config_snv               = WorkflowMain.getGenomeAttribute(params, 'score_config_snv')
params.score_config_sv                = WorkflowMain.getGenomeAttribute(params, 'score_config_sv')
params.svdb_query_dbs                 = WorkflowMain.getGenomeAttribute(params, 'svdb_query_dbs')
params.target_bed                     = WorkflowMain.getGenomeAttribute(params, 'target_bed')
params.variant_catalog                = WorkflowMain.getGenomeAttribute(params, 'variant_catalog')
params.vep_filters                    = WorkflowMain.getGenomeAttribute(params, 'vep_filters')
params.vcfanno_resources              = WorkflowMain.getGenomeAttribute(params, 'vcfanno_resources')
params.vcfanno_toml                   = WorkflowMain.getGenomeAttribute(params, 'vcfanno_toml')
params.vcfanno_lua                    = WorkflowMain.getGenomeAttribute(params, 'vcfanno_lua')
params.vep_cache                      = WorkflowMain.getGenomeAttribute(params, 'vep_cache')
params.vep_cache_version              = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.dict                           = WorkflowMain.getGenomeAttribute(params, 'dict')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RAREDISEASE             } from './workflows/raredisease'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_raredisease_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_RAREDISEASE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    RAREDISEASE (
        samplesheet
    )

    emit:
    multiqc_report = RAREDISEASE.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RAREDISEASE (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RAREDISEASE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}
