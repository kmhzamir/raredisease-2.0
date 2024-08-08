/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_raredisease_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRaredisease.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK MANDATORY PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def mandatoryParams = [
    "aligner",
    "fasta",
    "input",

]
def missingParamsCount = 0

for (param in mandatoryParams.unique()) {
    if (params[param] == null) {
        println("params." + param + " not set.")
        missingParamsCount += 1
    }
}

if (!params.skip_vep_filter) {
    if (!params.vep_filters && !params.vep_filters_scout_fmt) {
        println("params.vep_filters or params.vep_filters_scout_fmt should be set.")
        missingParamsCount += 1
    } else if (params.vep_filters && params.vep_filters_scout_fmt) {
        println("Either params.vep_filters or params.vep_filters_scout_fmt should be set.")
        missingParamsCount += 1
    }
}

if (missingParamsCount>0) {
    error("\nSet missing parameters and restart the run. For more information please check usage documentation on github.")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config              ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo                ? Channel.fromPath( params.multiqc_logo, checkIfExists: true )   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true)  : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOWS
//

include { CHECK_INPUT                           } from '../subworkflows/local/input_check'
include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { FASTP                                 } from '../modules/nf-core/fastp/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { ALIGN                                 } from '../subworkflows/local/align'
include { PROCESS_MT                            } from '../subworkflows/local/process_MT'
include { PREPARE_REFERENCES                    } from '../subworkflows/local/prepare_references'
include { SCATTER_GENOME                        } from '../subworkflows/local/scatter_genome'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { BAM_MARKDUPLICATES                    } from '../subworkflows/local/bam_markduplicates/main'
include { BAM_BASERECALIBRATOR                  } from '../subworkflows/local/bam_baserecalibrator/main'
include { BAM_APPLYBQSR                         } from '../subworkflows/local/bam_applybqsr/main'
include { QC_BAM                                } from '../subworkflows/local/qc_bam'
include { GATK4_HAPLOTYPECALLER                 } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { SMNCOPYNUMBERCALLER                   } from '../modules/nf-core/smncopynumbercaller/main'
include { FILTER_VEP as FILTER_VEP_SNV          } from '../modules/local/filter_vep'
include { CREATE_PEDIGREE_FILE                  } from '../modules/local/create_pedigree_file'
include { CALL_SNV                              } from '../subworkflows/local/call_snv'
include { GATK4_MERGEVCFS as MERGE_VCFS         } from '../modules/nf-core/gatk4/merge_vcfs/main'
include { CREATE_HGNCIDS_FILE                   } from '../modules/local/create_hgncids_file'
include { ANNOTATE_MT_SNVS                                   } from '../subworkflows/local/annotate_mt_snvs'
include { GENERATE_CLINICAL_SET as GENERATE_CLINICAL_SET_MT  } from '../subworkflows/local/generate_clinical_set'
include { RANK_VARIANTS as RANK_VARIANTS_MT                  } from '../subworkflows/local/rank_variants'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_MT                 } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_GENOME_SNVS                               } from '../subworkflows/local/annotate_genome_snvs'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV                } from '../subworkflows/local/annotate_consequence_pli'
include { GENERATE_CLINICAL_SET as GENERATE_CLINICAL_SET_SNV } from '../subworkflows/local/generate_clinical_set'
include { RANK_VARIANTS as RANK_VARIANTS_SNV                 } from '../subworkflows/local/rank_variants'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RAREDISEASE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_reports  = Channel.empty()

    ch_samples   = ch_samplesheet.map { meta, fastqs -> meta}
    ch_case_info = ch_samples.toList().map { CustomFunctions.createCaseChannel(it) }

    // Initialize file channels for PREPARE_REFERENCES subworkflow
    ch_genome_fasta             = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
    ch_genome_fai               = params.fai            ? Channel.fromPath(params.fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                        : Channel.empty()
    ch_target_bed_unprocessed   = params.target_bed     ? Channel.fromPath(params.target_bed).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]]) 
    known_indels                = params.known_indels   ? Channel.fromPath(params.known_indels).collect()            : Channel.value([])
    known_indels_2              = params.known_indels_2 ? Channel.fromPath(params.known_indels_2).collect()          : Channel.value([])
    dbsnp                       = params.dbsnp          ? Channel.fromPath(params.dbsnp).collect()                   : Channel.value([])
    ch_intervals_wgs            = params.intervals_wgs  ? Channel.fromPath(params.intervals_wgs).collect()
                                                        : Channel.empty()
    ch_intervals_y              = params.intervals_y    ? Channel.fromPath(params.intervals_y).collect()
                                                        : Channel.empty()
    chh_target_bed = params.target_bed ? Channel.fromPath(params.target_bed).map{ it -> [it] }.collect() : Channel.value([[],[]])
    ch_gnomad_af_tab            = params.gnomad_af      ? Channel.fromPath(params.gnomad_af).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_vep_cache_unprocessed    = params.vep_cache      ? Channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_mt_fasta                 = params.mt_fasta       ? Channel.fromPath(params.mt_fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.empty()

    // Prepare references and indices.
    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_genome_fai,
        ch_mt_fasta,
        ch_target_bed_unprocessed,
        dbsnp,
        known_indels,
        known_indels_2,
        ch_gnomad_af_tab,
        ch_vep_cache_unprocessed
    )
    .set { ch_references }

    // Gather built indices or get them from the params
    ch_cadd_header              = Channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_variant_consequences     = Channel.fromPath("$projectDir/assets/variant_consequences_v1.txt", checkIfExists: true).collect()
    ch_foundin_header           = Channel.fromPath("$projectDir/assets/foundin.hdr", checkIfExists: true).collect()
    ch_cadd_resources           = params.cadd_resources                    ? Channel.fromPath(params.cadd_resources).collect()
                                                                           : Channel.value([])
    ch_bait_intervals           = ch_references.bait_intervals
    ch_genome_bwaindex          = params.bwa                                ? Channel.fromPath(params.bwa).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                            ? Channel.fromPath(params.bwamem2).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_bwamem2_index
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = params.sequence_dictionary                ? Channel.fromPath(params.sequence_dictionary).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_dict
    ch_target_bed               = ch_references.target_bed
    known_indels_tbi            = params.known_indels   ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_REFERENCES.out.known_indels_tbi      : Channel.value([])
    known_indels_2_tbi          = params.known_indels_2 ? params.known_indels_2_tbi    ? Channel.fromPath(params.known_indels_2_tbi).collect()    : PREPARE_REFERENCES.out.known_indels_2_tbi    : Channel.value([])
    dbsnp_tbi                   = params.dbsnp          ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_REFERENCES.out.dbsnp_tbi             : Channel.value([])
    dict                        = params.dict           ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                                        : PREPARE_REFERENCES.out.genome_dict
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_target_intervals         = ch_references.target_intervals
    ch_gnomad_afidx             = params.gnomad_af_idx                     ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                           : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                         ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                           : Channel.empty()
    ch_reduced_penetrance       = params.reduced_penetrance                ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                           : Channel.value([])
    ch_score_config_snv         = params.score_config_snv                  ? Channel.fromPath(params.score_config_snv).collect()
                                                                           : Channel.value([])
    ch_score_config_mt          = params.score_config_mt                    ? Channel.fromPath(params.score_config_mt).collect()
                                                                            : Channel.value([])
    ch_vcfanno_resources        = params.vcfanno_resources                 ? Channel.fromPath(params.vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                                           : Channel.value([])
    ch_vcfanno_lua              = params.vcfanno_lua                       ? Channel.fromPath(params.vcfanno_lua).collect()
                                                                           : Channel.value([])
    ch_vcfanno_toml             = params.vcfanno_toml                      ? Channel.fromPath(params.vcfanno_toml).collect()
                                                                           : Channel.value([])
    ch_vep_cache                = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                           : ( params.vep_cache    ? Channel.fromPath(params.vep_cache).collect() : Channel.value([]) )
    ch_vep_filters              = params.vep_filters                       ? Channel.fromPath(params.vep_filters).collect()
                                                                           : Channel.value([])
    ch_vep_extra_files_unsplit  = params.vep_plugin_files                  ? Channel.fromPath(params.vep_plugin_files).collect()
                                                                           : Channel.value([])
    ch_vep_filters_std_fmt      = params.vep_filters                        ? Channel.fromPath(params.vep_filters).map { it -> [[id:'standard'],it]}.collect()
                                                                            : Channel.empty()
    ch_vep_filters_scout_fmt    = params.vep_filters_scout_fmt              ? Channel.fromPath(params.vep_filters_scout_fmt).map { it -> [[id:'scout'],it]}.collect()
                                                                            : Channel.empty()
    ch_variant_consequences_snv = params.variant_consequences_snv           ? Channel.fromPath(params.variant_consequences_snv).collect()
                                                                            : Channel.value([])
    ch_mtshift_bwamem2index     = ch_references.mtshift_bwamem2_index
    ch_mtshift_dictionary       = ch_references.mtshift_dict
    ch_mtshift_fai              = ch_references.mtshift_fai
    ch_mtshift_fasta            = ch_references.mtshift_fasta
    ch_mtshift_intervals        = ch_references.mtshift_intervals
    ch_mt_intervals             = ch_references.mt_intervals
    ch_mtshift_backchain        = ch_references.mtshift_backchain                         
    ch_versions                 = ch_versions.mix(ch_references.versions)

    // Read and store paths in the vep_plugin_files file
    if (params.vep_plugin_files) {
        ch_vep_extra_files_unsplit.splitCsv ( header:true )
            .map { row ->
                f = file(row.vep_files[0])
                if(f.isFile() || f.isDirectory()){
                    return [f]
                } else {
                    error("\nVep database file ${f} does not exist.")
                }
            }
            .collect()
            .set {ch_vep_extra_files}
    }

    // Read and store hgnc ids in a channel
    ch_vep_filters_scout_fmt
        .mix (ch_vep_filters_std_fmt)
        .set {ch_vep_filters}

    CREATE_HGNCIDS_FILE(ch_vep_filters)
        .txt
        .set {ch_hgnc_ids}

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    known_sites_indels     = known_indels_2.concat(known_indels).collect()
    known_sites_indels_tbi = known_indels_2_tbi.concat(known_indels_tbi).collect()

    // Generate pedigree file
    ch_pedfile   = CREATE_PEDIGREE_FILE(ch_samples.toList()).ped
    ch_versions = ch_versions.mix(CREATE_PEDIGREE_FILE.out.versions)

    
    // CREATE CHROMOSOME BED AND INTERVALS
    SCATTER_GENOME (
        ch_genome_dictionary,
        ch_genome_fai,
        ch_genome_fasta
    )
    .set { ch_scatter }

    ch_scatter_split_intervals  = ch_scatter.split_intervals  ?: Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_reports = ch_reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    
    //
    // ALIGNING READS, FETCH STATS, AND MERGE.
    //
    ALIGN (
        ch_samplesheet,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        ch_genome_dictionary,
        ch_mtshift_bwamem2index,
        ch_mtshift_fasta,
        ch_mtshift_dictionary,
        ch_mtshift_fai,
        params.platform
    )
    .set { ch_aligned }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)
    mapped_bam_bai = ALIGN.out.genome_bam_bai

    //BASERECALIBRATOR
    BAM_BASERECALIBRATOR(
                mapped_bam_bai,
                dict,
                ch_genome_fasta,
                ch_genome_fai,
                known_sites_indels,
                known_sites_indels_tbi
    )
    ch_versions = ch_versions.mix(BAM_BASERECALIBRATOR.out.versions)
    ch_table_bqsr = BAM_BASERECALIBRATOR.out.table_bqsr
    bam_applybqsr = mapped_bam_bai.join(ch_table_bqsr)

    //APPLYBQSR
    BAM_APPLYBQSR(
                bam_applybqsr,
                dict,
                ch_genome_fasta,
                ch_genome_fai,
    )
    .set{ch_recalibrated}

    ch_versions = ch_versions.mix(BAM_APPLYBQSR.out.versions)

    // BAM QUALITY CHECK
    QC_BAM (
        ch_recalibrated.bqsr_bam,
        ch_recalibrated.bqsr_bai,
        ch_recalibrated.a_bqsr_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_bait_intervals,
        ch_target_intervals,
        ch_genome_chrsizes,
        ch_intervals_wgs,
        ch_intervals_y,
        dict
    )

    ch_versions = ch_versions.mix(QC_BAM.out.versions)

    input_bam = ch_recalibrated.a_bqsr_bam_bai

    // Assign BAM
    input_bam.combine(chh_target_bed)
                .map{ meta, bam, bai, interval -> [ meta, bam, bai, interval]
            }.set{ch_haplotypecaller_interval_bam}

    ch_haplotypecaller_interval_bam
        .collect { it[1] } // BAM files are at index 1 in the list
        .toList()
        .set { ch_bam_list }

    ch_haplotypecaller_interval_bam
        .collect { it[2] } // BAI files are at index 1 in the list
        .toList()
        .set { ch_bai_list }

    ch_case_info
        .combine(ch_bam_list)
        .combine(ch_bai_list)
        .set { ch_bams_bais }

    SMNCOPYNUMBERCALLER (
        ch_bams_bais
    )
    ch_versions = ch_versions.mix(SMNCOPYNUMBERCALLER.out.versions)        

    //
    // SNV Variant calling
    //
    CALL_SNV (
        ch_haplotypecaller_interval_bam,
        ch_aligned.mt_bam_bai,
        ch_aligned.mtshift_bam_bai,
        ch_genome_chrsizes,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary,
        ch_mt_intervals,
        ch_mtshift_fasta,
        ch_mtshift_fai,
        ch_mtshift_dictionary,
        ch_mtshift_intervals,
        ch_mtshift_backchain,
        dbsnp,
        dbsnp_tbi,
        ch_case_info,
        ch_foundin_header,
    )
    ch_versions = ch_versions.mix(CALL_SNV.out.versions)

    //
    // ANNOTATE GENOME SNVs
    //
    if (!params.skip_snv_annotation) {

        ANNOTATE_GENOME_SNVS (
            CALL_SNV.out.genome_vcf_tabix,
            params.analysis_type,
            ch_cadd_header,
            ch_cadd_resources,
            ch_vcfanno_resources,
            ch_vcfanno_lua,
            ch_vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_gnomad_af,
            ch_samples,
            ch_scatter_split_intervals,
            ch_vep_extra_files,
            ch_genome_chrsizes
        ).set { ch_snv_annotate }
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)

        GENERATE_CLINICAL_SET_SNV(
            ch_snv_annotate.vcf_ann,
            ch_hgnc_ids
        )
        ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_SNV.out.versions)

        ANN_CSQ_PLI_SNV (
            GENERATE_CLINICAL_SET_SNV.out.vcf,
            ch_variant_consequences_snv
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

        RANK_VARIANTS_SNV (
            ANN_CSQ_PLI_SNV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_snv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

    }

    //
    // ANNOTATE MT SNVs
    //
    if (!params.skip_mt_annotation && (params.run_mt_for_wes || params.analysis_type.equals("wgs"))) {

        ANNOTATE_MT_SNVS (
            CALL_SNV.out.mt_vcf,
            CALL_SNV.out.mt_tabix,
            ch_cadd_header,
            ch_cadd_resources,
            ch_genome_fasta,
            ch_vcfanno_resources,
            ch_vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_vep_extra_files
        ).set { ch_mt_annotate }
        ch_versions = ch_versions.mix(ch_mt_annotate.versions)

        GENERATE_CLINICAL_SET_MT(
            ch_mt_annotate.vcf_ann,
            ch_hgnc_ids
        )
        ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_MT.out.versions)

        ANN_CSQ_PLI_MT(
            GENERATE_CLINICAL_SET_MT.out.vcf,
            ch_variant_consequences_snv
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_MT.out.versions)

        RANK_VARIANTS_MT (
            ANN_CSQ_PLI_MT.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_mt
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_MT.out.versions)

    }

    MERGE_VCFS(
        RANK_VARIANTS_SNV.out.vcf,
        RANK_VARIANTS_MT.out.vcf 
    ) 


    //
    // MODULE: Pipeline reporting
    //

    // The template v2.7.1 template update introduced: ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    // This caused the pipeline to stall
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRaredisease.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRaredisease.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.ped.map{it[1]}.collect().ifEmpty([]))
//    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.csv.map{it[1]}.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
