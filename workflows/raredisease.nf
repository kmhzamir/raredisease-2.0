/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

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
include { PREPARE_REFERENCES                    } from '../subworkflows/local/prepare_references'
include { SCATTER_GENOME                        } from '../subworkflows/local/scatter_genome'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { BAM_MARKDUPLICATES                    } from '../subworkflows/local/bam_markduplicates/main'
include { BAM_BASERECALIBRATOR                  } from '../subworkflows/local/bam_baserecalibrator/main'
include { BAM_APPLYBQSR                         } from '../subworkflows/local/bam_applybqsr/main'
include { QC_BAM                                } from '../subworkflows/local/qc_bam'
include { BAM_VARIANT_CALLING_HAPLOTYPECALLER   } from '../subworkflows/local/bam_variant_calling_haplotypecaller'
include { GATK4_HAPLOTYPECALLER                 } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { SMNCOPYNUMBERCALLER                   } from '../modules/nf-core/smncopynumbercaller/main'
include { FILTER_VEP as FILTER_VEP_SNV          } from '../modules/local/filter_vep'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV   } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_SNVS                         } from '../subworkflows/local/annotate_snvs'
include { RANK_VARIANTS as RANK_VARIANTS_SNV    } from '../subworkflows/local/rank_variants'
include { CREATE_PEDIGREE_FILE                  } from '../modules/local/create_pedigree_file'

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

    // Prepare references and indices.
    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_genome_fai,
        ch_target_bed_unprocessed,
        dbsnp,
        known_indels,
        ch_gnomad_af_tab,
        ch_vep_cache_unprocessed
    )
    .set { ch_references }

    // Gather built indices or get them from the params
    ch_cadd_header              = Channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_cadd_resources           = params.cadd_resources                    ? Channel.fromPath(params.cadd_resources).collect()
                                                                           : Channel.value([])
    ch_bait_intervals           = ch_references.bait_intervals
    ch_genome_bwaindex          = params.bwa                                ? Channel.fromPath(params.bwa).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                            ? Channel.fromPath(params.bwamem2).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_bwamem2_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = params.sequence_dictionary                ? Channel.fromPath(params.sequence_dictionary).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_dict
    ch_target_bed               = ch_references.target_bed
    known_indels_tbi            = params.known_indels   ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_REFERENCES.out.known_indels_tbi      : Channel.value([])
    dbsnp_tbi                   = params.dbsnp          ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_REFERENCES.out.dbsnp_tbi             : Channel.value([])
    dict                        = params.dict           ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                                        : PREPARE_REFERENCES.out.genome_dict
    ch_genome_chrsizes          = PREPARE_REFERENCES.out.genome_chrom_sizes
    ch_target_intervals         = ch_references.target_intervals
    ch_gnomad_afidx             = params.gnomad_af_idx                     ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                           : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                         ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                           : Channel.empty()
    ch_reduced_penetrance       = params.reduced_penetrance                ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                           : Channel.value([])
    ch_score_config_snv         = params.score_config_snv                  ? Channel.fromPath(params.score_config_snv).collect()
                                                                           : Channel.value([])
    ch_variant_consequences     = Channel.fromPath("$projectDir/assets/variant_consequences_v1.txt", checkIfExists: true).collect()
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
    chh_pedfile                 = params.ped_file                          ? Channel.fromPath(params.ped_file).collect()
                                                                           : Channel.value([])                                                           
    ch_versions                 = ch_versions.mix(ch_references.versions)

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

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
    // MODULE: Run Fastp trimming
    //
    // Trimming and/or splitting
    if (params.trim_fastq || params.split_fastq > 0) {
        save_trimmed_fail = false
        save_merged = false
        FASTP(
            ch_samplesheet,
            [],
            save_trimmed_fail,
            save_merged
        )

        ch_reports = ch_reports.mix(FASTP.out.json.collect{ meta, json -> json })
        ch_reports = ch_reports.mix(FASTP.out.html.collect{ meta, html -> html })

        if (params.split_fastq) {
            reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                [ meta + [ size:read_files.size() ], read_files ]
            }.transpose()
        } else reads_for_alignment = FASTP.out.reads

        ch_versions = ch_versions.mix(FASTP.out.versions)

    } 
    else {
       reads_for_alignment = ch_samplesheet
    }
    
    //
    // ALIGNING READS, FETCH STATS, AND MERGE.
    //
    ALIGN (
        reads_for_alignment,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        params.platform
    )
    .set { ch_aligned }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)


    //MARKDUPLICATES AND REMOVE DUPLICATES
    BAM_MARKDUPLICATES(
                ALIGN.out.ch_bam_bai,
                ch_genome_fasta,
                ch_genome_fai,
    )
    .set {ch_mapped}

    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES.out.versions)
    mapped_bam_bai = BAM_MARKDUPLICATES.out.a_genome_bam_bai

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
    if (params.analysis_type == "wes") {
        BAM_VARIANT_CALLING_HAPLOTYPECALLER(
            ch_haplotypecaller_interval_bam, 
            ch_genome_fasta, 
            ch_genome_fai, 
            dict, 
            dbsnp, 
            dbsnp_tbi
        )
    }
    ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.ch_versions)

    // VARIANT ANNOTATION

    if (!params.skip_snv_annotation) {
        ANNOTATE_SNVS (
            BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.vcf_tbi,
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
            ch_scatter_split_intervals
        ).set {ch_snv_annotate}
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)

        ch_snv_annotate = ANNOTATE_SNVS.out.vcf_ann

        ANN_CSQ_PLI_SNV (
            ch_snv_annotate,
            ch_variant_consequences
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

        RANK_VARIANTS_SNV (
            ANN_CSQ_PLI_SNV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_snv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

        FILTER_VEP_SNV(
            RANK_VARIANTS_SNV.out.vcf,
            ch_vep_filters
        )
        ch_versions = ch_versions.mix(FILTER_VEP_SNV.out.versions)

    } 

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
