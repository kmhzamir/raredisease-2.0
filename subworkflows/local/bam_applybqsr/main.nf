//
// RECALIBRATE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_APPLYBQSR           } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_APPLYBQSR {
    take:
    bam          // channel: [mandatory] [ meta, cram, crai, recal ]
    dict          // channel: [mandatory] [ dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]

    main:
    versions = Channel.empty()
    ch_bqsr_bam     = Channel.empty()
    ch_bqsr_bai     = Channel.empty()
    ch_bqsr_bam_bai = Channel.empty()


    // RUN APPLYBQSR
    GATK4_APPLYBQSR(bam, fasta, fasta_fai, dict.map{ meta, it -> [ it ] })
    SAMTOOLS_INDEX_MARKDUP ( GATK4_APPLYBQSR.out.bam )

    ch_bqsr_bam = GATK4_APPLYBQSR.out.bam
    ch_bqsr_bai = SAMTOOLS_INDEX_MARKDUP.out.bai
    ch_bqsr_bam_bai = ch_bqsr_bam.join(ch_bqsr_bai, failOnMismatch:true, failOnDuplicate:true)
    

    // Gather versions of all tools used
    versions = versions.mix(GATK4_APPLYBQSR.out.versions)

    emit:
    bqsr_bam  = ch_bqsr_bam         // channel: [ val(meta), path(bam) ]
    bqsr_bai  = ch_bqsr_bai         // channel: [ val(meta), path(bai) ]
    a_bqsr_bam_bai = ch_bqsr_bam_bai

    versions          // channel: [ versions.yml ]
}
