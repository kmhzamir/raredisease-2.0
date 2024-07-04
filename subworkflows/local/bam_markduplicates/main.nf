//
// MARKDUPLICATES AND QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { PICARD_MARKDUPLICATES as MARKDUPLICATES  } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES {
    take:
    bam                    // channel: [mandatory] [ meta, bam, bai ]
    fasta                  // channel: [mandatory] [ fasta ]
    fasta_fai              // channel: [mandatory] [ fasta_fai ]

    main:
    ch_versions       = Channel.empty()
    ch_marked_bam     = Channel.empty()
    ch_marked_bai     = Channel.empty()
    ch_genome_bam_bai = Channel.empty()


    // Marking duplicates
    MARKDUPLICATES ( bam , fasta, fasta_fai )
    SAMTOOLS_INDEX_MARKDUP ( MARKDUPLICATES.out.bam )
    ch_marked_bam = MARKDUPLICATES.out.bam
    ch_marked_bai = SAMTOOLS_INDEX_MARKDUP.out.bai
    ch_genome_bam_bai = ch_marked_bam.join(ch_marked_bai, failOnMismatch:true, failOnDuplicate:true)


    ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MARKDUP.out.versions.first())

    emit:
    metrics     = MARKDUPLICATES.out.metrics     // channel: [ val(meta), path(metrics) ]
    marked_bam  = ch_marked_bam         // channel: [ val(meta), path(bam) ]
    marked_bai  = ch_marked_bai // channel: [ val(meta), path(bai) ]
    a_genome_bam_bai = ch_genome_bam_bai
    versions    = ch_versions
}
