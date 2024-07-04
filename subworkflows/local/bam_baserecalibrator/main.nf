//
// PREPARE RECALIBRATION
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
//include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'

workflow BAM_BASERECALIBRATOR {
    take:
    bam            // channel: [mandatory] [ meta, cram_markduplicates, crai ]
    dict            // channel: [mandatory] [ dict ]
    fasta           // channel: [mandatory] [ fasta ]
    fasta_fai       // channel: [mandatory] [ fasta_fai ]
    known_sites     // channel: [optional]  [ known_sites ]
    known_sites_tbi // channel: [optional]  [ known_sites_tbi ]

    main:
    versions = Channel.empty()


    // RUN BASERECALIBRATOR
    GATK4_BASERECALIBRATOR(bam, fasta, fasta_fai, dict.map{ meta, it -> [ it ] }, known_sites, known_sites_tbi)

    // Gather versions of all tools used
    versions = versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    emit:
    table_bqsr = GATK4_BASERECALIBRATOR.out.table  // channel: [ meta, table ]

    versions   // channel: [ versions.yml ]
}
