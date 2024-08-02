//
// Align MT
//


include { BWAMEM2_MEM as BWAMEM2_MEM_MT                                     } from '../../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MERGEBAMALIGNMENT as GATK4_MERGEBAMALIGNMENT_MT             } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_MT } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT                 } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT                               } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MT                                 } from '../../../modules/nf-core/samtools/sort/main'

workflow ALIGN_MT {
    take:
        ch_fastq        // channel: [mandatory] [ val(meta), [ path(reads) ] ]
        ch_ubam         // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bwamem2index // channel: [mandatory for bwamem2] [ val(meta), path(index) ]
        ch_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_dict         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_fai          // channel: [mandatory] [ val(meta), path(fai) ]

    main:
        ch_versions     = Channel.empty()

        if (params.aligner.equals("bwamem2")) {
            BWAMEM2_MEM_MT (ch_fastq, ch_bwamem2index, ch_fasta, true)
            ch_align       = BWAMEM2_MEM_MT.out.bam
            ch_versions    = ch_versions.mix(BWAMEM2_MEM_MT.out.versions.first())
        }
        ch_align
            .join(ch_ubam, failOnMismatch:true, failOnDuplicate:true)
            .set {ch_bam_ubam}

        GATK4_MERGEBAMALIGNMENT_MT (ch_bam_ubam, ch_fasta, ch_dict)

        PICARD_ADDORREPLACEREADGROUPS_MT (GATK4_MERGEBAMALIGNMENT_MT.out.bam, [[:],[]], [[:],[]])

        PICARD_MARKDUPLICATES_MT (PICARD_ADDORREPLACEREADGROUPS_MT.out.bam, ch_fasta, ch_fai)

        SAMTOOLS_SORT_MT (PICARD_MARKDUPLICATES_MT.out.bam, [[:],[]])

        SAMTOOLS_INDEX_MT(SAMTOOLS_SORT_MT.out.bam)

        ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MT.out.versions.first())

    emit:
        marked_bam  = SAMTOOLS_SORT_MT.out.bam   // channel: [ val(meta), path(bam) ]
        marked_bai  = SAMTOOLS_INDEX_MT.out.bai  // channel: [ val(meta), path(bai) ]
        versions    = ch_versions                // channel: [ path(versions.yml) ]
}
