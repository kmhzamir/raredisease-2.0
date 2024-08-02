//
// Map to reference
//

include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'
include { ALIGN_MT                   } from './alignment/align_MT'
include { ALIGN_MT as ALIGN_MT_SHIFT } from './alignment/align_MT'
include { CONVERT_MT_BAM_TO_FASTQ    } from './mitochondria/convert_mt_bam_to_fastq'

workflow PROCESS_MT {
    take:
        ch_reads_input     // channel: [mandatory] [ val(meta), [path(reads)]  ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_bwamem2_index   // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_dictionary     // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_bwamem2index  // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mtshift_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_fai           // channel: [mandatory] [ val(meta), path(fai) ]

    main:
        ch_mt_bam_bai         = Channel.empty()
        ch_mt_marked_bam      = Channel.empty()
        ch_mt_marked_bai      = Channel.empty()
        ch_mtshift_bam_bai    = Channel.empty()
        ch_mtshift_marked_bam = Channel.empty()
        ch_mtshift_marked_bai = Channel.empty()
        ch_versions           = Channel.empty()

        // PREPARING READS FOR MT ALIGNMENT

        if (params.analysis_type.equals("wgs") || params.run_mt_for_wes) {
            CONVERT_MT_BAM_TO_FASTQ (
                ch_reads_input,
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_dictionary
            )

            ALIGN_MT (
                CONVERT_MT_BAM_TO_FASTQ.out.fastq,
                CONVERT_MT_BAM_TO_FASTQ.out.bam,
                ch_bwamem2_index,
                ch_genome_fasta,
                ch_genome_dictionary,
                ch_genome_fai
            )

            ALIGN_MT_SHIFT (
                CONVERT_MT_BAM_TO_FASTQ.out.fastq,
                CONVERT_MT_BAM_TO_FASTQ.out.bam,
                ch_mtshift_bwamem2index,
                ch_mtshift_fasta,
                ch_mtshift_dictionary,
                ch_mtshift_fai
            )

            ch_mt_marked_bam      = ALIGN_MT.out.marked_bam
            ch_mt_marked_bai      = ALIGN_MT.out.marked_bai
            ch_mt_bam_bai         = ch_mt_marked_bam.join(ch_mt_marked_bai, failOnMismatch:true, failOnDuplicate:true)
            ch_mtshift_marked_bam = ALIGN_MT_SHIFT.out.marked_bam
            ch_mtshift_marked_bai = ALIGN_MT_SHIFT.out.marked_bai
            ch_mtshift_bam_bai    = ch_mtshift_marked_bam.join(ch_mtshift_marked_bai, failOnMismatch:true, failOnDuplicate:true)
            ch_versions           = ch_versions.mix(ALIGN_MT.out.versions,
                                        ALIGN_MT_SHIFT.out.versions,
                                        CONVERT_MT_BAM_TO_FASTQ.out.versions)
        }


    emit:
        mt_marked_bam      = ch_mt_marked_bam      // channel: [ val(meta), path(bam) ]
        mt_marked_bai      = ch_mt_marked_bai      // channel: [ val(meta), path(bai) ]
        mt_bam_bai         = ch_mt_bam_bai         // channel: [ val(meta), path(bam), path(bai) ]
        mtshift_marked_bam = ch_mtshift_marked_bam // channel: [ val(meta), path(bam) ]
        mtshift_marked_bai = ch_mtshift_marked_bai // channel: [ val(meta), path(bai) ]
        mtshift_bam_bai    = ch_mtshift_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
        versions           = ch_versions           // channel: [ path(versions.yml) ]
}
