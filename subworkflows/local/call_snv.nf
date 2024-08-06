//
// call Single-nucleotide Varinats
//

include { CALL_SNV_MT                      } from './variant_calling/call_snv_MT'
include { CALL_SNV_MT as CALL_SNV_MT_SHIFT } from './variant_calling/call_snv_MT'
include { POSTPROCESS_MT_CALLS             } from './variant_calling/postprocess_MT_calls'
include { GATK4_SELECTVARIANTS             } from '../../modules/nf-core/gatk4/selectvariants/main'
include { BAM_VARIANT_CALLING_HAPLOTYPECALLER   } from './variant_calling/bam_variant_calling_haplotypecaller_2'


workflow CALL_SNV {
    take:
        ch_genome_bam_bai     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mt_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mtshift_bam_bai    // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_chrsizes    // channel: [mandatory] [ path(sizes) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_dictionary  // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mt_intervals       // channel: [optional] [ path(interval_list) ]
        ch_mtshift_fasta      // channel: [optional] [ val(meta), path(fasta) ]
        ch_mtshift_fai        // channel: [optional] [ val(meta), path(fai) ]
        ch_mtshift_dictionary // channel: [optional] [ val(meta), path(dict) ]
        ch_mtshift_intervals  // channel: [optional] [ path(interval_list) ]
        ch_mtshift_backchain  // channel: [mandatory] [ val(meta), path(back_chain) ]
        ch_dbsnp              // channel: [optional] [ val(meta), path(vcf) ]
        ch_dbsnp_tbi          // channel: [optional] [ val(meta), path(tbi) ]
        ch_case_info          // channel: [mandatory] [ val(case_info) ]
        ch_foundin_header     // channel: [mandatory] [ path(header) ]

    main:
        ch_versions      = Channel.empty()
        ch_mt_vcf        = Channel.empty()
        ch_mt_tabix      = Channel.empty()
        ch_genome_vcf_tabix    = Channel.empty()


        BAM_VARIANT_CALLING_HAPLOTYPECALLER(
            ch_genome_bam_bai, 
            ch_genome_fasta, 
            ch_genome_fai, 
            ch_genome_dictionary, 
            ch_dbsnp, 
            ch_dbsnp_tbi,
            ch_case_info
        )
        ch_vcf              = BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.vcf
        ch_tabix            = BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.tbi

        ch_vcf
            .join(ch_tabix, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set {ch_selvar_in}
        GATK4_SELECTVARIANTS(ch_selvar_in) // remove mitochondrial variants

        ch_genome_vcf       = GATK4_SELECTVARIANTS.out.vcf
        ch_genome_tabix     = GATK4_SELECTVARIANTS.out.tbi
        ch_genome_vcf_tabix = ch_genome_vcf.join(ch_genome_tabix, failOnMismatch:true, failOnDuplicate:true)

        if (params.analysis_type.equals("wgs") || params.run_mt_for_wes) {
            CALL_SNV_MT(
                ch_mt_bam_bai,
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_dictionary,
                ch_mt_intervals
            )

            CALL_SNV_MT_SHIFT(
                ch_mtshift_bam_bai,
                ch_mtshift_fasta,
                ch_mtshift_fai,
                ch_mtshift_dictionary,
                ch_mtshift_intervals
            )

            POSTPROCESS_MT_CALLS(
                CALL_SNV_MT.out.vcf,
                CALL_SNV_MT_SHIFT.out.vcf,
                ch_genome_fasta,
                ch_genome_dictionary,
                ch_genome_fai,
                ch_mtshift_backchain,
                ch_case_info,
                ch_foundin_header,
                ch_genome_chrsizes
            )
            ch_mt_vcf   = POSTPROCESS_MT_CALLS.out.vcf
            ch_mt_tabix = POSTPROCESS_MT_CALLS.out.tbi
            ch_versions = ch_versions.mix(CALL_SNV_MT.out.versions)
            ch_versions = ch_versions.mix(CALL_SNV_MT_SHIFT.out.versions)
            ch_versions = ch_versions.mix(POSTPROCESS_MT_CALLS.out.versions)
            //ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)
            ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.ch_versions)
        }

    emit:
        genome_vcf       = ch_genome_vcf       // channel: [ val(meta), path(vcf) ]
        genome_tabix     = ch_genome_tabix     // channel: [ val(meta), path(tbi) ]
        genome_vcf_tabix = ch_genome_vcf_tabix // channel: [ val(meta), path(vcf), path(tbi) ]
        mt_vcf           = ch_mt_vcf           // channel: [ val(meta), path(vcf) ]
        mt_tabix         = ch_mt_tabix         // channel: [ val(meta), path(tbi) ]
        versions         = ch_versions         // channel: [ path(versions.yml) ]
}
