//
// GATK4 HAPLOTYPACALLER germline variant calling:
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
//

// include { BAM_MERGE_INDEX_SAMTOOLS    } from '../bam_merge_index_samtools/main'
include { GATK4_HAPLOTYPECALLER       } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS             } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_HAPLOTYPECALLER {
    take:
    ch_bam_bai          // channel: [mandatory] [ val(meta), path(bam), path(bai), path(interval), path(dragstr_model) ]
    ch_fasta            // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai              // channel: [mandatory] [ val(meta), path(fai) ]
    ch_dict             // channel: [mandatory] [ val(meta), path(sequence dictionary) ]
    ch_dbsnp            // channel: [mandatory] [ val(meta), path(dbsnp) ]
    ch_dbsnp_tbi        // channel: [mandatory] [ val(meta), path(dbsnp_tbi) ]
    ch_case_info        // channel: [mandatory] [ val(case_info) ]

    main:
    ch_versions = Channel.empty()
    vcf_tbi     = Channel.empty()

    // Run haplotypecaller
    GATK4_HAPLOTYPECALLER(ch_bam_bai, ch_fasta, ch_fai, ch_dict)

    // For joint genotyping
    gvcf_tbi_intervals = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true)
        .join(ch_bam_bai, failOnMismatch: true)
        .map{ meta, gvcf, tbi, cram, crai, intervals -> [ meta, gvcf, tbi, intervals ] }

    // Figuring out if there is one or more vcf(s) from the same sample
    haplotypecaller_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{
            meta, vcf -> [ meta - meta.subMap('interval_name'), vcf]
        }
        .branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Figuring out if there is one or more tbi(s) from the same sample
    haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi.map{
            meta, tbi -> [ meta - meta.subMap('interval_name'), tbi]
        }.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Combine ch_case_info with the VCF output before merging
    combined_vcf_intervals = haplotypecaller_vcf.intervals
        .combine(ch_case_info)
        .map{ vcf, case_info -> [ case_info, vcf[1] ] }
    
    // Only when using intervals
    GATK4_MERGEVCFS(combined_vcf_intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple(), ch_dict)

    haplotypecaller_vcf = Channel.empty().mix(
            GATK4_MERGEVCFS.out.vcf,
            haplotypecaller_vcf.no_intervals)

    haplotypecaller_tbi = Channel.empty().mix(
            GATK4_MERGEVCFS.out.tbi,
            haplotypecaller_tbi.no_intervals)

    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

    // Remove no longer necessary field: num_intervals
    vcf = haplotypecaller_vcf.map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }
    tbi = haplotypecaller_tbi.map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }

    emit:
    gvcf_tbi_intervals // For joint genotyping
    vcf   // channel: [ val(meta), path(vcf) ]
    tbi   // channel: [ val(meta), path(tbi) ]
    vcf_tbi = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true)
    ch_versions                   // channel: [ path(versions.yml) ]
}
