//
// Prepare reference files
//

include { BWA_INDEX as BWA_INDEX_GENOME                      } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME              } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_MT_SHIFT            } from '../../modules/nf-core/bwamem2/index/main'
include { CAT_CAT as CAT_CAT_BAIT                            } from '../../modules/nf-core/cat/cat/main'
include { GATK4_BEDTOINTERVALLIST as GATK_BILT               } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD          } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD_MT_SHIFT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS as GATK_ILT                } from '../../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_PREPROCESSINTERVALS as GATK_PREPROCESS_WGS   } from '../../modules/nf-core/gatk4/preprocessintervals/main.nf'
include { GATK4_PREPROCESSINTERVALS as GATK_PREPROCESS_WES   } from '../../modules/nf-core/gatk4/preprocessintervals/main.nf'
include { GATK4_SHIFTFASTA as GATK_SHIFTFASTA                } from '../../modules/nf-core/gatk4/shiftfasta/main'
include { GET_CHROM_SIZES                                    } from '../../modules/local/get_chrom_sizes'
include { SAMTOOLS_FAIDX as SAMTOOLS_EXTRACT_MT              } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME            } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MT_SHIFT          } from '../../modules/nf-core/samtools/faidx/main'
//include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_GENOME      } from '../../modules/nf-core/sentieon/bwaindex/main'
//include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_MT_SHIFT    } from '../../modules/nf-core/sentieon/bwaindex/main'
include { TABIX_BGZIPTABIX as TABIX_PBT                      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_GNOMAD_AF                     } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PT                            } from '../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS                  } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS_2                } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DBSNP                         } from '../../modules/nf-core/tabix/tabix/main'

workflow PREPARE_REFERENCES {
    take:
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mt_fasta        // channel: [mandatory for dedicated mt analysis] [ val(meta), path(fasta) ]
        ch_target_bed      // channel: [mandatory for WES] [ path(bed) ]
        dbsnp             // channel: [optional]  dbsnp
        known_indels      // channel: [optional]  known_indels
        known_indels_2      // channel: [optional]  known_indels
        ch_gnomad_af_tab   // channel: [optional; used in for snv annotation] [ val(meta), path(tab) ]
        ch_vep_cache       // channel: [mandatory for annotation] [ path(cache) ]



    main:
        ch_versions    = Channel.empty()
        ch_tbi         = Channel.empty()
        ch_bgzip_tbi   = Channel.empty()
        ch_bwa         = Channel.empty()
        //ch_sentieonbwa = Channel.empty()

        // Genome indices
        SAMTOOLS_FAIDX_GENOME(ch_genome_fasta, [[],[]])
        GATK_SD(ch_genome_fasta)
        ch_fai = Channel.empty().mix(ch_genome_fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        GET_CHROM_SIZES( ch_fai )
        ch_genome_fasta.map { meta, fasta -> return [meta, fasta, [], [] ] }
            .set {ch_rtgformat_in}

       // MT genome indices
        SAMTOOLS_EXTRACT_MT(ch_genome_fasta, ch_fai)
        ch_mt_fasta_in = Channel.empty().mix(ch_mt_fasta, SAMTOOLS_EXTRACT_MT.out.fa).collect()
        SAMTOOLS_FAIDX_MT_SHIFT(ch_mt_fasta_in, [[],[]])
        GATK_SD_MT_SHIFT(ch_mt_fasta_in)
        GATK_SHIFTFASTA(ch_mt_fasta_in, SAMTOOLS_FAIDX_MT_SHIFT.out.fai, GATK_SD_MT_SHIFT.out.dict)

        // MT alignment indices
        BWAMEM2_INDEX_MT_SHIFT(GATK_SHIFTFASTA.out.shift_fa)
        GATK_SHIFTFASTA.out.intervals
            .multiMap{ meta, files ->
                    shift_intervals:
                        ind = files.findIndexValues {it.toString().endsWith("shifted.intervals")}
                        files[ind]
                    intervals:
                        ind = files.findIndexValues {!(it.toString().endsWith("shifted.intervals"))}
                        files[ind]
            }
            .set {ch_shiftfasta_mtintervals}

        //bed indices
        //TABIX_PT(ch_target_bed).tbi.set { ch_tbi }
        TABIX_PBT(ch_target_bed).gz_tbi.set { ch_bgzip_tbi }
        TABIX_GNOMAD_AF(ch_gnomad_af_tab)

        // Generate bait and target intervals
        GATK_BILT(ch_target_bed, GATK_SD.out.dict).interval_list
        GATK_ILT(GATK_BILT.out.interval_list)
        GATK_ILT.out.interval_list
            .collect{ it[1] }
            .map { it ->
                meta = it[0].toString().split("_split")[0].split("/")[-1] + "_bait.intervals_list"
                return [[id:meta], it]
            }
            .set { ch_bait_intervals_cat_in } 
        CAT_CAT_BAIT ( ch_bait_intervals_cat_in )
        TABIX_DBSNP(dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
        TABIX_KNOWN_INDELS(known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
        TABIX_KNOWN_INDELS_2(known_indels_2.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
        UNTAR_VEP_CACHE (ch_vep_cache)

        // Gather versions
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_versions = ch_versions.mix(GATK_SD.out.versions)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(GATK_BILT.out.versions)
        ch_versions = ch_versions.mix(CAT_CAT_BAIT.out.versions)
        ch_versions = ch_versions.mix(GATK_ILT.out.versions)
        ch_versions = ch_versions.mix(TABIX_GNOMAD_AF.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_EXTRACT_MT.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_MT_SHIFT.out.versions)
        ch_versions = ch_versions.mix(GATK_SD_MT_SHIFT.out.versions)
        ch_versions = ch_versions.mix(GATK_SHIFTFASTA.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_MT_SHIFT.out.versions)

    emit:
        genome_bwa_index      = Channel.empty().mix(ch_bwa).collect()            // channel: [ val(meta), path(index) ]
        //genome_bwamem2_index  = BWAMEM2_INDEX_GENOME.out.index.collect()                         // channel: [ val(meta), path(index) ]
        genome_chrom_sizes    = GET_CHROM_SIZES.out.sizes.collect()                              // channel: [ path(sizes) ]
        genome_fai            = ch_fai                                                           // channel: [ val(meta), path(fai) ]
        genome_dict           = GATK_SD.out.dict.collect()                                       // channel: [ path(dict) ]
        target_bed            = Channel.empty().mix(ch_bgzip_tbi).collect()              // channel: [ val(meta), path(bed), path(tbi) ]
        target_intervals      = GATK_BILT.out.interval_list.map{ meta, inter -> inter}.collect() // channel: [ path(interval_list) ]
        known_indels_tbi      = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        known_indels_2_tbi    = TABIX_KNOWN_INDELS_2.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        dbsnp_tbi             = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
        bait_intervals        = CAT_CAT_BAIT.out.file_out.map{ meta, inter -> inter}.collect()   // channel: [ path(intervals) ]
        gnomad_af_idx         = TABIX_GNOMAD_AF.out.tbi.collect()                                // channel: [ val(meta), path(fasta) ]
        vep_resources         = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()  // channel: [ path(cache) ]
        mt_intervals          = ch_shiftfasta_mtintervals.intervals.collect()                    // channel: [ path(intervals) ]
        mtshift_intervals     = ch_shiftfasta_mtintervals.shift_intervals.collect()              // channel: [ path(intervals) ]
        mtshift_backchain     = GATK_SHIFTFASTA.out.shift_back_chain.collect()                   // channel: [ val(meta), path(backchain) ]
        mtshift_fai           = GATK_SHIFTFASTA.out.shift_fai.collect()                          // channel: [ val(meta), path(fai) ]
        mtshift_fasta         = GATK_SHIFTFASTA.out.shift_fa.collect()                           // channel: [ val(meta), path(fai) ]
        mtshift_dict          = GATK_SHIFTFASTA.out.dict.collect()                               // channel: [ path(dict) ]
        mtshift_bwamem2_index = BWAMEM2_INDEX_MT_SHIFT.out.index.collect()                       // channel: [ val(meta), path(index) ]
        versions              = ch_versions                                                      // channel: [ path(versions.yml) ]

}
