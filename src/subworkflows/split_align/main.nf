/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN } from '../../modules/custom/minimap2/align/main.nf'
include { FXSPLIT } from '../../modules/custom/fxsplit/main.nf'
include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_FUSION as ISOTOOLS_FUSION_DETECTOR } from '../../modules/custom/isotools/fusion/main.nf'
include { SAMTOOLS_BAM } from '../../modules/custom/samtools/bam/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE } from '../../modules/custom/samtools/merge/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPLIT_ALIGN_CLEAN_CHUNKS {
    take:
      ch_reads                 // [ meta, reads ]
      ch_minimap2_index        // [ meta, index ]
      ch_reference_transcripts // [ file ]
      ch_splice_scores         // [ meta, scores ]
      cluster_mode             // string [ per_sample, multi_sample, both ]
      entrypoint               // string [ isoseq, map ]
      minimap2_use_junc_bed    // bool
      ch_versions              // [ meta, versions.yml ]

    main:
      FXSPLIT(ch_reads)
      FXSPLIT.out.fastx_gz
          .flatMap { 
              meta, fa ->
              def fas = fa instanceof List ? fa : [fa]
              fas.collect { it ->
                  def parts = it.baseName.split('_')
                  def chunk = parts.size() > 2 ? parts[2] : 0
                  [ meta + [ chunk: chunk ], it ]
              }
          }
          .set { ch_fastx_gz }
 
      if (minimap2_use_junc_bed) {
          MINIMAP2_ALIGN(
            ch_fastx_gz, 
            ch_minimap2_index,
            ch_splice_scores,
            ch_reference_transcripts
              .map { 
                junctions -> 
                [ [id:junctions.baseName], junctions ] 
              }
          )
      } else {
          MINIMAP2_ALIGN(
            ch_fastx_gz, 
            ch_minimap2_index,
            ch_splice_scores,
            Channel.value([[:], []])
          )
      }

      ch_aligned_bam = Channel.empty()
      ch_aligned_bai = Channel.empty()
      if (entrypoint == "isoseq") {
        // INFO: reads were already merged in isoseq clustering
        SAMTOOLS_BAM(MINIMAP2_ALIGN.out.sam)
        ch_aligned_bam = SAMTOOLS_BAM.out.bam
        ch_aligned_bai = SAMTOOLS_BAM.out.bai

        ch_versions = ch_versions.mix(SAMTOOLS_BAM.out.versions)
      } else if (entrypoint == "map") {
        // INFO: fastq reads need to be merged if multi-sample
        if (cluster_mode == "per_sample") {
          SAMTOOLS_BAM(MINIMAP2_ALIGN.out.sam)
          ch_aligned_bam = SAMTOOLS_BAM.out.bam
          ch_aligned_bai = SAMTOOLS_BAM.out.bai

          ch_versions = ch_versions.mix(SAMTOOLS_BAM.out.versions)
        } else if (cluster_mode == "multi_sample") {
          MINIMAP2_ALIGN.out.sam
              .map { meta, sam -> sam }
              .collect()
              .map { sams -> [ 
                [
                  id: 'pooled', 
                  single_end: false,
                  singleton: false,
                  chunk: 0
                ], sams ] }
              .set { ch_joined_sam }

          SAMTOOLS_MERGE_BAM_MULTI_SAMPLE(ch_joined_sam)
          ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bam
          ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bai

          ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.versions)
        } else if (cluster_mode == "both") {
          SAMTOOLS_BAM(MINIMAP2_ALIGN.out.sam)
          ch_aligned_bam = SAMTOOLS_BAM.out.bam
          ch_aligned_bai = SAMTOOLS_BAM.out.bai

          MINIMAP2_ALIGN.out.sam
              .map { meta, sam -> sam }
              .collect()
              .map { sams -> [ [id:'pooled'], sams ] }
              .set { ch_joined_sam }
          SAMTOOLS_MERGE_BAM_MULTI_SAMPLE(ch_joined_sam)

          ch_aligned_bam = ch_aligned_bam.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bam)
          ch_aligned_bai = ch_aligned_bai.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bai)

          ch_versions = ch_versions.mix(SAMTOOLS_BAM.out.versions)
          ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.versions)
        }
      }

      ISOTOOLS_SEGMENT_POLYA(ch_aligned_bam, ch_aligned_bai) // INFO: polyA tails
      ISOTOOLS_SEGMENT_POLYA.out.hq_bed
          .flatMap { 
            meta, bed -> { 
              def beds = bed instanceof List ? bed : [bed]
              beds.collect { it ->
                  [ meta + [ chr: it.name.split('@')[0] ], it ]
              }
            } 
          }
          .map { meta, bed -> [ [ meta.chr, meta.id ], meta, bed ] }   // INFO: key by parent
          .groupTuple(by: 0)
          .map { chr, metas, beds ->
              // Reconstruct a clean meta for the merged output
              def meta = metas[0]
              def group_meta = [ 
                    id: meta.id, 
                    single_end: true, 
                    chr: meta.chr
                ]
              [ group_meta, beds ]                           
          }
          .set { ch_hq_bed_per_chr }

      ISOTOOLS_FUSION_DETECTOR(ch_hq_bed_per_chr, ch_reference_transcripts)

      ch_versions = ch_versions.mix(FXSPLIT.out.versions)
      ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_SEGMENT_POLYA.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_FUSION_DETECTOR.out.versions)

    emit:
      reads    = ISOTOOLS_FUSION_DETECTOR.out.free_fusion
      fusions  = ISOTOOLS_FUSION_DETECTOR.out.fusion
      versions = ch_versions
}
