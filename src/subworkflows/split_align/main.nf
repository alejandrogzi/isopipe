/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN } from '../../modules/custom/minimap2/align/main.nf'
include { FXSPLIT } from '../../modules/custom/fxsplit/main.nf'
include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_FUSION as ISOTOOLS_FUSION_DETECTOR } from '../../modules/custom/isotools/fusion/main.nf'

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
 
      MINIMAP2_ALIGN(ch_fastx_gz, ch_minimap2_index) // INFO: align hq reads to reference
      ISOTOOLS_SEGMENT_POLYA(MINIMAP2_ALIGN.out.bam, MINIMAP2_ALIGN.out.bai) // INFO: polyA tails
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
