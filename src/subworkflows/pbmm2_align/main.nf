/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PBMM2_ALIGN } from '../../modules/nf-core/pbmm2/align/main.nf'

include { FXSPLIT } from '../../modules/custom/fxsplit/main.nf'

include { SAMTOOLS_INDEX as SAMTOOLS_BAM_INDEX } from '../../modules/custom/samtools/index/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE } from '../../modules/custom/samtools/merge/main.nf'

include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_FUSION as ISOTOOLS_FUSION_DETECTOR } from '../../modules/custom/isotools/fusion/main.nf'
include { ISOTOOLS_CIGAR as ISOTOOLS_CIGAR_EXTENSION } from '../../modules/custom/isotools/cigar/main.nf'
include { ISOTOOLS_ADAPTER as ISOTOOLS_REMOVE_ADAPTERS } from '../../modules/custom/isotools/adapter/main.nf'

include { COLLAPSE as COLLAPSE_TWINS } from '../../modules/custom/collapse/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPLIT_PBMM2_ALIGN_CLEAN_CHUNKS {
    take:
      ch_reads                 // [ meta, reads ]
      ch_genome                // [ genome ]
      ch_genome_index          // [ meta, index ]
      ch_reference_transcripts // [ meta, bed ]
      ch_splice_scores         // [ meta, scores ]
      cluster_mode             // string [ per_sample, multi_sample, both ]
      entrypoint               // string [ isoseq, map ]
      remove_adapters          // bool
      prefix                   // string
      collapse_twins           // bool
      cigar_extension          // bool
      ch_versions              // [ meta, versions.yml ]

    main:
      // INFO: Chunking step ///////////////////////////////////

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

      // INFO: Alignment step ///////////////////////////////////
      
      PBMM2_ALIGN(
        ch_fastx_gz,
        ch_genome_index
      )

      // INFO: Clustering and branching step ///////////////////

      ch_aligned_bam = Channel.empty()
      ch_aligned_bai = Channel.empty()
      ch_pooled_reads = Channel.empty()
      if (entrypoint == "isoseq") {
        // INFO: reads were already merged in isoseq clustering
        SAMTOOLS_BAM_INDEX(PBMM2_ALIGN.out.bam)
        ch_aligned_bam = SAMTOOLS_BAM_INDEX.out.bam
        ch_aligned_bai = SAMTOOLS_BAM_INDEX.out.bai
        ch_pooled_reads = ch_reads

        ch_versions = ch_versions.mix(SAMTOOLS_BAM_INDEX.out.versions)
      } else if (entrypoint == "map") {
        // INFO: fastq reads need to be merged if multi-sample
        if (cluster_mode == "per_sample") {
          SAMTOOLS_BAM_INDEX(PBMM2_ALIGN.out.bam)
          ch_aligned_bam = SAMTOOLS_BAM_INDEX.out.bam
          ch_aligned_bai = SAMTOOLS_BAM_INDEX.out.bai
          ch_pooled_reads = ch_reads

          ch_versions = ch_versions.mix(SAMTOOLS_BAM_INDEX.out.versions)
        } else if (cluster_mode == "multi_sample") {
          SAMTOOLS_BAM_INDEX(PBMM2_ALIGN.out.bam)
          SAMTOOLS_BAM_INDEX.out.bam
              .map { meta, bam -> bam }
              .collect()
              .map { bams -> [
                [
                  id: prefix,
                  single_end: false,
                  singleton: false,
                  chunk: 0
                ], bams ] }
              .set { ch_joined_bam }

          SAMTOOLS_MERGE_BAM_MULTI_SAMPLE(ch_joined_bam)
          ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bam
          ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bai

          ch_reads
            .map { meta, reads -> reads }
            .collect()
            .map { reads -> [ [ id: 'pooled.reads' ], reads ] }
            .set { ch_pooled_reads }

          ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.versions)
        } else if (cluster_mode == "both") {
          SAMTOOLS_BAM_INDEX(PBMM2_ALIGN.out.bam)
          ch_aligned_bam = SAMTOOLS_BAM_INDEX.out.bam
          ch_aligned_bai = SAMTOOLS_BAM_INDEX.out.bai

          SAMTOOLS_BAM_INDEX.out.bam
              .map { meta, bam -> bam }
              .collect()
              .map { bams -> [ [ id: prefix ], bams ] }
              .set { ch_joined_bam }
          SAMTOOLS_MERGE_BAM_MULTI_SAMPLE(ch_joined_bam)

          ch_aligned_bam = ch_aligned_bam.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bam)
          ch_aligned_bai = ch_aligned_bai.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.bai)

          ch_reads
            .map { meta, reads -> reads }
            .collect()
            .map { reads -> [ [ id: 'pooled.reads' ], reads ] }
            .set { ch_pooled_reads }

          ch_versions = ch_versions.mix(SAMTOOLS_BAM_INDEX.out.versions)
          ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE.out.versions)
        }
      }

      // INFO: Adapter removal step ///////////////////////////

      if (remove_adapters) {
        ISOTOOLS_REMOVE_ADAPTERS(
          ch_aligned_bam,
          ch_aligned_bai,
        )
        ch_aligned_bam = ISOTOOLS_REMOVE_ADAPTERS.out.bam
        ch_aligned_bai = ISOTOOLS_REMOVE_ADAPTERS.out.bai
        ch_versions = ch_versions.mix(ISOTOOLS_REMOVE_ADAPTERS.out.versions)
      }

      // INFO: Cigar extension + fragment detection step //////

      if (cigar_extension) {
        ISOTOOLS_CIGAR_EXTENSION(
          ch_aligned_bam,
          ch_aligned_bai,
          ch_genome.map { genome -> [ [id:genome.baseName], genome ] },
          ch_reference_transcripts
        )

        // WARN: replacing ch_aligned_bam with [ meta, bam, bai ]
        ch_aligned_bam = ISOTOOLS_CIGAR_EXTENSION.out.extended
        ch_versions = ch_versions.mix(ISOTOOLS_CIGAR_EXTENSION.out.versions)
      } 

      // INFO: Segmentation step ///////////////////////////////////

      ISOTOOLS_SEGMENT_POLYA(ch_aligned_bam) // INFO: polyA tails + cigar 
      ISOTOOLS_SEGMENT_POLYA.out.hq_bed
          .set { ch_aligned_segmented }

      ch_aligned_segmented
          .flatMap { meta, bed ->   // ← no extra braces
              def beds = bed instanceof List ? bed : [bed]
              beds.collect { it ->
                  [ meta + [ chr: it.name.split('@')[0] ], it ]
              }
          }
          // INFO: fmt -> chr1@pooled.hq.bed OR chr1@pooled.extended.hq.bed
          .map { meta, bed ->
              def sampleId = bed.name.split('@')[1].split('\\.')[0]  // → "pooled"
              [ [ meta.chr, sampleId ], meta, bed ]
          }
          .groupTuple(by: 0)
          .map { chr, metas, beds ->
              def meta = metas[0]
              def sampleId = beds[0].name.split('@')[1].split('\\.')[0]  // → "pooled"
              def group_meta = [ id: sampleId, single_end: true, chr: meta.chr ]
              [ group_meta, beds ]
          }
          .set { ch_aligned_segmented_hq_per_chr }

      // INFO: Collapse + fusion detection step ///////////////////////////

      ch_aligned_segmented_collapsed = Channel.empty()
      if (collapse_twins) {
        COLLAPSE_TWINS(ch_aligned_segmented_hq_per_chr)
        ch_aligned_segmented_collapsed = COLLAPSE_TWINS.out.collapsed
      } else {
        ch_aligned_segmented_collapsed = ch_aligned_segmented_hq_per_chr
      }

      ISOTOOLS_FUSION_DETECTOR(
        ch_aligned_segmented_collapsed, 
        ch_reference_transcripts
      )

      ch_versions = ch_versions.mix(FXSPLIT.out.versions)
      ch_versions = ch_versions.mix(PBMM2_ALIGN.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_SEGMENT_POLYA.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_FUSION_DETECTOR.out.versions)

    emit:
      reads    = ISOTOOLS_FUSION_DETECTOR.out.free_fusion
      fusions  = ISOTOOLS_FUSION_DETECTOR.out.fusion
      versions = ch_versions
}
