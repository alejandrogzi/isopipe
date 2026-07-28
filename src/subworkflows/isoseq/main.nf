/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PBCCS } from '../../modules/nf-core/pbccs/main.nf'
include { PBTK_PBINDEX as PBINDEX } from '../../modules/nf-core/pbtk/pbindex/main.nf'
include { PBTK_PBMERGE as PBMERGE } from '../../modules/nf-core/pbtk/pbmerge/main.nf'
include { PBTK_PBMERGE as PBMERGE_MULTI_SAMPLE } from '../../modules/nf-core/pbtk/pbmerge/main.nf'
include { PBTK_PBMERGE as PBMERGE_MULTI_LIMA } from '../../modules/nf-core/pbtk/pbmerge/main.nf'
include { LIMA } from '../../modules/nf-core/lima/main.nf'
include { ISOSEQ_REFINE } from '../../modules/nf-core/isoseq/refine/main.nf'
include { ISOSEQ_CLUSTER2 } from '../../modules/custom/isoseq/cluster2/main.nf'
include { ISOSEQ_CLUSTER2 as ISOSEQ_CLUSTER2_MULTI_SAMPLE } from '../../modules/custom/isoseq/cluster2/main.nf'
include { BAM_TO_FA } from '../../modules/custom/bamtofa/main.nf'
include { BAM_TO_FA as BAM_TO_FA_MULTI_SAMPLE } from '../../modules/custom/bamtofa/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
workflow ISOSEQ {
    take:
      global_input_dir       // path
      global_primers         // path
      ccs_chunk              // int
      isoseq_cluster2_mode   // string
      prefix                 // string
      entrypoint             // [ subreads, ccs ] (flnc unreachable)

    main:
      ch_versions = Channel.empty()
      ch_primers = Channel.value(file(global_primers))

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHANNELING/INDEXING [ SUBREADS, CCS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      Channel
          .fromPath("${global_input_dir}/*.bam", checkIfExists: true)
          .map { bam ->
              def pbi = bam + '.pbi'
              return [
                  [
                      id:         bam.baseName,
                      single_end: true,
                      indexed:    file(pbi).exists()
                  ],
                  bam,
                  file(pbi).exists() ? file(pbi) : []   // placeholder if missing
              ]
          }
          .set { ch_bam }


      ch_bam
          .branch {
              indexed:     it[0].indexed
              not_indexed: true
          }
          .set { ch_bam_branched }


      PBINDEX(
        ch_bam_branched.not_indexed
        .map { 
          meta, bam, _pbi -> [ meta, bam ] 
        }
      )

       PBINDEX.out.pbi
          .join(
              ch_bam_branched.not_indexed.map { meta, bam, _pbi -> [ meta, bam ] },
              by: 0   // join on meta
          )
          .map { meta, pbi, bam ->
              def meta_updated = meta + [ indexed: true ]
              [ meta_updated, bam, pbi ]
          }.set { ch_bam_reindexed }

      ch_bam = ch_bam_branched.indexed.mix(ch_bam_reindexed)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ENTRYPOINT BRANCHING  [ SUBREADS, CCS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */
  
      ch_lima_bams = Channel.empty()
      switch (entrypoint) {
        case 'subreads':

          ch_bam
              .combine(Channel.of(1..ccs_chunk))   // INFO: cartesian product: N_bam × chunk combos
              .map { meta, bam, pbi, chunk_idx ->
                  [ meta + [chunk: chunk_idx], bam, pbi ]
              }
              .set { ch_chunks }

          PBCCS(ch_chunks, ccs_chunk) // INFO: generate CCS from raw reads
          PBCCS.out.bam // INFO: update meta: update id (+chunkX) and store former id
          .map {
              def chunk   = it[0].chunk
              def parent  = it[0].id
              def child   = it[0].id + "." + chunk
              return [ [id:child, parent:parent, single_end:true, chunk:chunk], it[1] ]
          }
          .set { ch_pbccs_bam_updated }

          // INFO: group all chunks belonging to the same parent sample
          ch_pbccs_bam_updated
              .map { meta, bam -> [ meta.parent, meta, bam ] }   // INFO: key by parent
              .groupTuple(by: 0, size: ccs_chunk)              // INFO: wait for all chunks
              .map { parent, metas, bams ->
                  // Reconstruct a clean meta for the merged output
                  def meta_merged = [ id: parent, single_end: true ]
                  [ meta_merged, bams ]                           // INFO: bams is now a List
              }
              .set { ch_pbccs_merged }

          PBMERGE(ch_pbccs_merged) // INFO: merge chunks
          ch_lima_bams = PBMERGE.out.bam

        break

        case 'ccs':
          ch_lima_bams = ch_bam
        break

        default:
          error """
          ERROR: Unknown entrypoint -> options are at this step are: ccs, flnc
          """.stripIndent()
          System.exit(1)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          LIMA [ REMOVE PRIMERS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      LIMA(ch_lima_bams, ch_primers)  // INFO: remove primers from CCS
      LIMA.out.bam
          .map { meta, bams ->
              [ meta, bams instanceof List ? bams : [bams] ]
          }
          .branch {
              meta, bams ->
              merge:  bams.size() > 1
              single: bams.size() == 1
          }
          .set { ch_lima_branched }

      PBMERGE_MULTI_LIMA(ch_lima_branched.merge)
      ch_lima_bams = ch_lima_bams.mix(PBMERGE_MULTI_LIMA.out.bam)
      ch_lima_bams = ch_lima_bams.mix(ch_lima_branched.single.map { meta, bams -> [ meta, bams[0] ] })

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ISOSEQ_REFINE [ DISCARD CCS WITHOUT POLYA TAILS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOSEQ_REFINE(ch_lima_bams, ch_primers) // INFO: discard CCS without polyA tails, remove it from the other

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ISOSEQ_CLUSTER2 [ CLUSTER READS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: three possible modes: per_tissue, pan_tissue, both
      ch_pbccs_merged_flnc_clustered_fa = Channel.empty()
      if (isoseq_cluster2_mode == "per_sample") {
        ISOSEQ_CLUSTER2(ISOSEQ_REFINE.out.bam) // INFO: cluster reads
        BAM_TO_FA(ISOSEQ_CLUSTER2.out.bam)

        ch_pbccs_merged_flnc_clustered_fa  = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.singletons)
        ch_pbccs_merged_flnc_clustered_fa  = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.hq)

        ch_versions = ch_versions.mix(ISOSEQ_CLUSTER2.out.versions)
        ch_versions = ch_versions.mix(BAM_TO_FA.out.versions)
      } else if (isoseq_cluster2_mode == "multi_sample") {
        // INFO: collect all refined reads into a single channel for clustering
        ISOSEQ_REFINE.out.bam
          .map { meta, bam -> bam  }
          .collect()
          .map { bams -> [ [ id: prefix ], bams ] }
          .set { ch_pooled_bams }

        PBMERGE_MULTI_SAMPLE(ch_pooled_bams)
        ISOSEQ_CLUSTER2_MULTI_SAMPLE(PBMERGE_MULTI_SAMPLE.out.bam)
        BAM_TO_FA_MULTI_SAMPLE(ISOSEQ_CLUSTER2_MULTI_SAMPLE.out.bam)

        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA_MULTI_SAMPLE.out.singletons)
        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA_MULTI_SAMPLE.out.hq)

        ch_versions = ch_versions.mix(BAM_TO_FA_MULTI_SAMPLE.out.versions)
        ch_versions = ch_versions.mix(ISOSEQ_CLUSTER2_MULTI_SAMPLE.out.versions)
        ch_versions = ch_versions.mix(PBMERGE_MULTI_SAMPLE.out.versions)
      } else if (isoseq_cluster2_mode == "both") {
        ISOSEQ_CLUSTER2(ISOSEQ_REFINE.out.bam)
        BAM_TO_FA(ISOSEQ_CLUSTER2.out.bam)
        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.singletons)
        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.hq)

        ISOSEQ_REFINE.out.bam
          .map { meta, bam -> bam  }
          .collect()
          .map { bams -> [ [id: prefix, single_end:true], bams ] }
          .set { ch_pooled_bams }

        PBMERGE_MULTI_SAMPLE(ch_pooled_bams)
        ISOSEQ_CLUSTER2_MULTI_SAMPLE(PBMERGE_MULTI_SAMPLE.out.bam)
        BAM_TO_FA_MULTI_SAMPLE(ISOSEQ_CLUSTER2_MULTI_SAMPLE.out.bam)
        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA_MULTI_SAMPLE.out.singletons)
        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA_MULTI_SAMPLE.out.hq)

        ch_versions = ch_versions.mix(ISOSEQ_CLUSTER2.out.versions)
        ch_versions = ch_versions.mix(BAM_TO_FA_MULTI_SAMPLE.out.versions)
        ch_versions = ch_versions.mix(BAM_TO_FA.out.versions)
        ch_versions = ch_versions.mix(ISOSEQ_CLUSTER2_MULTI_SAMPLE.out.versions)
        ch_versions = ch_versions.mix(PBMERGE_MULTI_SAMPLE.out.versions)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VERSIONING
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_versions = ch_versions.mix(ISOSEQ_REFINE.out.versions)
      ch_versions = ch_versions.mix(LIMA.out.versions)
      ch_versions = ch_versions.mix(PBMERGE.out.versions)
      ch_versions = ch_versions.mix(PBINDEX.out.versions)
      ch_versions = ch_versions.mix(PBCCS.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
        reads   = ch_pbccs_merged_flnc_clustered_fa
        versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
