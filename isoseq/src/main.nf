#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PBCCS } from './modules/nf-core/pbccs/main.nf'
include { PBTK_PBINDEX as PBINDEX } from './modules/nf-core/pbtk/pbindex/main.nf'
include { PBTK_PBMERGE as PBMERGE } from './modules/nf-core/pbtk/pbmerge/main.nf'
include { LIMA } from './modules/nf-core/lima/main.nf'
include { ISOSEQ_REFINE } from './modules/nf-core/isoseq/refine/main.nf'
include { ISOSEQ_CLUSTER2 } from './modules/custom/isoseq/cluster2/main.nf'
include { ISOSEQ_CLUSTER2 as ISOSEQ_CLUSTER2_POOLED } from './modules/custom/isoseq/cluster2/main.nf'
include { BAM_TO_FQ } from './modules/custom/bamtofq/main.nf'
include { BAM_TO_FQ as BAM_TO_FQ_POOLED } from './modules/custom/bamtofq/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISOPIPE {
    main: 

    ch_versions = Channel.empty()
    // ch_genome = PREPARE_GENOME(params.genome)
    // ch_minimap2_index = MINIMAP2_INDEX(ch_genome)

    // INFO: isoseq entrypoint
    if (params.entrypoint == "isoseq") {
      ch_primers = Channel.value(file(params.primers))

      Channel
          .fromPath("${params.input_dir}/*.bam", checkIfExists: true)
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

      ch_bam
          .combine(Channel.of(1..params.chunk))   // INFO: cartesian product: N_bam × chunk combos
          .map { meta, bam, pbi, chunk_idx ->
              [ meta + [chunk: chunk_idx], bam, pbi ]
          }
          .set { ch_chunks }

      PBCCS(ch_chunks, params.chunk) // INFO: generate CCS from raw reads
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
          .groupTuple(by: 0, size: params.chunk)              // INFO: wait for all chunks
          .map { parent, metas, bams ->
              // Reconstruct a clean meta for the merged output
              def meta_merged = [ id: parent, single_end: true ]
              [ meta_merged, bams ]                           // INFO: bams is now a List
          }
          .set { ch_pbccs_merged }

      PBMERGE(ch_pbccs_merged) // INFO: merge chunks
      LIMA(PBMERGE.out.bam, ch_primers)  // INFO: remove primers from CCS
      ISOSEQ_REFINE(LIMA.out.bam, ch_primers) // INFO: discard CCS without polyA tails, remove it from the other

      // INFO: three possible modes: per_tissue, pan_tissue, both
      ch_pbccs_merged_flnc_clustered_fqs = Channel.empty()
      if (params.cluster_mode == "per_tissue") {
        ISOSEQ_CLUSTER2(ISOSEQ_REFINE.out.bam) // INFO: cluster reads
        BAM_TO_FQ(ISOSEQ_CLUSTER2.out.bam)
      } else if (params.cluster_mode == "pan_tissue") {
        // INFO: collect all refined reads into a single channel for clustering
        ISOSEQ_REFINE.out.bam
          .map { meta, bam -> bam  }
          .collect()
          //.flatten()
          .map { bams -> [ [id:'pooled'], bams ] }
          .set { ch_pooled_bams }
        ISOSEQ_CLUSTER2_POOLED(ch_pooled_bams)
        BAM_TO_FQ(ISOSEQ_CLUSTER2_POOLED.out.bam)
      } else if (params.cluster_mode == "both") {
        ISOSEQ_CLUSTER2(ISOSEQ_REFINE.out.bam)
        BAM_TO_FQ(ISOSEQ_CLUSTER2.out.bam)

        ISOSEQ_REFINE.out.bam
          .map { meta, bam -> bam  }
          .collect()
          //.flatten()
          .map { bams -> [ [id:'pooled'], bams ] }
          .set { ch_pooled_bams }
        ISOSEQ_CLUSTER2_POOLED(ch_pooled_bams)
        BAM_TO_FQ_POOLED(ISOSEQ_CLUSTER2_POOLED.out.bam)
      }

      // INFO: should end up emitting a channel with fastqs (hq and singletons in this branch)
    }

      //MINIMAP2_ALIGN(BAM_TO_FQ.out.singletons, ch_genome) // INFO: align singleton reads to reference
      //MINIMAP2_ALIGN(BAM_TO_FQ.out.hq, ch_genome) // INFO: align hq reads to reference
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    ISOPIPE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
