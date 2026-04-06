#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING } from './subworkflows/preprocessing/main.nf'
include { SPLIT_ALIGN_CLEAN_CHUNKS } from './subworkflows/split_align/main.nf'
include { XORF as XORF_PREDICT_ORFS } from '../modules/xorf/src/subworkflows/xorf/main.nf'
include { XORF as XORF_PREDICT_FUSION_ORFS } from '../modules/xorf/src/subworkflows/xorf/main.nf'
include { ISOTOOLS_NMD as ISOTOOLS_NMD_FILTER } from './modules/custom/isotools/nmd/main.nf'
include { PREPOLISH as ISOTOOLS_PREPOLISH } from './subworkflows/prepolish/main.nf'
include { POLISH as ISOTOOLS_POLISH } from './subworkflows/polish/main.nf'

include { LOAD_TRACK as LOAD_PASS_TRACK } from './subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_TRASH_TRACK } from './subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_RETENTIONS_TRACK } from './subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_TRUNCATIONS_TRACK } from './subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_INTRAPRIMMING_TRACK } from './subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_FUSIONS_TRACK } from './subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_NMD_TRACK } from './subworkflows/track/main.nf'

include { BEDTOBIGBED as BEDTOBIGBED_FUSIONS } from './modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_NMD } from './modules/custom/bigtools/bedtobigbed/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISOPIPE {
    main: 
      ch_versions = Channel.empty()
      ch_reads = Channel.empty()

      PREPROCESSING(
        params.entrypoint,
        params.global_input_dir,
        params.global_primers,
        params.global_genome,
        params.global_annotation,
        params.ccs_chunk,
        params.isoseq_cluster2_mode,
        params.xorf_protein_database,
        params.minimap2_index_path,
        params.minimap2_align_use_splice_scores,
        params.minimap2_align_splicing_algorithm,
        params.spliceai_bigwigs_dir,
        params.minisplice_scores_path,
        params.spliceai_scores_path,
        ch_versions
      )

      SPLIT_ALIGN_CLEAN_CHUNKS(
        PREPROCESSING.out.reads, 
        PREPROCESSING.out.minimap2_index,
        PREPROCESSING.out.reference_transcripts,
        PREPROCESSING.out.splice_scores,
        params.isoseq_cluster2_mode,
        params.entrypoint,
        params.minimap2_align_use_junc_bed,
        ch_versions
      )

      XORF_PREDICT_ORFS(
          SPLIT_ALIGN_CLEAN_CHUNKS.out.reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.database,
          params.global_output_dir,
          params.xorf_chunk_size,
          params.xorf_samba_weights,
          params.xorf_predict_keep_raw
      )
      XORF_PREDICT_ORFS.out.files
          .map { meta, bed, tsv -> [ meta, bed ] }
          .set { ch_orf_predictions_bed }
      
      XORF_PREDICT_FUSION_ORFS(
          SPLIT_ALIGN_CLEAN_CHUNKS.out.fusions,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.database,
          params.global_output_dir,
          params.xorf_chunk_size,
          params.xorf_samba_weights,
          params.xorf_predict_keep_raw
      )
      XORF_PREDICT_FUSION_ORFS.out.files
          .map { meta, bed, tsv -> [ meta, bed ] }
          .set { ch_fusion_orf_predictions_bed }

      ISOTOOLS_NMD_FILTER(ch_orf_predictions_bed)
      ISOTOOLS_PREPOLISH(
          ISOTOOLS_NMD_FILTER.out.reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.chrom_sizes,
          params.global_repeats,
          PREPROCESSING.out.reference_transcripts,
          PREPROCESSING.out.bigwigs,
          params.aparent_predict_weights,
          ch_versions
      )

      ISOTOOLS_NMD_FILTER.out.reads
        .map { meta, read -> tuple(meta.id, meta, read) }
        .join(
          ISOTOOLS_PREPOLISH.out.introns
            .map { meta, introns -> tuple(meta.id, introns) }
        )
        .map { id, meta, read, introns ->
          tuple(id, meta, read, introns)
        }
        .join(
          XORF_PREDICT_ORFS.out.files
              .map { meta, bed, tsv -> tuple(meta.id, tsv) }
        )
        .map { id, meta, read, introns, tsv ->
          tuple(meta, read, introns, tsv)
        }
        .set { ch_full_length_reads }

      ISOTOOLS_POLISH(
          ch_full_length_reads,
          PREPROCESSING.out.reference_transcripts,
          ISOTOOLS_PREPOLISH.out.aparent_plus,
          ISOTOOLS_PREPOLISH.out.aparent_minus,
          PREPROCESSING.out.chrom_sizes,
          ch_versions
      )

      autosql = Channel.value(file('${projectDir}/../../assets/as/base.as', checkIfExists: true))
      BEDTOBIGBED_FUSIONS(ch_fusion_orf_predictions_bed, PREPROCESSING.out.chrom_sizes, autosql)
      BEDTOBIGBED_NMD(ISOTOOLS_NMD_FILTER.out.nmd, PREPROCESSING.out.chrom_sizes, autosql)

      // INFO: ch_fusion_orf_predictions_bed + ISOTOOLS_NMD_FILTER.out.nmd
      if (params.load_track) {
          LOAD_PASS_TRACK(
            ISOTOOLS_POLISH.out.pass,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            ch_versions
          )

          LOAD_TRASH_TRACK(
            ISOTOOLS_POLISH.out.trash,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            ch_versions
          )

          LOAD_RETENTIONS_TRACK(
            ISOTOOLS_POLISH.out.retentions,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            ch_versions
          )

          LOAD_TRUNCATIONS_TRACK(
            ISOTOOLS_POLISH.out.truncations,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            ch_versions
          )

          LOAD_INTRAPRIMMING_TRACK(
            ISOTOOLS_POLISH.out.intraprimming,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            ch_versions
          )

          LOAD_FUSIONS_TRACK(
            BEDTOBIGBED_FUSIONS.out.bigbed,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            ch_versions
          )

          LOAD_NMD_TRACK(
            ISOTOOLS_NMD_FILTER.out.nmd,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            ch_versions
          )
      }
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
