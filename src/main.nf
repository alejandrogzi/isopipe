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
          params.bigtools_bedtobigbed_autosql,
          PREPROCESSING.out.chrom_sizes,
          ch_versions
      )

      // INFO: ch_fusion_orf_predictions_bed + ISOTOOLS_NMD_FILTER.out.nmd
      // LOAD_TRACKS()
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
