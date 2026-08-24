#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING } from '../subworkflows/preprocessing/main.nf'

include { SPLIT_ALIGN_CLEAN_CHUNKS } from '../subworkflows/split_align/main.nf'

include { PREPOLISH as ISOTOOLS_PREPOLISH } from '../subworkflows/prepolish/main.nf'
include { POLISH as ISOTOOLS_POLISH } from '../subworkflows/polish/main.nf'

include { LOAD_TRACK as LOAD_PASS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_DUPLICATES_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_ORPHANS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_TRASH_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_RETENTIONS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_TRUNCATIONS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_INTRAPRIMMING_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_FUSIONS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_NMD_TRACK } from '../subworkflows/track/main.nf'

include { XORF as XORF_PREDICT_ORFS } from '../../modules/xorf/src/subworkflows/xorf/main.nf'
include { XORF as XORF_PREDICT_FUSION_ORFS } from '../../modules/xorf/src/subworkflows/xorf/main.nf'

include { WGET as WGET_APARENT_WEIGHTS } from '../modules/nf-core/wget/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from '../modules/nf-core/wget/main.nf'

include { ISOTOOLS_NMD as ISOTOOLS_NMD_FILTER } from '../modules/custom/isotools/nmd/main.nf'

include { GAWK_JOIN as JOIN_FUSIONS } from '../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_NMD } from '../modules/custom/gawk/join/main.nf'

include { BEDTOBIGBED as BEDTOBIGBED_FUSIONS } from '../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_NMD } from '../modules/custom/bigtools/bedtobigbed/main.nf'

include { PUBLISH as PUBLISH_ADDITIONAL_BIGBEDS } from '../modules/custom/publish/main.nf'
include { TRACKDB } from '../modules/custom/track/main.nf'

include { AUTOSQL_BASE } from '../modules/custom/autosql/base/main.nf'
include { AUTOSQL_SCHEMA } from '../modules/custom/autosql/schema/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ARK {
    main:
      ch_versions = Channel.empty()
      ch_reads = Channel.empty()

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          AUTOSQL
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      AUTOSQL_BASE()
      AUTOSQL_SCHEMA()

      autosql = AUTOSQL_BASE.out.autosql
      schema = AUTOSQL_SCHEMA.out.autosql
      track = Channel.value(file('${projectDir}/../../assets/as/track.as'))

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PRE-PROCESSING [ INDEXES, SPLICE SCORES, READS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      PREPROCESSING(
        params.entrypoint,
        params.global_input_dir,
        params.global_primers,
        params.global_genome,
        params.global_annotation,
        params.ccs_chunk,
        params.isoseq_cluster2_mode,
        params.xorf_protein_database,
        params.xorf_custom_database,
        params.xorf_raw_database,
        params.minimap2_index_path,
        params.minimap2_align_use_splice_scores,
        params.minimap2_align_splicing_algorithm,
        params.spliceai_bigwigs_dir,
        params.minisplice_scores_path,
        params.spliceai_scores_path,
        params.spliceai_chunk_compression,
        params.global_prefix,
        params.aligner,
        params.ultra_use_annotation,
        params.ultra_index,
        params.desalt_index,
        params.pbmm2_index,
        params.skera_is_kinnex_library,
        ch_versions
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ALIGNMENT [ SPLIT_ALIGN_CLEAN_CHUNKS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */
      
      if (params.aligner in ['mm2', 'ultra', 'pbmm2', 'desalt', 'ark', 'flair']) {
        SPLIT_ALIGN_CLEAN_CHUNKS(
          PREPROCESSING.out.reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.genome_index,
          PREPROCESSING.out.reference_transcripts,
          params.global_prefix,
          PREPROCESSING.out.splice_scores,
          params.isoseq_cluster2_mode,
          params.entrypoint,
          params.aligner,
          params.minimap2_align_use_junc_bed,
          params.isotools_adapter_remove_adapters,
          params.collapse_shrink_twins,
          params.isotools_cigar_extension_extend,
          params.minimap2_align_do_second_pass,
          ch_versions
        )

        ch_reads = SPLIT_ALIGN_CLEAN_CHUNKS.out.reads
        ch_fusions = SPLIT_ALIGN_CLEAN_CHUNKS.out.fusions
      } else {
        error """
        ERROR: Unsupported aligner: '${params.aligner}'.
        Valid aligners are: mm2, pbmm2, desalt, ultra, ark, flair
        """.stripIndent()
        System.exit(1)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          WEIGHTS [ XORF ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_samba_weights = Channel.empty()
      if (params.xorf_samba_local_weights) {
        ch_samba_weights = Channel.value(
          file(params.xorf_samba_local_weights, checkIfExists: true)
        ).map { path -> [ [id : path.baseName ], path ] }
      } else {
        WGET_SAMBA_WEIGHTS(
          Channel.value(
            params.xorf_samba_weights
          ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
        )
        ch_samba_weights = WGET_SAMBA_WEIGHTS.out.outfile
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ORF CALLING [ XORF ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // channel comes per chromosome previously splitted at segmentation step
      // metadata -> [ id: sampleId, single_end: true, chr: meta.chr ]
      XORF_PREDICT_ORFS(
          ch_reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.database,
          params.global_output_dir,
          params.xorf_chunk_size,
          ch_samba_weights,
          params.xorf_predict_keep_raw,
          params.xorf_selenocysteine_codons,
          params.xorf_skip_netstart,
          params.xorf_rename_deactivate,
          false, // xorf_do_polishing
          true,  // xorf_skip_joined_concat
          false, // xorf_run_only_on
          null,  // xorf_run_only_mode
          null,  // xorf_run_only_target
          Channel.empty() // xorf_database_versions
      )
      XORF_PREDICT_ORFS.out.files
          .map { meta, bed, tsv -> [ meta, bed ] }
          .set { ch_orf_predictions_bed }

      XORF_PREDICT_FUSION_ORFS(
          ch_fusions,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.database,
          params.global_output_dir,
          params.xorf_chunk_size,
          ch_samba_weights,
          params.xorf_predict_keep_raw,
          params.xorf_selenocysteine_codons,
          params.xorf_skip_netstart,
          params.xorf_rename_deactivate,
          false,
          true,
          false, // xorf_run_only_on
          null,  // xorf_run_only_mode
          null,  // xorf_run_only_target
          Channel.empty() // xorf_database_versions
      )
      XORF_PREDICT_FUSION_ORFS.out.files
          .map { meta, bed, tsv -> [ meta.name, meta, bed ] }
          .groupTuple()
          .map { name, metas, files ->
              [ [ id: name + '.fusions', name: name ], files ]
          }
          .set { ch_fusion_orf_predictions_bed }
      JOIN_FUSIONS(ch_fusion_orf_predictions_bed, 'bed')
      BEDTOBIGBED_FUSIONS(JOIN_FUSIONS.out.output, PREPROCESSING.out.chrom_sizes, autosql)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          NMD CALLING [ ISOTOOLS_NMD_FILTER ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOTOOLS_NMD_FILTER(ch_orf_predictions_bed)

      ISOTOOLS_NMD_FILTER.out.nmd
          .map { meta, bed -> [ meta.name, meta, bed ] }
          .groupTuple()
          .map { name, metas, files ->
              [ [ id: name + '.nmd', name: name ], files ]
          }
          .set { ch_nmd_bed }
      JOIN_NMD(ch_nmd_bed, 'bed')
      BEDTOBIGBED_NMD(JOIN_NMD.out.output, PREPROCESSING.out.chrom_sizes, autosql)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ADDITIONAL BIGBEDS [ BEDTOBIGBED_FUSIONS, BEDTOBIGBED_NMD ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_additional_bbs = Channel.empty()
      ch_additional_bbs = ch_additional_bbs.mix(BEDTOBIGBED_FUSIONS.out.bigbed)
      ch_additional_bbs = ch_additional_bbs.mix(BEDTOBIGBED_NMD.out.bigbed)
      ch_additional_bbs.map { meta, file -> [meta.name, meta, file] }
         .groupTuple()
         .map { name, metas, files -> [ [ id: name ], files] }
         .set { ch_additional_bbs }
      PUBLISH_ADDITIONAL_BIGBEDS(ch_additional_bbs)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          APARENT WEIGHTS [ POLYA PEAKS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_aparent_weights = Channel.empty()
      if (params.aparent_predict_weights_local_path) {
          ch_aparent_weights = Channel.value(
            file(params.aparent_predict_weights_local_path, checkIfExists: true)
          ).map { path -> [ [id : path.baseName ], path ] }
      } else {
          WGET_APARENT_WEIGHTS(
              Channel.value(
                params.aparent_predict_weights
              ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
          )
          ch_aparent_weights = WGET_APARENT_WEIGHTS.out.outfile
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PRE-POLISHING [ APARENT, INTRON, RETENTION, INTRAPRIMMING, BIGWIG  ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOTOOLS_PREPOLISH(
          ISOTOOLS_NMD_FILTER.out.reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.chrom_sizes,
          params.global_repeats,
          PREPROCESSING.out.reference_transcripts,
          PREPROCESSING.out.bigwigs,
          ch_aparent_weights,
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

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          POLISHING [ ISOTOOLS_POLISH ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOTOOLS_POLISH(
          ch_full_length_reads,
          PREPROCESSING.out.reference_transcripts,
          ISOTOOLS_PREPOLISH.out.aparent_plus,
          ISOTOOLS_PREPOLISH.out.aparent_minus,
          PREPROCESSING.out.chrom_sizes,
          PREPROCESSING.out.bigwigs,
          schema,
          ch_versions
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          TRACKING [ TRACKDB ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      if (params.load_track) {
          TRACKDB(
            track,
            params.load_track_browser,
            params.global_species_name,
            params.load_track_name,
            ISOTOOLS_POLISH.out.additional_columns,
            ISOTOOLS_POLISH.out.sample,
          )

          LOAD_PASS_TRACK(
            ISOTOOLS_POLISH.out.pass,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_DUPLICATES_TRACK(
            ISOTOOLS_POLISH.out.duplicates,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_ORPHANS_TRACK(
            ISOTOOLS_POLISH.out.orphans,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_TRASH_TRACK(
            ISOTOOLS_POLISH.out.trash,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_RETENTIONS_TRACK(
            ISOTOOLS_POLISH.out.retentions,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_TRUNCATIONS_TRACK(
            ISOTOOLS_POLISH.out.truncations,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_INTRAPRIMMING_TRACK(
            ISOTOOLS_POLISH.out.intraprimming,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_FUSIONS_TRACK(
            BEDTOBIGBED_FUSIONS.out.bigbed,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_NMD_TRACK(
            ISOTOOLS_NMD_FILTER.out.nmd,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          ch_versions = ch_versions.mix(TRACKDB.out.versions)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VERSIONING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */


      ch_versions = ch_versions.mix(ISOTOOLS_NMD_FILTER.out.versions)
      ch_versions = ch_versions.mix(JOIN_FUSIONS.out.versions)
      ch_versions = ch_versions.mix(JOIN_NMD.out.versions)
      ch_versions = ch_versions.mix(BEDTOBIGBED_FUSIONS.out.versions)
      ch_versions = ch_versions.mix(BEDTOBIGBED_NMD.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
