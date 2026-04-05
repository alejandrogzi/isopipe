/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOTOOLS_INTRON_RETENTION } from '../../modules/custom/isotools/intron/main.nf'
include { ISOTOOLS_PAS_CALLER } from '../../modules/custom/isotools/pas/main.nf'
include { ISOTOOLS_TRUNCATION_DETECTOR } from '../../modules/custom/isotools/utr/main.nf'
include { VEREDICT as ISOTOOLS_PLUGIN_VEREDICT } from '../../modules/custom/veredict/main.nf'

include { GAWK_JOIN as JOIN_VEREDICT_PASSES } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_TRASH } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_RETENTIONS } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_TRUNCATIONS } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_INTRAPRIMMING } from '../../modules/custom/gawk/join/main.nf'

include { BEDTOBIGBED as BEDTOBIGBED_PASSES } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_TRASH } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_RETENTIONS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_TRUNCATIONS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_INTRAPRIMMING } from '../../modules/custom/bigtools/bedtobigbed/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POLISH {
    take:
      reads                  // channel: [ val(meta), [ reads ], [ introns ], [ orfs ] ]
      annotation             // Channel.value(path)
      forward_peaks          // channel: [ val(meta), [ forward_peaks ] ]
      reverse_peaks          // channel: [ val(meta), [ reverse_peaks ] ]
      autosql                // path
      chrom_sizes            // path
      ch_versions            // [ meta, versions.yml ]

    main:
      ch_reference_transcripts = annotation.map { annotation -> [ [id:annotation.baseName], annotation ] }

      ISOTOOLS_INTRON_RETENTION(
        reads
      )

      ISOTOOLS_PAS_CALLER(
        reads,
        ch_reference_transcripts,
        forward_peaks,
        reverse_peaks
      )

      ISOTOOLS_TRUNCATION_DETECTOR(
        reads
      )

      reads
        .join(ISOTOOLS_INTRON_RETENTION.out.descriptor)
        .join(ISOTOOLS_PAS_CALLER.out.descriptor)
        .join(ISOTOOLS_TRUNCATION_DETECTOR.out.descriptor)
        .set { ch_polished_reads }

      ISOTOOLS_PLUGIN_VEREDICT(
        ch_polished_reads
      )

      // INFO: collect all veredict results per category and join
      ISOTOOLS_PLUGIN_VEREDICT.out.pass
        .collect()
        .set { ch_pass }
      JOIN_VEREDICT_PASSES(ch_pass, 'bed')
      BEDTOBIGBED_PASSES(JOIN_VEREDICT_PASSES.out.output, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.trash
        .collect()
        .set { ch_trash }
      JOIN_VEREDICT_TRASH(ch_trash, 'bed')
      BEDTOBIGBED_TRASH(JOIN_VEREDICT_TRASH.out.output, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.retentions
        .collect()
        .set { ch_retentions }
      JOIN_VEREDICT_RETENTIONS(ch_retentions, 'bed')
      BEDTOBIGBED_RETENTIONS(JOIN_VEREDICT_RETENTIONS.out.output, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.truncations
        .collect()
        .set { ch_truncations }
      JOIN_VEREDICT_TRUNCATIONS(ch_truncations, 'bed')
      BEDTOBIGBED_TRUNCATIONS(JOIN_VEREDICT_TRUNCATIONS.out.output, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.intraprimming
        .collect()
        .set { ch_intraprimming }
      JOIN_VEREDICT_INTRAPRIMMING(ch_intraprimming, 'bed')
      BEDTOBIGBED_INTRAPRIMMING(JOIN_VEREDICT_INTRAPRIMMING.out.output, chrom_sizes, autosql)

      ch_bbs = Channel.empty()
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_PASSES.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_TRASH.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_RETENTIONS.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_TRUNCATIONS.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_INTRAPRIMMING.out.bigbed)

    emit:
      bigbeds               = ch_bbs
      versions              = ch_versions
}
