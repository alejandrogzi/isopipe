/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOTOOLS_INTRON_RETENTION } from '../../modules/custom/isotools/intron/main.nf'
include { ISOTOOLS_PAS_CALLER } from '../../modules/custom/isotools/pas/main.nf'
include { ISOTOOLS_TRUNCATION_DETECTOR } from '../../modules/custom/isotools/utr/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POLISH {
    take:
      reads                  // channel: [ val(meta), [ reads ] ]
      annotation             // Channel.value(path)
      introns                // channel: [ val(meta), [ introns ] ]
      forward_peaks          // channel: [ val(meta), [ forward_peaks ] ]
      reverse_peaks          // channel: [ val(meta), [ reverse_peaks ] ]
      ch_versions            // [ meta, versions.yml ]

    main:
      ch_reference_transcripts = annotation.map { annotation -> [ [id:annotation.baseName], annotation ] }

      ISOTOOLS_INTRON_RETENTION(
        reads,
        introns
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

    emit:
      versions              = ch_versions
}
