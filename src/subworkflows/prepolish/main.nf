/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XLOC_INTRON as XLOCI_EXTRACT_INTRONS } from '../../modules/custom/xloci/intron/main.nf'
include { INTRONIC as IIC_PREDICT_SPLICEOSOME } from '../../modules/custom/intronic/main.nf'
include { ISOTOOLS_CLASSIFY_INTRON } from '../../modules/custom/isotools/classify/intron/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPOLISH {
    take:
      reads                  // channel: [ val(meta), [ reads ] ]
      genome                 // Channel.value(path)
      repeats                // path
      annotation             // Channel.value(path)
      bigwigs                // path
      ch_versions            // [ meta, versions.yml ]

    main:
      ch_genome = genome.map { genome -> [ [id:genome.baseName], genome ] }
      ch_reference_transcripts = annotation.map { annotation -> [ [id:annotation.baseName], annotation ] }

      XLOCI_EXTRACT_INTRONS(ch_genome, reads)
      IIC_PREDICT_SPLICEOSOME(XLOCI_EXTRACT_INTRONS.out.tsv)
    
      if (bigwigs) {
          Channel.value([
              [ id: "spliceai" ],
              file(bigwigs, checkIfExists: true)
          ]).set { ch_spliceai_bigwigs }
      } else {
          ch_spliceai_bigwigs = Channel.empty()
      }

      if (repeats) {
          Channel.value([
              [ id: "repeats" ],
              file(repeats, checkIfExists: true)
          ]).set { ch_repeats }
      } else {
          ch_repeats = Channel.empty()
      }

      ISOTOOLS_CLASSIFY_INTRON(
        reads,
        ch_genome,
        ch_reference_transcripts,
        ch_repeats,
        ch_spliceai_bigwigs,
        IIC_PREDICT_SPLICEOSOME.out.iic
      )

     // XISO_CHUNK_FOR_APARENT()
     // APARENT()

      
      ch_versions = ch_versions.mix(ch_genome.versions)

    emit:
      genome                = ch_genome.genome
      database              = ch_database
      reads                 = ch_reads
      minimap2_index        = ch_minimap2_index
      reference_transcripts = ch_reference_transcripts
      splice_scores         = ch_splice_scores
      versions              = ch_versions
}
