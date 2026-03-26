/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XLOCI_INTRON as XLOCI_EXTRACT_INTRONS } from '../../modules/custom/xloci/intron/main.nf'
include { INTRONIC as IIC_PREDICT_SPLICEOSOME } from '../../modules/custom/intronic/main.nf'
include { ISOTOOLS_CLASSIFY_INTRON } from '../../modules/custom/isotools/classify/intron/main.nf'
include { APARENT_CHUNK as XISO_APARENT_CHUNK } from '../../modules/custom/aparent/chunk/main.nf'
include { WGET as WGET_APARENT_WEIGHTS } from '../../modules/nf-core/wget/main.nf'

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
      bigwigs                // channel: [ val(meta), [ bigwigs ] ]
      aparent_weights        // path
      ch_versions            // [ meta, versions.yml ]

    main:
      ch_genome = genome.map { genome -> [ [id:genome.baseName], genome ] }
      ch_reference_transcripts = annotation.map { annotation -> [ [id:annotation.baseName], annotation ] }

      ch_aparent_weights = WGET_APARENT_WEIGHTS(
          Channel.value(
            aparent_weights
          ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
      )

      XLOCI_EXTRACT_INTRONS(ch_genome, reads)
      IIC_PREDICT_SPLICEOSOME(XLOCI_EXTRACT_INTRONS.out.tsv)
    
      if (repeats) {
          Channel.value([
              [ id: "repeats" ],
              file(repeats, checkIfExists: true)
          ]).set { ch_repeats }
      } else {
          ch_repeats = Channel.value([[:], []])
      }

      ISOTOOLS_CLASSIFY_INTRON(
        reads,
        ch_genome,
        ch_reference_transcripts,
        ch_repeats,
        bigwigs,
        IIC_PREDICT_SPLICEOSOME.out.iic
      )

     XISO_APARENT_CHUNK(
        reads,
        ch_genome
     )

     XISO_APARENT_CHUNK.out.chunks
        .flatMap { 
            meta, bed ->
            def beds = bed instanceof List ? bed : [bed]
            beds.withIndex().collect { it, idx ->
                [ meta + [ chunk: idx ], it ]
            }
        }
        .set { ch_aparent_chunks }

      APARENT_PREDICT(
        ch_aparent_chunks, 
        ch_aparent_weights
      )
        .map { meta, bed -> bed }
        .collect()
        .map { beds -> [ [id:'aparent'], beds ] }
        .set { ch_joined_aparent_beds }

      ch_versions = ch_versions.mix(XLOCI_EXTRACT_INTRONS.out.versions)
      ch_versions = ch_versions.mix(IIC_PREDICT_SPLICEOSOME.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_CLASSIFY_INTRON.out.versions)
      ch_versions = ch_versions.mix(XISO_APARENT_CHUNK.out.versions)
      ch_versions = ch_versions.mix(APARENT_PREDICT.out.versions)

    emit:
      introns               = ISOTOOLS_CLASSIFY_INTRON.out.tsv
      versions              = ch_versions
}
