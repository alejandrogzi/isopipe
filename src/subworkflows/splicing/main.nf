/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINISPLICE_DOWNLOAD } from '../../modules/custom/minisplice/download/main.nf'
include { MINISPLICE_PREDICT } from '../../modules/custom/minisplice/predict/main.nf'
include { SPLICEAI_DERIVE } from '../../modules/custom/spliceai/derive/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow SPLICING {
    take:
      genome
      annotation
      bigwigs
      algorithm
      ch_versions

    main:
      if (algorithm == "spliceai") {
          if (bigwigs) {
            Channel.value([
                    [ id: "spliceai" ],
                    file(bigwigs, checkIfExists: true)
            ]).set { ch_spliceai_bigwigs }

            SPLICEAI_DERIVE(
                genome.map { genome -> [ [id:genome.baseName], genome ] },
                annotation.map { annotation -> [ [id:annotation.baseName], annotation ] },
                ch_spliceai_bigwigs
            )
            ch_scores = SPLICEAI_DERIVE.out.scores
            ch_versions = ch_versions.mix(SPLICEAI_DERIVE.out.versions)
          } else {
              // SPLICEAI_RUN()

              // ch_scores = SPLICEAI_RUN.out.scores
              // ch_versions = ch_versions.mix(SPLICEAI_RUN.out.versions)
          }
      } else if (algorithm == "minisplice") {
          MINISPLICE_DOWNLOAD()
          MINISPLICE_PREDICT(
              genome.map { genome -> [ [id:genome.baseName], genome ] },
              MINISPLICE_DOWNLOAD.out.model,
              MINISPLICE_DOWNLOAD.out.calibration
          )

          ch_scores = MINISPLICE_PREDICT.out.scores
          ch_versions = ch_versions.mix(MINISPLICE_DOWNLOAD.out.versions)
          ch_versions = ch_versions.mix(MINISPLICE_PREDICT.out.versions)
      }

    emit:
      scores = ch_scores
      versions = ch_versions
}
