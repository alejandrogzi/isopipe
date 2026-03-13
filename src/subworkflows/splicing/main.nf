/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINISPLICE_DOWNLOAD } from '../../modules/custom/minisplice/download/main.nf'
include { MINISPLICE_PREDICT } from '../../modules/custom/minisplice/predict/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow SPLICING {
    take:
      genome
      algorithm
      ch_versions

    main:
      if (algorithm == "spliceai") {
          Channel
              .fromPath("${params.spliceai_bigwigs_dir}/*.bw", checkIfExists: true)
              .collect()
              .map { bws -> [ [ id:"spliceai" ] , bws ] }
              .set { ch_spliceai_bigwigs }
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
