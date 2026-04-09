/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RSYNC_SSH } from '../../modules/custom/ssh/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOAD_TRACK {
    take:
      bigbed                 // channel: [ [ meta ], [ bigbed ] ]
      user                   // string
      server                 // string
      target_dir             // string
      web                    // string
      species                // string
      ch_versions            // [ meta, versions.yml ]

    main:
      RSYNC_SSH(
        bigbed,
        user,
        server,
        target_dir,
        web,
        species,
      )

      ch_versions = ch_versions.mix(RSYNC_SSH.out.versions)

    emit:
      versions              = ch_versions
}
