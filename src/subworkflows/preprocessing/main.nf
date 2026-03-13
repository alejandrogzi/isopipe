/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOSEQ } from '../isoseq/main.nf'
include { GENOME } from '../genome/main.nf'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { SPLICING as GENOMIC_SPLICE_SCORES } from '../splicing/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {
    take:
      entrypoint             // isoseq or fastq
      global_input_dir       // path
      genome                 // path
      toga2                  // path
      protein_database       // path
      minimap2_index         // path
      use_splice_scores      // bool
      splice_score_algorithm // string
      ch_versions            // [ meta, versions.yml ]

    main:
      ch_genome = GENOME(genome)
      ch_reference_transcripts = Channel.value(file(toga2, checkIfExists: true))
      ch_database = Channel.value(file(protein_database, checkIfExists: true))

      ch_minimap2_index = Channel.empty()
      if (minimap2_index) {
          ch_minimap2_index = Channel.value(file(minimap2_index, checkIfExists: true))
              .map { path -> [ [id:path.name], path ] }
      } else {
          ch_minimap2_index = MINIMAP2_INDEX(
              ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
          )
      }

      if (use_splice_scores) {
          if (splice_score_algorithm == "spliceai") {
              def spliceai_use_pre_computed_scores = params.spliceai_bigwigs_dir ? true : false
              // if (spliceai_use_pre_computed_scores) {
              //    Channel
              //      .fromPath("${params.spliceai_bigwigs_dir}/*.bw", checkIfExists: true)
              //    .collect()
              //     .map { bws -> [ [ id:"spliceai" , bws ] }
              //      .set { ch_spliceai_bigwigs }
              // } 
              // else {
                // GENOMIC_SPLICE_SCORES(
                //    ch_genome.genome,  
                // )
              // }
          } else if (splice_score_algorithm == "minisplice") {
              GENOMIC_SPLICE_SCORES(
                  ch_genome.genome,
                  splice_score_algorithm,
                  ch_versions
              )             
          }
      }

      // INFO: isoseq entrypoint
      if (entrypoint == "isoseq") {
          ISOSEQ()
          ch_reads = ch_reads.mix(ISOSEQ.out.reads)
          ch_versions = ch_versions.mix(ISOSEQ.out.versions)
      } else if (entrypoint == "map") {
          Channel
              .fromPath("${global_input_dir}/*.fast*", checkIfExists: true)
              .map { fastx ->
                  def singleton = fastx.baseName.contains("singleton")
                  return [
                      [
                          id:         fastx.baseName,
                          single_end: true,
                          singleton:  singleton
                      ],
                      fastx,
                  ]
              }
              .set { ch_reads }
      }
      
      ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
      ch_versions = ch_versions.mix(ch_genome.versions)

    emit:
      genome         = ch_genome.genome
      database       = ch_database
      reads          = ch_reads
      minimap2_index = ch_minimap2_index.index
      reference_transcripts = ch_reference_transcripts
      versions = ch_versions
}
