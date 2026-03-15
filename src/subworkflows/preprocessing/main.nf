/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOSEQ } from '../isoseq/main.nf'
include { GENOME } from '../genome/main.nf'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { SPLICING as MINISPLICE_GENOMIC_SPLICE_SCORES } from '../splicing/main.nf'
include { SPLICING as SPLICEAI_GENOMIC_SPLICE_SCORES } from '../splicing/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {
    take:
      entrypoint             // isoseq or fastq
      global_input_dir       // path
      global_primers         // path
      genome                 // path
      annotation             // path
      ccs_chunk              // int
      isoseq_cluster2_mode   // string
      protein_database       // path
      minimap2_index         // path
      use_splice_scores      // bool
      splice_score_algorithm // string
      bigwigs                // path
      minisplice             // path
      ch_versions            // [ meta, versions.yml ]

    main:
      ch_genome = GENOME(genome)
      ch_reference_transcripts = Channel.value(file(annotation, checkIfExists: true))
      ch_database = Channel.value(file(protein_database, checkIfExists: true))

      ch_reads = Channel.empty()

      ch_minimap2_index = Channel.empty()
      if (minimap2_index) {
          ch_minimap2_index = Channel.value(file(minimap2_index, checkIfExists: true))
              .map { path -> [ [id:path.name], path ] }
      } else {
          MINIMAP2_INDEX(
              ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
          )
          ch_minimap2_index = MINIMAP2_INDEX.out.index
          ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
      }

      ch_splice_scores = Channel.empty()
      if (use_splice_scores) {
          if (splice_score_algorithm == "spliceai") {
              SPLICEAI_GENOMIC_SPLICE_SCORES(
                  ch_genome.genome,
                  ch_reference_transcripts,
                  bigwigs,
                  splice_score_algorithm,
                  minisplice,
                  ch_versions
              )
              ch_splice_scores = SPLICEAI_GENOMIC_SPLICE_SCORES.out.scores
          } else if (splice_score_algorithm == "minisplice") {
              MINISPLICE_GENOMIC_SPLICE_SCORES(
                  ch_genome.genome,
                  ch_reference_transcripts,
                  bigwigs,
                  splice_score_algorithm,
                  minisplice,
                  ch_versions
              )             
              ch_splice_scores = MINISPLICE_GENOMIC_SPLICE_SCORES.out.scores
          }
      } else {
          // INFO: keep MINIMAP2_ALIGN schedulable when splice scores are disabled
          ch_splice_scores = Channel.value([[:], []])
      }

      // INFO: isoseq entrypoint
      if (entrypoint == "isoseq") {
          ISOSEQ(
              global_input_dir,
              global_primers,
              ccs_chunk,
              isoseq_cluster2_mode
          )

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
