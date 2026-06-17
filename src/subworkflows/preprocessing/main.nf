/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOSEQ } from '../isoseq/main.nf'
include { GENOME } from '../genome/main.nf'
include { GXF2BED } from '../../modules/custom/gxf2bed/main.nf'
include { BED2GTF } from '../../modules/custom/bed2gtf/main.nf'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { ULTRA_INDEX } from '../../modules/nf-core/ultra/index/main.nf'
include { GENEPRED_LINT } from '../../modules/custom/genepred/lint/main.nf'
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
      spliceai               // path
      compression            // bool
      global_prefix          // string
      aligner                // string [ minimap2, ultra ]
      ultra_use_annotation   // bool
      ultra_index            // path
      ch_versions            // [ meta, versions.yml ]

    main:
      ch_genome = GENOME(genome)
      ch_database = Channel.value(file(protein_database, checkIfExists: true))

      ch_reference_transcripts = Channel.value(file(annotation, checkIfExists: true))
      ch_reference_transcripts_gtf = Channel.empty()
      GENEPRED_LINT(ch_reference_transcripts.map { file -> [ [ id:file.baseName ], file ] })

      if (annotation.endsWith(".gtf.gz")
        | annotation.endsWith(".gtf") 
        | annotation.endsWith(".gff.gz") 
        | annotation.endsWith(".gff")) {
          ch_reference_transcripts_gtf = ch_reference_transcripts
            .map { gtf -> [ [id:gtf.baseName], gtf ] }

          GXF2BED(ch_reference_transcripts_gtf)
          ch_reference_transcripts = GXF2BED.out.bed
      } else {
          BED2GTF(
            ch_reference_transcripts,
            Channel.of([[:], []])
          )

          ch_reference_transcripts_gtf = BED2GTF.out.gtf
          ch_reference_transcripts_gtf = ch_reference_transcripts.map { bed -> [ [id:bed.baseName], bed ] }
      }

      ch_reads = Channel.empty()
      ch_genome_index = Channel.empty()
      if (aligner == "minimap2") {
        if (minimap2_index) {
            ch_genome_index = Channel.value(file(minimap2_index, checkIfExists: true))
                .map { path -> [ [id:path.name], path ] }
        } else {
            MINIMAP2_INDEX(
                ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
            )
            ch_genome_index = MINIMAP2_INDEX.out.index
            ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        }
      } else if (aligner == "ultra") {
          if (ultra_index) {
              ch_genome_index = Channel.value(file(ultra_index, checkIfExists: true))
                  .map { path -> [ [id:path.name], path ] }
          } else {
              if (ultra_use_annotation) {
                ULTRA_INDEX(
                    ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] },
                    ch_reference_transcripts_gtf
                )
              } else {
                ULTRA_INDEX(
                    ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] },
                    Channel.of([[:], []])
                )
              }

              ch_genome_index = ULTRA_INDEX.out.index
              ch_versions = ch_versions.mix(ULTRA_INDEX.out.versions)
          }
      }

      ch_splice_scores = Channel.empty()
      ch_spliceai_bigwigs = Channel.value([[:], []])
      if (use_splice_scores) {
          if (splice_score_algorithm == "spliceai") {
              SPLICEAI_GENOMIC_SPLICE_SCORES(
                  ch_genome.genome,
                  ch_reference_transcripts,
                  bigwigs,
                  splice_score_algorithm,
                  minisplice,
                  spliceai,
                  ch_genome.chrom_sizes,
                  compression,
                  ch_versions
              )
              ch_splice_scores = SPLICEAI_GENOMIC_SPLICE_SCORES.out.scores
              ch_spliceai_bigwigs = SPLICEAI_GENOMIC_SPLICE_SCORES.out.bigwigs
          } else if (splice_score_algorithm == "minisplice") {
              MINISPLICE_GENOMIC_SPLICE_SCORES(
                  ch_genome.genome,
                  ch_reference_transcripts,
                  bigwigs,
                  splice_score_algorithm,
                  minisplice,
                  spliceai,
                  ch_genome.chrom_sizes,
                  compression,
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
              isoseq_cluster2_mode,
              global_prefix
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
      genome                    = ch_genome.genome
      database                  = ch_database
      reads                     = ch_reads
      genome_index              = ch_genome_index
      reference_transcripts     = ch_reference_transcripts
      reference_transcripts_gtf = ch_reference_transcripts_gtf
      splice_scores             = ch_splice_scores
      bigwigs                   = ch_spliceai_bigwigs
      chrom_sizes               = ch_genome.chrom_sizes
      versions                  = ch_versions
}
