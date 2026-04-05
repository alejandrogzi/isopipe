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
include { APARENT_PREDICT } from '../../modules/custom/aparent/predict/main.nf'
include { BEDGRAPHTOBIGWIG as BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD } from '../../modules/custom/bigtools/bedgraphtobigwig/main.nf'
include { BEDGRAPHTOBIGWIG as BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE } from '../../modules/custom/bigtools/bedgraphtobigwig/main.nf'
include { GAWK_JOIN as GAWK_JOIN_BEDGRAPH_FORWARD } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as GAWK_JOIN_BEDGRAPH_REVERSE } from '../../modules/custom/gawk/join/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPOLISH {
    take:
      reads                  // channel: [ val(meta), [ reads ] ]
      genome                 // Channel.value(path)
      chrom_sizes            // Channel.value(path)
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

      reads
        .map { meta, read -> tuple(meta.id, meta, read) }
        .join(
          IIC_PREDICT_SPLICEOSOME.out.iic
            .map { meta, iic -> tuple(meta.id, meta, iic) }
        )
        .map { id, read_meta, read, iic_meta, iic ->
          tuple(read_meta, read, iic)
        }
        .set { ch_classify_inputs }

      ISOTOOLS_CLASSIFY_INTRON(
        ch_classify_inputs,
        ch_genome,
        ch_reference_transcripts,
        ch_repeats,
        bigwigs,
      )

     XISO_APARENT_CHUNK(
        reads,
        ch_genome
     )

     XISO_APARENT_CHUNK.out.chunks
        .flatMap { 
            meta, chunk_tsv ->
            def chunk_tsvs = chunk_tsv instanceof List ? chunk_tsv : [chunk_tsv]
            chunk_tsvs.withIndex().collect { it, idx ->
                [ meta + [ chunk: idx ], it ]
            }
        }
        .set { ch_aparent_chunks }

      APARENT_PREDICT(
        ch_aparent_chunks, 
        ch_aparent_weights.outfile
      )

      APARENT_PREDICT.out.bg_forward
        .map { meta, bg -> bg }
        .collect()
        .map { bgs -> [ [id:'aparent.forward', strand:'forward'], bgs ] }
        .set { ch_joined_aparent_bgs_forward }

      GAWK_JOIN_BEDGRAPH_FORWARD(ch_joined_aparent_bgs_forward, 'bg')
      BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD(GAWK_JOIN_BEDGRAPH_FORWARD.out.output, chrom_sizes)

      APARENT_PREDICT.out.bg_reverse
        .map { meta, bg -> bg  }
        .collect()
        .map { bgs -> [ [id:'aparent.reverse', strand:'reverse'], bgs ] }
        .set { ch_joined_aparent_bgs_reverse }

      GAWK_JOIN_BEDGRAPH_REVERSE(ch_joined_aparent_bgs_reverse, 'bg')
      BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE(GAWK_JOIN_BEDGRAPH_REVERSE.out.output, chrom_sizes)

      ch_versions = ch_versions.mix(XLOCI_EXTRACT_INTRONS.out.versions)
      ch_versions = ch_versions.mix(IIC_PREDICT_SPLICEOSOME.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_CLASSIFY_INTRON.out.versions)
      ch_versions = ch_versions.mix(XISO_APARENT_CHUNK.out.versions)
      ch_versions = ch_versions.mix(APARENT_PREDICT.out.versions)
      ch_versions = ch_versions.mix(GAWK_JOIN_BEDGRAPH_FORWARD.out.versions)
      ch_versions = ch_versions.mix(GAWK_JOIN_BEDGRAPH_REVERSE.out.versions)
      ch_versions = ch_versions.mix(BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD.out.versions)
      ch_versions = ch_versions.mix(BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE.out.versions)

    emit:
      introns               = ISOTOOLS_CLASSIFY_INTRON.out.tsv
      aparent_plus          = BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD.out.bigwig
      aparent_minus         = BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE.out.bigwig
      versions              = ch_versions
}
