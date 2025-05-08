#!/usr/bin/env nextflow
params.chunkSize = 200
params.workDir = "TEMP_TRANSLATIONAI_NF"

workflow {
    Channel
        .fromPath(params.query)
        .splitFasta(by: params.chunkSize, file:true)
        .set { ch_fasta }

    ch_hits = predict(ch_fasta)
    ch_hits.collectFile(name: params.out)
}


process predict {
    module 'translationai/920ff06'
    errorStrategy = { task.attempt <= 3 ? 'retry' : 'ignore' } 
    memory = "15GB"
    

    input:
    path 'query.fa'

    output:
    path 'predictions'

    """
    translationai  -I query.fa -t 0.01,0.01 -o predictions
    """
}
