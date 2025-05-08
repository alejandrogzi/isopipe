#!/usr/bin/env nextflow
params.chunkSize = 10000
params.workDir = "TEMP_DIAMOND_NF"

db_name = file(params.db).name

workflow {
    Channel
        .fromPath(params.query)
        .splitFasta(by: params.chunkSize, file:true)
        .set { ch_fasta }

    ch_hits = blast(ch_fasta, params.db)
    ch_hits.collectFile(name: params.out)
}


process blast {
    module 'diamond'
    errorStrategy = { task.attempt <= 3 ? 'retry' : 'ignore' } 
    memory = "35GB"
    

    input:
    path 'query.fa'
    path db

    output:
    path 'blast_result'

    """
    diamond blastp --query query.fa --db $db_name > blast_result
    """
}
