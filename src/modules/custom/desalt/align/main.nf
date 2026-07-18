/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DESALT_ALIGN — Align reads to a genome using deSALT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DESALT_ALIGN {
    tag "${meta.id}.${meta.chunk}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/desalt:latest' }"


    input:
    tuple val(meta), path(reads)
    tuple val(meta1), path(index)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.sam")        , emit: sam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def annotation = gtf ? "--gtf ${gtf}" : ''
    def is_singleton = meta.singleton ? "singleton" : "hq"
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}.${is_singleton}"

    """
    deSALT \\
        aln \\
        ${args} \\
        ${annotation} \\
        --thread ${task.cpus} \\
        --output ${prefix}.sam \\
        ${index} \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deSALT: \$(deSALT | grep 'Version' | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def is_singleton = meta.singleton ? "singleton" : "hq"
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}.${is_singleton}"
    """
    touch ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deSALT: \$(deSALT | grep 'Version' | awk '{print \$2}')
    END_VERSIONS
    """
}
