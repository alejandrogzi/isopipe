/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the GNU General Public License, Version 3.0.
*/

process FASTA_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    tuple val(meta), path(custom), path(raw)

    output:
    tuple val(meta), path("${prefix}.fa"), emit: fasta
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: 'merged'
    def custom_is_gz = custom.toString().endsWith('.gz')
    def raw_is_gz = raw.toString().endsWith('.gz')
    """
    ${custom_is_gz ? "gzip -c -d ${custom}" : "cat ${custom}"} > custom.fa
    ${raw_is_gz ? "gzip -c -d ${raw}" : "cat ${raw}"} > raw.fa
    cat raw.fa custom.fa > ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //' )
        gzip: \$(gzip --version | head -n1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'merged'
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //' )
        gzip: \$(gzip --version | head -n1 | sed 's/^.* //')
    END_VERSIONS
    """
}
