process GAWK_JOIN_BEDGRAPH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0':
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(bedgraphs, stageAs: "bedgraph/*")

    output:
    tuple val(meta), path("*.bg"), emit: bedgraph
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gawk 1 bedgraph/* > ${prefix}.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed 's/^.*awk version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed 's/^.*awk version //; s/ .*//')
    END_VERSIONS
    """
}               
