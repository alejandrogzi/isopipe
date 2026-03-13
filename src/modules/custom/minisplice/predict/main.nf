process MINISPLICE_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minisplice:0.4--h577a1d6_0' :
        'quay.io/biocontainers/minisplice:0.4--h577a1d6_0' }"


    input:
    tuple val(meta), path(genome)
    tuple val(meta1), path(model)
    tuple val(meta2), path(calibration)

    output:
    tuple val(meta), path("*_splice_scores.tsv"), emit: scores
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minisplice \\
        predict \\
        $args \\
        -t ${task.cpus} \\
        -c ${calibration} ${model} ${genome} \\
        > ${prefix}_splice_scores.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minisplice: \$( minisplice --version | sed 's/minisplice //g' )
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}_splice_scores.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minisplice: \$( minisplice --version | sed 's/minisplice //g' )
    END_VERSIONS
    """
}
