process COLLAPSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isox-rs:latest' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("collapsed/*.tsv"), emit: collapsed
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    collapse \\
        $args \\
        -t ${task.cpus} \\
        --bed ${bed} \\
        --extend \\
        --collapse-mode cds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse: \$( collapse --version | sed 's/collapse //g' )
    END_VERSIONS
    """

    stub:
    """
    touch *.derived.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse: \$( collapse --version | sed 's/collapse //g' )
    END_VERSIONS
    """
}
