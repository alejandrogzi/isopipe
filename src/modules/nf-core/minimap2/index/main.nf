process MINIMAP2_INDEX {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.31--h118bc1c_0' :
        'biocontainers/minimap2:2.31--h118bc1c_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    tuple val("${task.process}"), val("minimap2"), eval("minimap2 --version"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
