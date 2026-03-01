process PBCCS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbccs:6.4.0--h9ee0642_0' :
        'biocontainers/pbccs:6.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(pbi)
    val chunk_on

    output:
    tuple val(meta), path("*.chunk*.bam")     , emit: bam
    tuple val(meta), path("*.chunk*.bam.pbi") , emit: pbi
    tuple val(meta), path("*.report.txt" )    , emit: report_txt
    tuple val(meta), path("*.report.json" )   , emit: report_json
    tuple val(meta), path("*.metrics.json.gz"), emit: metrics
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ccs \\
        $bam \\
        ${prefix}.chunk${meta.chunk}.bam \\
        --report-file ${prefix}.chunk${meta.chunk}.report.txt \\
        --report-json ${prefix}.chunk${meta.chunk}.report.json \\
        --metrics-json ${prefix}.chunk${meta.chunk}.metrics.json.gz \\
        --chunk $meta.chunk/$chunk_on \\
        -j $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbccs: \$(echo \$(ccs --version 2>&1) | grep 'ccs' | sed 's/^.*ccs //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch *.chunk*.bam
    touch *.chunk*.bam.pbi
    touch *.report.txt
    touch *.report.json
    gzip *.metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbccs: \$(echo \$(ccs --version 2>&1) | grep 'ccs' | sed 's/^.*ccs //; s/ .*\$//')
    END_VERSIONS
    """
}
