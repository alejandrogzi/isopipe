process APARENT_PREDICT {
    tag "$meta.id $meta.chunk"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/alejandrogzi/isox-py:latest' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta1), path(weights)

    output:
    tuple val(meta), path("*.aparent.bed")  , optional: true, emit: bed
    tuple val(meta), path("*.aparent.bg")   , optional: true, emit: bg
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}.${meta.chunk}"
    """
    aparent \\
        --bed $bed \\
        --outdir aparent \\
        --prefix $prefix \\
        --model $aparent_weights

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aparent: \$(aparent --version 2>&1 | sed 's/aparent //')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}.${meta.chunk}"
    """
    touch aparent/*.aparent.bed 
    touch aparent/*.aparent.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aparent: \$(aparent --version 2>&1 | sed 's/aparent //')
    END_VERSIONS
    """
} 
