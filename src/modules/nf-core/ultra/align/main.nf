process ULTRA_ALIGN {
    tag "$meta.id chunk $meta.chunk"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0':
        'quay.io/biocontainers/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(genome)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def singleton = meta.singleton ? ".singleton" : ""
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}${singleton}"
    """
    uLTRA \\
        align \\
        --t $task.cpus \\
        --prefix $prefix \\
        --isoseq \\
        --index $index \\
        $args \\
        $genome \\
        $reads \\
        ./ 

    samtools \\
        sort \\
        --threads $task.cpus \\
        -o ${prefix}.bam \\
        -O BAM \\
        $args2 \\
        ${prefix}.sam

    samtools index ${prefix}.bam

    rm ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}${singleton}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
