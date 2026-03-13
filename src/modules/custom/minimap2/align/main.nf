process MINIMAP2_ALIGN {
    tag "$meta.id chunk $meta.chunk"
    label 'process_high'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.bam")                       , optional: true, emit: bam
    tuple val(meta), path("*.bai")                       , optional: true, emit: bai
    tuple val(meta), path("*.sam")                       , optional: true, emit: sam
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def keep_temp = task.ext.keep_temp ?: false
    
    def singleton = meta.singleton ? ".singleton" : ""
    def sam = "${meta.id}.${meta.chunk}${singleton}.sam"
    def bam = "${meta.id}.${meta.chunk}${singleton}.bam"
    def bam_index = bam + '.bai'
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        ${reference} \\
        ${reads} \\
        -o $sam

    samtools view -@ ${task.cpus} -b $sam | samtools sort -@ ${task.cpus} -o $bam
    samtools index -@ {task.cpus} $bam

    if [ $keep_temp == false ]; then
        rm $sam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}${singleton}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai
    touch ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
