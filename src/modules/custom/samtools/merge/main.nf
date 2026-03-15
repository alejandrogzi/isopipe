process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23--h96c455f_0' :
        'biocontainers/samtools:1.23--h96c455f_0' }"

    input:
    tuple val(meta), path(sams, stageAs: "sam/*")

    output:
    tuple val(meta), path("${meta.id}.bam")              , optional: true, emit: bam
    tuple val(meta), path("${meta.id}.bam.bai")          , optional: true, emit: bai
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def keep_temp = task.ext.keep_temp ?: false
    """
    for sam in ${sams}; do
        base=\$(basename "\$sam" .sam)
        bam="\${base}.sorted.bam"

        samtools view \\
            -@ ${task.cpus} \\
            -b $args \\
            "\$sam" | \\
        samtools sort \\
            -@ ${task.cpus} \\
            -o "\$bam"

        samtools index \\
            -@ ${task.cpus} \\
            "\$bam"

        if [ "${keep_temp}" == "false" ]; then
            if [ -L "\$sam" ]; then
                realpath=\$(readlink -f "\$sam")
                rm -f "\$realpath"
            else
                rm -f "\$sam"
            fi
        fi
    done

    # INFO: merge all sorted BAMs into a single BAM
    samtools merge \\
        -@ ${task.cpus} \\
        -f $args2 \\
        ${meta.id}.bam \\
        *.sorted.bam

    # INFO: index the merged BAM
    samtools index \\
        -@ ${task.cpus} \\
        ${meta.id}.bam

    # INFO: clean up intermediate sorted BAMs and their indices
    rm -f *.sorted.bam *.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
