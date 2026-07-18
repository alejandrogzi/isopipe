process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def chunk = meta.chunk ? ".${meta.chunk}" : ""
    def singleton = meta.singleton ? ".singleton" : ""
    def prefix = task.ext.prefix ?: "${meta.id}${chunk}${singleton}"
    """
    samtools \\
        sort \\
        -@ ${task.cpus} \\
        $input \\
        > ${prefix}.sorted.bam

    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${args} \\
        ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def chunk = meta.chunk ? ".${meta.chunk}" : ""
    def singleton = meta.singleton ? ".singleton" : ""
    def prefix = task.ext.prefix ?: "${meta.id}${chunk}${singleton}"
    """
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
