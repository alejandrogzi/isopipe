process DETACH_DUPLICATES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43bff7_4':
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.deduplicated.bed"),   optional: true, emit: pass
    tuple val(meta1), path("*.duplicates.bed"),    optional: true, emit: duplicates
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    meta1 = meta.clone()
    meta1.id  = "${meta.id}.duplicates"
    """
    grep -v '#DU' ${bed} > ${prefix}.deduplicated.bed || [[ \$? == 1 ]]
    grep    '#DU' ${bed} > ${prefix}.duplicates.bed   || [[ \$? == 1 ]]

    if [[ ! -s ${prefix}.deduplicated.bed ]]; then
        rm ${prefix}.deduplicated.bed
    fi
    if [[ ! -s ${prefix}.duplicates.bed ]]; then
        rm ${prefix}.duplicates.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | sed 's/grep (GNU grep) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.deduplicated.bed
    touch ${prefix}.duplicates.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | sed 's/grep (GNU grep) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
