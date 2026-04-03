process ISOTOOLS_INTRON_RETENTION {
    tag "$meta.id:$meta.chr"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isotools:latest' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta1), path(introns)

    output:
    tuple val(meta), path("*.tsv")       , optional: true, emit: descriptor
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}_${meta.chr}"
    """
    iso-intron \\
        $args \\
        --introns $introns \\
        --query $bed \\
        --threads ${task.cpus} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-intron: \$( iso-intron --version | sed 's/iso-intron //g' )
    END_VERSIONS
    """

    stub:
    """
    touch *.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-intron: \$( iso-intron --version | sed 's/iso-intron //g' )
    END_VERSIONS
    """
}
