process PUBLISH {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path(files)

    script:
    """
    """
}
