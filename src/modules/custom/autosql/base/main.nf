process AUTOSQL_BASE {
    tag "autosql base"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/hillerlab/git:latest' }"

    output:
    path "*.as",          emit: autosql
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat <<-EOF > base.as
    table ark
    "Schema for base ark-produced long-read transcriptomic reads [fusions + nmd]"
    (
        string  chrom;              "Reference sequence chromosome or scaffold"
        uint    chromStart;         "Start position of feature on chromosome"
        uint    chromEnd;           "End position of feature on chromosome"
        string  name;               "Feature ID"
        uint    score;              "Score (0-1000)"
        char[1] strand;             "+ or - for strand"
        uint    thickStart;         "Coding region start"
        uint    thickEnd;           "Coding region end"
        uint    itemRgb;           "RGB color value"
        int     blockCount;         "Number of exons"
        int[blockCount] blockSizes; "Exon sizes"
        int[blockCount] chromStarts;"Exon starts relative to chromStart"
    )
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -n1 | sed 's/Bash (GNU bash) //g; s/  .*//')
    END_VERSIONS
    """

    stub:
    """
    touch base.as

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -n1 | sed 's/Bash (GNU bash) //g; s/  .*//')
    END_VERSIONS
    """
}
