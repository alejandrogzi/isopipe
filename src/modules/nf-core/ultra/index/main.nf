process ULTRA_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0':
        'quay.io/biocontainers/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("index"), emit: index
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def disable_infer = gtf ? "--disable_infer" : ""
    """
    mkdir index
    uLTRA \\
        index \\
        ${disable_infer} \\
        ${args} \\
        ${fasta} \\
        ${gtf} \\
        index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
        gffutils: \$(python -c "import gffutils; print(gffutils.__version__)")
        sqlite: \$(python -c "import sqlite3; print(sqlite3.sqlite_version)")
    END_VERSIONS
    """

    stub:
    """
    mkdir index
    touch index/database.db
    touch index/all_splice_pairs_annotations.pickle
    touch index/all_splice_sites_annotations.pickle
    touch index/chr_to_id.pickle
    touch index/exon_choordinates_to_id.pickle
    touch index/flank_choordinates.pickle
    touch index/gene_to_small_segments.pickle
    touch index/id_to_chr.pickle
    touch index/max_intron_chr.pickle
    touch index/parts_to_segments.pickle
    touch index/ref_exon_sequences.pickle
    touch index/ref_flank_sequences.pickle
    touch index/ref_part_sequences.pickle
    touch index/ref_segment_sequences.pickle
    touch index/refs_id_lengths.pickle
    touch index/refs_lengths.pickle
    touch index/segment_id_to_choordinates.pickle
    touch index/segment_to_gene.pickle
    touch index/segment_to_ref.pickle
    touch index/splices_to_transcripts.pickle
    touch index/transcripts_to_splices.pickle

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
        gffutils: \$(python -c "import gffutils; print(gffutils.__version__)")
        sqlite: \$(python -c "import sqlite3; print(sqlite3.sqlite_version)")
    END_VERSIONS
    """
}
