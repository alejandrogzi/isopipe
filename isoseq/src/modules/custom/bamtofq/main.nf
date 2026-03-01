process BAM_TO_FQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def singleton_bam = prefix + ".singletons.bam"
    def singleton_fq  = prefix + ".singletons.fastq.gz"
    def hq_bam        = prefix + ".hq.bam"
    def hq_fq         = prefix + ".hq.fastq.gz"
    """
    samtools view -h -@ 4 ${bam} \\
        | tee >(awk '{{if(\$1 ~ /^@/){{print; next}} for(i=12;i<=NF;i++){{if(\$i ~ /^is:i:1\$/){{print; break}}}}}}' \\
               | samtools view -@ 4 -bo ${singleton_bam} -) \\
        | awk '{{if(\$1 ~ /^@/){{print; next}} for(i=12;i<=NF;i++){{if(\$i ~ /^is:i:/){{split(\$i,a,\":\"); if(a[3]!=1){{print; break}}}}}}}}' \\
        | samtools view -@ 4 -bo ${hq_bam} -

    samtools fasta -@ 8 ${hq_bam}        | gzip -9 > ${hq_fq}
    samtools fasta -@ 8 ${singleton_bam} | gzip -9 > ${singleton_fq}

    rm ${hq_bam} ${singleton_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed 's/samtools //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed 's/samtools //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
