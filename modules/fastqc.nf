/*
========================================================================================
    FASTQC Module
========================================================================================
    FastQC - A quality control tool for high throughput sequence data
*/

process FASTQC {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.html"), emit: html
    tuple val(sample), path("*.zip") , emit: zip
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: sample
    """
    fastqc \\
        --threads $task.cpus \\
        --outdir . \\
        $args \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: sample
    """
    touch ${prefix}_1_fastqc.html
    touch ${prefix}_1_fastqc.zip
    touch ${prefix}_2_fastqc.html
    touch ${prefix}_2_fastqc.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}
