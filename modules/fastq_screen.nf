/*
========================================================================================
    FASTQ_SCREEN Module
========================================================================================
    FastQ Screen - A tool for multi-genome mapping and species identification
*/

process FASTQ_SCREEN {
    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/fastq_screen", mode: 'copy', saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(sample), path(reads)
    path config

    output:
    tuple val(sample), path("*_screen.txt") , emit: txt
    tuple val(sample), path("*_screen.html"), emit: html
    tuple val(sample), path("*_screen.png") , emit: png, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: sample
    // Handle both single-end and paired-end reads
    def read_files = reads instanceof List ? reads.join(' ') : reads
    """
    fastq_screen \\
        --threads $task.cpus \\
        --conf $config \\
        --outdir . \\
        --aligner bowtie2 \\
        $args \\
        $read_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_screen: \$(fastq_screen --version 2>&1 | sed 's/fastq_screen v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: sample
    """
    touch ${prefix}_1_screen.txt
    touch ${prefix}_1_screen.html
    touch ${prefix}_2_screen.txt
    touch ${prefix}_2_screen.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_screen: 0.16.0
    END_VERSIONS
    """
}
