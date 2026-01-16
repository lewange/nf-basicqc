/*
========================================================================================
    KRAKEN2 Module
========================================================================================
    Kraken2 - Taxonomic sequence classification system
    Based on nf-core/modules kraken2 module
*/

process KRAKEN2 {
    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path db

    output:
    tuple val(sample), path("*.kraken2.classifiedreads.txt"), emit: classified_reads, optional: true
    tuple val(sample), path("*.kraken2.unclassifiedreads.txt"), emit: unclassified_reads, optional: true
    tuple val(sample), path("*.kraken2.report.txt")         , emit: report
    tuple val(sample), path("*.kraken2.output.txt")         , emit: output
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: sample
    // Determine if paired-end or single-end
    def paired = reads instanceof List && reads.size() == 2 ? '--paired' : ''
    def classified_option = params.kraken2_save_reads ? "--classified-out ${prefix}.kraken2.classifiedreads.txt" : ''
    def unclassified_option = params.kraken2_save_reads ? "--unclassified-out ${prefix}.kraken2.unclassifiedreads.txt" : ''
    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --output ${prefix}.kraken2.output.txt \\
        --gzip-compressed \\
        $paired \\
        $classified_option \\
        $unclassified_option \\
        $args \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -n1 | sed 's/Kraken version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: sample
    """
    touch ${prefix}.kraken2.report.txt
    touch ${prefix}.kraken2.output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: 2.1.3
    END_VERSIONS
    """
}
