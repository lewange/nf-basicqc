/*
========================================================================================
    KRAKEN2_BATCH Module
========================================================================================
    Kraken2 - Taxonomic sequence classification system
    Batch mode: processes all samples in a single job to avoid reloading the database
*/

process KRAKEN2_BATCH {
    tag "all_samples"
    label 'process_kraken'
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    path(reads)
    path(db)
    val(sample_names)

    output:
    path("*.kraken2.report.txt"), emit: reports
    path("*.kraken2.output.txt"), emit: outputs, optional: true
    path("versions.yml")        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Convert sample names to a bash array
    def samples_list = sample_names.join(' ')
    """
    # Array of sample names
    samples=($samples_list)

    # Process each sample
    for sample in "\${samples[@]}"; do
        echo "Processing sample: \${sample}"

        # Find the read files for this sample (handles both SE and PE)
        r1=\$(ls \${sample}*_1.fastq.gz 2>/dev/null || ls \${sample}*.fastq.gz 2>/dev/null | head -1)
        r2=\$(ls \${sample}*_2.fastq.gz 2>/dev/null || echo "")

        if [[ -n "\$r2" ]]; then
            # Paired-end
            kraken2 \\
                --db $db \\
                --threads $task.cpus \\
                --report \${sample}.kraken2.report.txt \\
                --output \${sample}.kraken2.output.txt \\
                --gzip-compressed \\
                --paired \\
                $args \\
                \$r1 \$r2
        else
            # Single-end
            kraken2 \\
                --db $db \\
                --threads $task.cpus \\
                --report \${sample}.kraken2.report.txt \\
                --output \${sample}.kraken2.output.txt \\
                --gzip-compressed \\
                $args \\
                \$r1
        fi

        echo "Completed sample: \${sample}"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -n1 | sed 's/Kraken version //')
    END_VERSIONS
    """

    stub:
    def samples_list = sample_names.join(' ')
    """
    samples=($samples_list)
    for sample in "\${samples[@]}"; do
        touch \${sample}.kraken2.report.txt
        touch \${sample}.kraken2.output.txt
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: 2.1.3
    END_VERSIONS
    """
}
