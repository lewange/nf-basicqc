/*
========================================================================================
    KRAKEN2_BATCH Module
========================================================================================
    Kraken2 - Taxonomic sequence classification system
    Batch mode: processes all samples in a single job to avoid reloading the database
*/

process KRAKEN2_BATCH {
    tag "${sample_names.size()} samples"
    label 'process_kraken'
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    // Dynamic resources based on number of samples
    // CPUs: 4 per sample for parallel processing, minimum 8, capped by max_cpus
    cpus   { Math.min(params.max_cpus as int, Math.max(8, sample_names.size() * 4)) }
    // Memory: 64GB base for DB + 2GB per sample
    memory { "${64 + (sample_names.size() * 2)} GB" }
    // Time: 1h base + 15min per sample
    time   { "${1 + (int)((sample_names.size() * 15) / 60)} h" }

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
    // Run multiple samples in parallel; each gets fewer threads but they share the memory-mapped DB
    def n_samples = sample_names.size()
    def parallel_jobs = Math.min(n_samples, Math.max(1, (task.cpus / 2).intValue()))
    def threads_per_job = Math.max(1, (task.cpus / parallel_jobs).intValue())
    """
    # Function to process a single sample
    process_sample() {
        sample=\$1
        threads=\$2
        db=\$3

        echo "Processing sample: \${sample} with \${threads} threads"

        # Find the read files for this sample (handles both SE and PE)
        r1=\$(ls \${sample}*_1.fastq.gz 2>/dev/null || ls \${sample}*.fastq.gz 2>/dev/null | head -1)
        r2=\$(ls \${sample}*_2.fastq.gz 2>/dev/null || echo "")

        if [[ -n "\$r2" ]]; then
            # Paired-end
            kraken2 \\
                --db \$db \\
                --threads \$threads \\
                --report \${sample}.kraken2.report.txt \\
                --output \${sample}.kraken2.output.txt \\
                --gzip-compressed \\
                --memory-mapping \\
                --paired \\
                $args \\
                \$r1 \$r2
        else
            # Single-end
            kraken2 \\
                --db \$db \\
                --threads \$threads \\
                --report \${sample}.kraken2.report.txt \\
                --output \${sample}.kraken2.output.txt \\
                --gzip-compressed \\
                --memory-mapping \\
                $args \\
                \$r1
        fi

        echo "Completed sample: \${sample}"
    }

    # Process samples in parallel using bash background jobs
    # Database is memory-mapped and shared between processes
    echo "Running $parallel_jobs samples in parallel with $threads_per_job threads each"

    # Track job PIDs and their sample names for error reporting
    declare -A pids_to_samples
    failed_samples=()

    samples=($samples_list)
    for sample in "\${samples[@]}"; do
        process_sample "\$sample" $threads_per_job $db &
        pids_to_samples[\$!]=\$sample

        # Limit concurrent jobs
        while [[ \$(jobs -r -p | wc -l) -ge $parallel_jobs ]]; do
            # Wait for any job to finish and check its exit status
            for pid in "\${!pids_to_samples[@]}"; do
                if ! kill -0 \$pid 2>/dev/null; then
                    wait \$pid || failed_samples+=("\${pids_to_samples[\$pid]}")
                    unset pids_to_samples[\$pid]
                fi
            done
            sleep 1
        done
    done

    # Wait for all remaining jobs and check their exit status
    for pid in "\${!pids_to_samples[@]}"; do
        wait \$pid || failed_samples+=("\${pids_to_samples[\$pid]}")
    done

    # Report failures and exit with error if any samples failed
    if [[ \${#failed_samples[@]} -gt 0 ]]; then
        echo "ERROR: The following samples failed: \${failed_samples[*]}" >&2
        exit 1
    fi

    echo "All samples completed successfully"

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
