/*
========================================================================================
    SUMMARIZE_FASTQ_SCREEN Module
========================================================================================
    Parses FastQ Screen .txt output files and generates a MultiQC custom content
    table with the top-mapped genome and its percentage for each sample.

    FastQ Screen runs per FLI (per read-pair), so this module aggregates results
    to the sample level using the fli→sample_name mapping.
*/

process SUMMARIZE_FASTQ_SCREEN {
    label 'process_single'
    executor 'local'

    input:
    path(screen_files)   // *_screen.txt files collected across all samples
    val(sample_info)     // List of maps: fli, sample_name, species

    output:
    path("fastq_screen_summary_mqc.txt"), emit: summary

    script:
    def fli_to_display = sample_info.collectEntries { info ->
        def display = info.sample_name
        if (info.species) {
            display = "${display} (${info.species})"
        }
        [(info.get('fli', info.sample_name)): display]
    }
    def fli_to_display_json = groovy.json.JsonOutput.toJson(fli_to_display)
    """
    #!/usr/bin/env python3
    import re, json, glob
    from collections import defaultdict

    fli_to_display = json.loads('${fli_to_display_json}')

    # Per-sample aggregation: accumulate total reads and per-genome mapped reads
    # across all FLIs/read-files belonging to the same sample.
    sample_genome_mapped = defaultdict(lambda: defaultdict(float))  # display -> genome -> sum(pct_mapped * reads)
    sample_reads         = defaultdict(float)                       # display -> sum(reads_processed)

    for screen_file in sorted(glob.glob("*_screen.txt")):
        # Derive the FLI from filename: strip _1_screen.txt / _2_screen.txt / _screen.txt
        fli = re.sub(r'(_[12])?_screen\\.txt\$', '', screen_file)
        display = fli_to_display.get(fli, fli)

        reads_processed = None
        genome_pcts = {}   # genome -> % mapped (one_hit + multi_hit for this library)

        with open(screen_file) as f:
            for line in f:
                line = line.rstrip('\\n')
                if line.startswith('#') or not line.strip():
                    continue
                if line.startswith('%Hit_no_genomes') or line.startswith('%Hit_no_libraries'):
                    continue
                parts = line.split('\\t')
                if parts[0] in ('Genome', 'Library'):
                    continue  # header row
                if len(parts) < 8:
                    continue
                library     = parts[0]
                try:
                    n_reads = int(parts[1])
                    pct_one  = float(parts[5])   # %One_hit_one_library
                    pct_mult = float(parts[7])   # %Multiple_hits_one_library
                except (ValueError, IndexError):
                    continue
                if reads_processed is None:
                    reads_processed = n_reads
                pct_mapped = pct_one + pct_mult
                genome_pcts[library] = pct_mapped

        if reads_processed and genome_pcts:
            sample_reads[display] += reads_processed
            for genome, pct in genome_pcts.items():
                sample_genome_mapped[display][genome] += pct * reads_processed

    # Write MultiQC custom content table
    with open("fastq_screen_summary_mqc.txt", 'w') as f:
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'FastQ Screen Summary'\\n")
        f.write("# description: 'Top-mapped genome and percentage from FastQ Screen'\\n")
        f.write("# pconfig:\\n")
        f.write("#     id: 'fastq_screen_summary_table'\\n")
        f.write("#     namespace: 'FastQ Screen'\\n")
        f.write("# headers:\\n")
        f.write("#     top_genome:\\n")
        f.write("#         title: 'Top Genome'\\n")
        f.write("#         description: 'Genome with the highest unique-mapping rate in FastQ Screen'\\n")
        f.write("#     pct_top_genome:\\n")
        f.write("#         title: '% Top Genome'\\n")
        f.write("#         description: 'Percentage of reads uniquely mapping to the top genome'\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.1f}'\\n")
        f.write("Sample\\ttop_genome\\tpct_top_genome\\n")

        for display in sorted(sample_reads):
            total = sample_reads[display]
            if total == 0:
                continue
            genomes = sample_genome_mapped[display]
            top_genome = max(genomes, key=lambda g: genomes[g])
            pct_top    = genomes[top_genome] / total
            f.write(f"{display}\\t{top_genome}\\t{pct_top:.1f}\\n")

    print(f"FastQ Screen summary: {len(sample_reads)} samples")
    """
}
