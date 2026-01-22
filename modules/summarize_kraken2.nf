/*
========================================================================================
    SUMMARIZE_KRAKEN2 Module
========================================================================================
    Parses Kraken2 reports and generates a summary file for MultiQC with:
    - Top species per sample
    - Percent of reads assigned to top species
    - Percent of reads classified (= % mitochondrial when using mtDNA database)
*/

process SUMMARIZE_KRAKEN2 {
    label 'process_single'
    executor 'local'

    input:
    path(reports)
    val(sample_info)  // List of maps with: sample_name, species

    output:
    path("kraken2_top_species_mqc.txt"), emit: summary

    script:
    // Build a lookup map from sample_name to display name (with species)
    def name_lookup = sample_info.collectEntries { info ->
        def display = info.sample_name
        if (info.species) {
            display = "${display} (${info.species})"
        }
        [(info.sample_name): display]
    }
    def name_lookup_json = groovy.json.JsonOutput.toJson(name_lookup)
    """
    #!/usr/bin/env python3
    import os
    import json
    import glob

    # Sample name lookup from samplesheet metadata
    name_lookup = json.loads('${name_lookup_json}')

    # Parse all kraken2 reports and extract top species
    results = []

    for report_file in glob.glob("*.kraken2.report.txt"):
        # Extract sample ID from filename (remove .kraken2.report.txt)
        sample_id = report_file.replace(".kraken2.report.txt", "")

        # Get display name from lookup, fall back to sample_id
        display_name = name_lookup.get(sample_id, sample_id)

        top_species = "Unknown"
        top_percent = 0.0
        percent_classified = 0.0
        percent_unclassified = 0.0

        with open(report_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) >= 6:
                    percent = float(parts[0].strip())
                    rank = parts[3].strip()
                    taxon = parts[5].strip()

                    # Get unclassified percentage (rank 'U')
                    if rank == 'U':
                        percent_unclassified = percent

                    # Look for species rank (S) with highest percentage
                    if rank == 'S' and percent > top_percent:
                        top_percent = percent
                        top_species = taxon

        # Calculate % classified (= % mitochondrial for mtDNA database)
        percent_classified = 100.0 - percent_unclassified

        results.append((display_name, top_species, top_percent, percent_classified))

    # Write MultiQC custom content file
    # This format adds columns to the General Stats table
    with open("kraken2_top_species_mqc.txt", 'w') as f:
        # Header with MultiQC configuration
        f.write("# plot_type: 'generalstats'\\n")
        f.write("# pconfig:\\n")
        f.write("#     - percent_classified:\\n")
        f.write("#         title: '% mtDNA'\\n")
        f.write("#         description: 'Percent of reads classified by Kraken2 (= % mitochondrial reads when using mtDNA database)'\\n")
        f.write("#         max: 100\\n")
        f.write("#         min: 0\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.2f}'\\n")
        f.write("#         scale: 'Blues'\\n")
        f.write("#     - percent_top_species:\\n")
        f.write("#         title: '% Top Species'\\n")
        f.write("#         description: 'Percent of reads assigned to top species by Kraken2'\\n")
        f.write("#         max: 100\\n")
        f.write("#         min: 0\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.1f}'\\n")
        f.write("#         scale: 'RdYlGn'\\n")
        f.write("#     - top_species:\\n")
        f.write("#         title: 'Top Species'\\n")
        f.write("#         description: 'Most abundant species detected by Kraken2'\\n")
        f.write("#         scale: False\\n")
        f.write("Sample\\tpercent_classified\\tpercent_top_species\\ttop_species\\n")

        for display_name, top_species, top_percent, percent_classified in results:
            f.write(f"{display_name}\\t{percent_classified:.2f}\\t{top_percent:.2f}\\t{top_species}\\n")

    print(f"Processed {len(results)} Kraken2 reports")
    """
}
