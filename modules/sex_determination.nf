/*
========================================================================================
    SEX_DETERMINATION Module
========================================================================================
    Determines genetic sex from sequencing reads by matching to sex-specific markers:
    - Mammals: Y chromosome markers (SRY, AMELY, DDX3Y, ZFY)
    - Birds: W chromosome markers (CHD1W, NIPBLW, SPINW)

    Uses the species classification from Kraken2 to interpret results correctly.
*/

process SEX_DETERMINATION {
    tag "$sample"
    label 'process_low'
    publishDir "${params.outdir}/sex_determination", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path(marker_db)
    val(species_class)  // 'mammal', 'bird', or 'unknown' from Kraken2

    output:
    tuple val(sample), path("${sample}.sex_markers.txt"), emit: results
    path "versions.yml", emit: versions

    script:
    def input_reads = reads instanceof List ? reads[0] : reads
    """
    #!/usr/bin/env python3
    import subprocess
    import gzip
    import os
    import re
    from collections import defaultdict

    sample = "${sample}"
    reads_file = "${input_reads}"
    marker_fasta = "${marker_db}"
    species_class = "${species_class}".lower()

    # Load marker sequences and metadata
    markers = {}
    current_id = None
    current_seq = []
    current_meta = {}

    with open(marker_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id and current_seq:
                    markers[current_id] = {
                        'seq': ''.join(current_seq).upper(),
                        'meta': current_meta
                    }
                current_id = line[1:].split()[0]
                current_seq = []
                current_meta = {}
                # Parse metadata from header
                if 'sex=male' in line:
                    current_meta['sex'] = 'male'
                elif 'sex=female' in line:
                    current_meta['sex'] = 'female'
                if 'class=mammal' in line:
                    current_meta['class'] = 'mammal'
                elif 'class=bird' in line:
                    current_meta['class'] = 'bird'
                # Extract gene name
                gene_match = re.search(r'gene=(\\S+)', line)
                if gene_match:
                    current_meta['gene'] = gene_match.group(1)
                else:
                    current_meta['gene'] = current_id.split('_')[0]
            else:
                current_seq.append(line.strip())
        if current_id and current_seq:
            markers[current_id] = {
                'seq': ''.join(current_seq).upper(),
                'meta': current_meta
            }

    print(f"Loaded {len(markers)} marker sequences")

    # Create k-mer index from markers (k=31 for specificity)
    K = 31
    kmer_to_marker = defaultdict(set)

    for marker_id, data in markers.items():
        seq = data['seq']
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i+K]
            if 'N' not in kmer:
                kmer_to_marker[kmer].add(marker_id)
                # Also add reverse complement
                rc = seq[i:i+K].translate(str.maketrans('ATCG', 'TAGC'))[::-1]
                kmer_to_marker[rc].add(marker_id)

    print(f"Built k-mer index with {len(kmer_to_marker)} unique {K}-mers")

    # Scan reads for marker k-mers
    marker_hits = defaultdict(int)
    total_reads = 0

    def process_read(seq):
        seq = seq.upper()
        hits = set()
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i+K]
            if kmer in kmer_to_marker:
                hits.update(kmer_to_marker[kmer])
        return hits

    # Read FASTQ file
    open_func = gzip.open if reads_file.endswith('.gz') else open
    with open_func(reads_file, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()

            total_reads += 1
            hits = process_read(seq)
            for marker_id in hits:
                marker_hits[marker_id] += 1

    print(f"Processed {total_reads} reads")

    # Aggregate results by sex and class
    male_mammal_hits = 0
    female_bird_hits = 0
    male_genes = defaultdict(int)
    female_genes = defaultdict(int)

    for marker_id, count in marker_hits.items():
        meta = markers[marker_id]['meta']
        sex = meta.get('sex', 'unknown')
        cls = meta.get('class', 'unknown')
        gene = meta.get('gene', marker_id)

        if sex == 'male' and cls == 'mammal':
            male_mammal_hits += count
            male_genes[gene] += count
        elif sex == 'female' and cls == 'bird':
            female_bird_hits += count
            female_genes[gene] += count

    # Determine sex based on species class
    inferred_sex = "Unknown"
    confidence = "low"
    evidence = ""

    if species_class == 'mammal' or (male_mammal_hits > 0 and female_bird_hits == 0):
        # Mammalian logic: Y markers = male
        if male_mammal_hits >= 10:
            inferred_sex = "Male"
            confidence = "high" if male_mammal_hits >= 50 else "medium"
            evidence = f"Y-chromosome markers: {male_mammal_hits} hits"
        elif male_mammal_hits > 0:
            inferred_sex = "Male"
            confidence = "low"
            evidence = f"Y-chromosome markers: {male_mammal_hits} hits (low count)"
        else:
            inferred_sex = "Female"
            confidence = "medium"
            evidence = "No Y-chromosome markers detected"

    elif species_class == 'bird' or (female_bird_hits > 0 and male_mammal_hits == 0):
        # Bird logic: W markers = female
        if female_bird_hits >= 10:
            inferred_sex = "Female"
            confidence = "high" if female_bird_hits >= 50 else "medium"
            evidence = f"W-chromosome markers: {female_bird_hits} hits"
        elif female_bird_hits > 0:
            inferred_sex = "Female"
            confidence = "low"
            evidence = f"W-chromosome markers: {female_bird_hits} hits (low count)"
        else:
            inferred_sex = "Male"
            confidence = "medium"
            evidence = "No W-chromosome markers detected"

    else:
        # Unknown class - report based on any hits
        if male_mammal_hits > female_bird_hits and male_mammal_hits >= 5:
            inferred_sex = "Male (mammal)"
            confidence = "low"
            evidence = f"Y markers: {male_mammal_hits}, W markers: {female_bird_hits}"
        elif female_bird_hits > male_mammal_hits and female_bird_hits >= 5:
            inferred_sex = "Female (bird)"
            confidence = "low"
            evidence = f"W markers: {female_bird_hits}, Y markers: {male_mammal_hits}"
        else:
            inferred_sex = "Unknown"
            confidence = "none"
            evidence = f"Y markers: {male_mammal_hits}, W markers: {female_bird_hits}"

    # Write results
    with open(f"{sample}.sex_markers.txt", 'w') as f:
        f.write(f"Sample: {sample}\\n")
        f.write(f"Total reads scanned: {total_reads}\\n")
        f.write(f"Species class: {species_class}\\n")
        f.write(f"\\n")
        f.write(f"Inferred sex: {inferred_sex}\\n")
        f.write(f"Confidence: {confidence}\\n")
        f.write(f"Evidence: {evidence}\\n")
        f.write(f"\\n")
        f.write(f"=== Marker hits ===\\n")
        f.write(f"Mammalian Y markers: {male_mammal_hits}\\n")
        for gene, count in sorted(male_genes.items(), key=lambda x: -x[1]):
            f.write(f"  {gene}: {count}\\n")
        f.write(f"Bird W markers: {female_bird_hits}\\n")
        for gene, count in sorted(female_genes.items(), key=lambda x: -x[1]):
            f.write(f"  {gene}: {count}\\n")

    print(f"\\nResult: {inferred_sex} ({confidence} confidence)")
    print(f"Evidence: {evidence}")

    # Write versions
    import sys
    with open("versions.yml", "w") as vf:
        vf.write('"${task.process}":\\n')
        vf.write(f"    python: {sys.version.split()[0]}\\n")
    """
}

/*
========================================================================================
    SUMMARIZE_SEX Module
========================================================================================
    Aggregates sex determination results for MultiQC
*/

process SUMMARIZE_SEX {
    label 'process_single'
    executor 'local'

    input:
    path(results)
    val(sample_info)  // List of maps with: sample_name, species

    output:
    path("sex_determination_mqc.txt"), emit: summary

    script:
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
    import re

    name_lookup = json.loads('${name_lookup_json}')

    results = []

    for result_file in glob.glob("*.sex_markers.txt"):
        sample_id = result_file.replace(".sex_markers.txt", "")
        display_name = name_lookup.get(sample_id, sample_id)

        inferred_sex = "Unknown"
        confidence = "none"
        y_hits = 0
        w_hits = 0

        with open(result_file, 'r') as f:
            for line in f:
                if line.startswith("Inferred sex:"):
                    inferred_sex = line.split(":", 1)[1].strip()
                elif line.startswith("Confidence:"):
                    confidence = line.split(":", 1)[1].strip()
                elif line.startswith("Mammalian Y markers:"):
                    y_hits = int(line.split(":", 1)[1].strip())
                elif line.startswith("Bird W markers:"):
                    w_hits = int(line.split(":", 1)[1].strip())

        results.append({
            'sample': display_name,
            'sex': inferred_sex,
            'confidence': confidence,
            'y_hits': y_hits,
            'w_hits': w_hits
        })

    # Write MultiQC custom content file as a table section
    with open("sex_determination_mqc.txt", 'w') as f:
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'Sex Determination'\\n")
        f.write("# description: 'Genetic sex inferred from sex-specific markers'\\n")
        f.write("# pconfig:\\n")
        f.write("#     id: 'sex_determination_table'\\n")
        f.write("#     namespace: 'Sex'\\n")
        f.write("# headers:\\n")
        f.write("#     inferred_sex:\\n")
        f.write("#         title: 'Inferred Sex'\\n")
        f.write("#         description: 'Genetic sex based on marker analysis'\\n")
        f.write("#     sex_confidence:\\n")
        f.write("#         title: 'Confidence'\\n")
        f.write("#         description: 'Confidence level (high/medium/low)'\\n")
        f.write("#     y_marker_hits:\\n")
        f.write("#         title: 'Y Marker Hits'\\n")
        f.write("#         description: 'Reads matching Y-chromosome markers (mammals)'\\n")
        f.write("#         format: '{:,.0f}'\\n")
        f.write("#     w_marker_hits:\\n")
        f.write("#         title: 'W Marker Hits'\\n")
        f.write("#         description: 'Reads matching W-chromosome markers (birds)'\\n")
        f.write("#         format: '{:,.0f}'\\n")
        f.write("Sample\\tinferred_sex\\tsex_confidence\\ty_marker_hits\\tw_marker_hits\\n")

        for r in results:
            f.write(f"{r['sample']}\\t{r['sex']}\\t{r['confidence']}\\t{r['y_hits']}\\t{r['w_hits']}\\n")

    print(f"Processed {len(results)} sex determination results")
    """
}
