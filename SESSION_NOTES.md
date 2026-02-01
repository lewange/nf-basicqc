# Session Notes - 2026-01-27

## Current Status

### Working Features
- **Pipeline runs successfully** with all modules (FastQC, FastQ Screen, Kraken2, Sex Determination, MultiQC)
- **Summary table output** (`results/summary/qc_summary.tsv`) - consolidated TSV with all QC metrics
- **Kraken2 plot without unclassified** - modified reports exclude unclassified reads, percentages recalculated
- **MultiQC version updated to 1.33**

### Known Issues to Fix

1. **General Stats table not showing in MultiQC report**
   - Data IS collected (visible in `multiqc_data.json`)
   - Log shows "Found 2 General Statistics columns" but we have 5 defined
   - Tried: module_order, table_columns_visible, table_columns_placement configs
   - May be an issue with custom content YAML format or MultiQC parsing
   - File: `modules/summarize_kraken2.nf` generates `kraken2_top_species_mqc.txt`

2. **Kraken plot only shows 5 species**
   - Need to configure `top_n` in MultiQC config
   - File: `modules/prepare_multiqc_config.nf`

3. **Sex determination not showing in MultiQC**
   - Check if `sex_determination_mqc.txt` is being passed to MultiQC
   - File: `modules/sex_determination.nf` (SUMMARIZE_SEX process)

## Key Files

### Pipeline Configuration
- `nextflow.config` - Container versions, resource settings
- `submit_pipeline.sh` - Production run script
- `test/submit_tests.sh` - Test scripts

### MultiQC Integration
- `modules/prepare_multiqc_config.nf` - Generates MultiQC config YAML
- `modules/summarize_kraken2.nf` - Kraken2 summary + modified reports
- `modules/sex_determination.nf` - Sex determination + SUMMARIZE_SEX
- `modules/summarize_results.nf` - Consolidated QC summary table

### Custom Content Files Generated
- `kraken2_top_species_mqc.txt` - General stats columns (plot_type: 'generalstats')
- `*_classified.kraken2.report.txt` - Modified Kraken reports without unclassified
- `sex_determination_mqc.txt` - Sex determination results

## Databases
- Kraken2 mtDNA: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/kraken/k2_mtdna`
- Sex markers: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/sex_markers/all_sex_markers.fasta`

### Callithrix species in mtDNA database
| Species | Taxon ID | K-mers |
|---------|----------|--------|
| C. aurita | 57375 | 2,240 |
| C. geoffroyi | 52231 | 853 |
| C. kuhlii | 867363 | 572 |
| C. jacchus | 9483 | 476 |
| C. penicillata | 57378 | 350 |

~39% of Callithrix k-mers are shared at genus level.

## Next Steps

1. **Fix Kraken top_n** - Add to prepare_multiqc_config.nf:
   ```yaml
   kraken:
       top_n: 10
   ```

2. **Debug sex determination in MultiQC**
   - Check if SUMMARIZE_SEX.out.summary is being added to ch_multiqc_files
   - Verify sex_determination_mqc.txt format matches MultiQC expectations

3. **Fix General Stats table**
   - Try simplifying the pconfig format in kraken2_top_species_mqc.txt
   - Or create a separate custom table section instead of generalstats

## Recent Commits
```
26876d0 Update MultiQC to version 1.33
2128c76 Fix Python boolean syntax in summarize_results.nf
5f3130b Add consolidated QC summary table output
95c52ed Update MultiQC container to version 1.27
815e0ee Update submit_pipeline.sh for production runs
2af1214 Fix General Stats table visibility in MultiQC report
7743828 Enhance Kraken2 MultiQC reporting with genus stats and unclassified-free plot
```

## Commands

### Run tests
```bash
sbatch test/submit_tests.sh --full        # Full pipeline test
sbatch test/submit_tests.sh --kraken-fresh # Kraken2 only, no resume
```

### Production run
```bash
sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name] [application]
```

### Check MultiQC log
```bash
cat results/multiqc/*_multiqc_report_data/multiqc.log
```
