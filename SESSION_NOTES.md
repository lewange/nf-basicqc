# Session Notes

## Current Status (2026-04-02)

### Working Features
- **Pipeline runs successfully** with all modules (FastQC, FastQ Screen, Kraken2, Sex Determination, SortMeRNA, RiboDetector, MultiQC)
- **Summary table output** (`results/summary/qc_summary.tsv`) - consolidated TSV with all QC metrics
- **Kraken2 plot without unclassified** - modified reports exclude unclassified reads, percentages recalculated
- **SortMeRNA rRNA quantification** - % rRNA in subsampled reads, MultiQC general stats integration
- **RiboDetector rRNA quantification** - % rRNA via deep learning, MultiQC general stats integration
- **rRNA Kraken2 (SILVA SSU)** - species ID from RiboDetector-extracted rRNA reads; custom animal rRNA database (16S/12S/18S/28S Metazoa)
- **Unified subsampling** - single `subsample_reads` (1M default) shared by all tools; separate override params available
- **MultiQC bargraph plots** - % mtDNA (per sample), % rRNA SortMeRNA, % rRNA RiboDetector
- **Non-Metazoa flagging** - Kraken2 top hits outside Metazoa flagged in Warning column
- **Reads-used columns** - number of reads submitted to each Kraken2 run shown in MultiQC tables
- **FastQ Screen summary** - top genome + % uniquely mapped per sample as General Stats columns
- **Experiment-type presets** - `wgs`, `rnaseq`, `amplicon` Nextflow profiles
- **Submit script auto-profile** - application argument maps to preset profile automatically

### rRNA Pipeline (updated 2026-04-02)
- RiboDetector saves rRNA reads when `--rrna_kraken2_db` is provided (`save_rrna=true`)
- KRAKEN2_RRNA runs on RiboDetector rRNA reads (not SortMeRNA)
- Custom animal rRNA database (`k2_animal_rrna`): 16S/12S (biomol_rrna, Metazoa), 18S + 28S (RefSeq NR_, Metazoa)
- Build script: `data/.../build_animal_rrna_db.sh`

### MultiQC Bargraph Format Notes
- Bargraph files must NOT have a file-level `# id:` header (let MultiQC derive from filename)
- Section names must be unique across the report
- `$` in regex patterns inside Nextflow `"""..."""` heredocs must be written as `\$`
- Triple-double-quote Python docstrings inside Nextflow heredocs will close the block — use `#` comments instead

## Key Files

### Pipeline Configuration
- `nextflow.config` - Container versions, resource settings, params
- `submit_pipeline.sh` - Production run script
- `test/submit_tests.sh` - Test scripts (use `--rrna_kraken2` to test new feature)

### Modules
- `modules/fastqc.nf` - FastQC
- `modules/fastq_screen.nf` - FastQ Screen
- `modules/kraken2.nf` - Kraken2 (used for both mtDNA and rRNA classification via alias)
- `modules/sex_determination.nf` - Sex determination + SUMMARIZE_SEX
- `modules/sortmerna.nf` - SortMeRNA (index + classification)
- `modules/ribodetector.nf` - RiboDetector (optionally saves rRNA reads when `save_rrna=true`)
- `modules/summarize_kraken2.nf` - SUMMARIZE_KRAKEN2 + SUMMARIZE_RRNA_KRAKEN2; bargraph + table
- `modules/summarize_rrna.nf` - SUMMARIZE_RRNA; generalstats + bargraph plots
- `modules/summarize_fastq_screen.nf` - SUMMARIZE_FASTQ_SCREEN; top genome + % mapped table
- `modules/summarize_results.nf` - Consolidated QC summary table (qc_summary.tsv)
- `modules/prepare_multiqc_config.nf` - Generates MultiQC config YAML
- `modules/multiqc.nf` - MultiQC

### Custom Content Files Generated
- `kraken2_top_species_mqc.txt` - Table: top species/genus, % mtDNA, reads used, non-Metazoa warning
- `kraken2_pct_mtdna_mqc.txt` - Bargraph: % mitochondrial DNA per sample
- `*_classified.kraken2.report.txt` - Modified Kraken reports without unclassified
- `sex_determination_mqc.txt` - Sex determination results
- `rrna_pct_mqc.txt` - General stats: % rRNA (SortMeRNA + RiboDetector)
- `rrna_pct_sortmerna_mqc.txt` - Bargraph: % rRNA (SortMeRNA)
- `rrna_pct_ribodetector_mqc.txt` - Bargraph: % rRNA (RiboDetector)
- `fastq_screen_summary_mqc.txt` - Table: top genome + % unique mapping per sample
- `rrna_kraken2_species_mqc.txt` - Table: rRNA species ID, reads used

## Databases
- Kraken2 mtDNA: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/kraken/k2_mtdna`
- Sex markers: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/sex_markers/all_sex_markers.fasta`
- SortMeRNA rRNA FASTAs: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/rRNA_indices`
- SILVA SSU Kraken2 (to download): `https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz`
  - Suggested path: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/kraken/k2_silva_ssu`
  - Update `SILVA_DB` in `test/submit_tests.sh` once downloaded

### Callithrix species in mtDNA database
| Species | Taxon ID | K-mers |
|---------|----------|--------|
| C. aurita | 57375 | 2,240 |
| C. geoffroyi | 52231 | 853 |
| C. kuhlii | 867363 | 572 |
| C. jacchus | 9483 | 476 |
| C. penicillata | 57378 | 350 |

~39% of Callithrix k-mers are shared at genus level.

## Commands

### Run tests
```bash
sbatch test/submit_tests.sh --full             # Full pipeline
sbatch test/submit_tests.sh --sortmerna        # SortMeRNA rRNA
sbatch test/submit_tests.sh --rrna_kraken2     # rRNA Kraken2 (SILVA SSU) — needs SILVA_DB set
sbatch test/submit_tests.sh --kraken-fresh     # Kraken2 only, no resume
```

### Production run
```bash
sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name] [application]
```

### Download SILVA SSU database
```bash
wget https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz
mkdir -p k2_silva_ssu && tar -xzf 16S_Silva138_20200326.tgz -C k2_silva_ssu/
```

## Recent Commits
```
Add FastQ Screen summary, bargraph plots, non-Metazoa flag, reads-used columns, experiment presets
Add % mtDNA and % rRNA bargraph plots to MultiQC report
Run rRNA tools per-sample; unify subsampling into single step
Shift rRNA Kraken2 to RiboDetector reads; add custom animal rRNA database support
```
