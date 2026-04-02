# nf-basicqc

A Nextflow pipeline for basic quality control analysis of Illumina FASTQ sequencing files.

## Overview

This pipeline performs quality control, contamination screening, and taxonomic classification on raw sequencing data. It aggregates results into an interactive MultiQC report.

**Key analyses:**
- **FastQC** - Sequence quality metrics (per-base quality, GC content, adapter detection)
- **FastQ Screen** - Multi-genome contamination screening
- **Kraken2** - Taxonomic classification against a custom mtDNA database
- **Sex determination** - Inferred sex from marker read mapping
- **SortMeRNA** - rRNA quantification (% rRNA in subsampled reads)
- **RiboDetector** - Deep-learning rRNA detection (% rRNA in subsampled reads)
- **Kraken2 (rRNA)** - Species ID from RiboDetector-extracted rRNA reads classified against a custom animal rRNA database (16S/18S/28S Metazoa)
- **MultiQC** - Aggregated interactive report
- **Consolidated summary table** - `qc_summary.tsv` with all key metrics per sample

## Requirements

- Nextflow ≥23.04.0
- Singularity or Docker
- FastQ Screen configuration file and genome database (if using FastQ Screen)
- Kraken2 database (if using Kraken2 mtDNA classification)
- SortMeRNA rRNA FASTA database directory (if using SortMeRNA)
- Custom animal rRNA Kraken2 database (if using rRNA-based species ID; build with `build_animal_rrna_db.sh`)

## Quick Start

```bash
# Minimal run (FastQC + MultiQC only)
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --skip_fastq_screen \
  --skip_kraken2

# Full pipeline
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --fastq_screen_conf /path/to/fastq_screen.conf \
  --kraken2_db /path/to/kraken2_db \
  -profile singularity
```

## Input

### Samplesheet (CSV)

```csv
sample,fastq_1,fastq_2,sample_name,species
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,SampleA,Homo sapiens
sample2,/path/to/sample2_R1.fastq.gz,,SampleB,Mus musculus
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample identifier |
| `fastq_1` | Yes | Path to R1/forward reads (gzipped) |
| `fastq_2` | No | Path to R2/reverse reads for paired-end |
| `sample_name` | No | Display name for reports |
| `species` | No | Species information for report grouping |

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to input samplesheet (CSV) |
| `--outdir` | Output directory |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastq_screen_conf` | - | FastQ Screen configuration file |
| `--kraken2_db` | - | Kraken2 mtDNA database path |
| `--subsample_reads` | 1000000 | Reads to subsample per sample for all tools |
| `--kraken2_subsample` | - | Override subsample depth for Kraken2 + sex det only |
| `--rrna_subsample` | - | Override subsample depth for SortMeRNA + RiboDetector only |
| `--sex_markers_db` | - | Sex marker FASTA for sex determination |
| `--sortmerna_db` | - | Directory of rRNA FASTA files for SortMeRNA |
| `--sortmerna_index` | - | Pre-built SortMeRNA index dir (skips rebuild) |
| `--read_length` | 150 | Read length in bp (for RiboDetector) |
| `--rrna_kraken2_db` | - | Kraken2 database for rRNA species ID (uses RiboDetector-extracted reads) |
| `--project_name` | - | Project name for MultiQC header |
| `--application` | - | Application type for MultiQC header |

### Skip Options

| Parameter | Description |
|-----------|-------------|
| `--skip_fastqc` | Skip FastQC analysis |
| `--skip_fastq_screen` | Skip FastQ Screen |
| `--skip_kraken2` | Skip Kraken2 mtDNA classification |
| `--skip_sex_determination` | Skip sex determination |
| `--skip_sortmerna` | Skip SortMeRNA rRNA quantification |
| `--skip_ribodetector` | Skip RiboDetector rRNA quantification |
| `--skip_rrna_kraken2` | Skip rRNA-based Kraken2 species ID |

## Output

```
results/
├── fastqc/           # FastQC reports (HTML + ZIP)
├── fastq_screen/     # FastQ Screen reports
├── kraken2/          # Kraken2 mtDNA taxonomy reports
├── kraken2_rrna/     # Kraken2 SILVA SSU reports (rRNA species ID)
├── sortmerna/        # SortMeRNA logs (% rRNA metric only)
├── ribodetector/     # RiboDetector logs (+ rRNA reads if --rrna_kraken2_db set)
├── sex_determination/# Sex determination results
├── summary/
│   └── qc_summary.tsv  # Consolidated metrics for all samples
├── multiqc/
│   └── basicqc_multiqc_report.html  # Main report
└── pipeline_info/    # Execution reports
```

### MultiQC report sections

In addition to the standard FastQC / FastQ Screen / Kraken2 interactive plots, the report includes custom sections:

| Section | Type | Description |
|---------|------|-------------|
| Kraken2 Species Identification | Table | Top species/genus, % mtDNA, reads used, non-Metazoa warning |
| % Mitochondrial DNA | Bargraph | % reads classified against the mtDNA database per sample |
| FastQ Screen Summary | Table | Top-mapped genome and its unique-mapping % per sample |
| Sex Determination | Table | Inferred sex and confidence score |
| % rRNA (SortMeRNA) | Bargraph | % rRNA detected by SortMeRNA per sample |
| % rRNA (RiboDetector) | Bargraph | % rRNA detected by RiboDetector per sample |
| rRNA Species ID (SILVA) | Table | Top species/genus from SILVA Kraken2, reads used |

> **Non-Metazoa warning**: samples where the top Kraken2 hit is not in Metazoa (e.g. fungal contamination) are flagged with `Non-Metazoa (Fungi)` in the Warning column.

### Summary table columns

`qc_summary.tsv` contains (columns added conditionally based on which modules ran):

| Column | Source |
|--------|--------|
| `sample_name`, `expected_species` | Samplesheet |
| `total_reads`, `read_length`, `percent_gc`, `percent_duplicates` | FastQC |
| `percent_mtdna`, `top_genus`, `percent_top_genus`, `top_species`, `percent_top_species` | Kraken2 mtDNA |
| `inferred_sex`, `sex_confidence` | Sex determination |
| `sortmerna_pct_rrna` | SortMeRNA |
| `ribodetector_pct_rrna` | RiboDetector |
| `rrna_top_genus`, `rrna_pct_top_genus`, `rrna_top_species`, `rrna_pct_top_species` | Kraken2 rRNA |

## Profiles

### Execution environment

```bash
-profile singularity   # Use Singularity containers
-profile docker        # Use Docker containers
-profile conda         # Use Conda environments
-profile test          # Test with reduced resources
```

### Experiment-type presets

Combine with an execution profile to set sensible defaults for each assay type:

```bash
-profile singularity,wgs        # WGS: subsample_reads = 5 M (rRNA/mtDNA expected low)
-profile singularity,rnaseq     # RNA-seq: subsample_reads = 1 M (default)
-profile singularity,amplicon   # Amplicon: 1 M reads, skips SortMeRNA + RiboDetector
```

| Preset | `subsample_reads` | Notes |
|--------|------------------|-------|
| `wgs` | 5,000,000 | Higher depth needed to reliably detect low-level rRNA/mtDNA |
| `rnaseq` | 1,000,000 | Default depth; rRNA contamination typically detectable at 1 M |
| `amplicon` | 1,000,000 | rRNA/SortMeRNA not meaningful for amplicon data; both are skipped |

### SLURM Configuration

For SLURM cluster execution, provide a custom config file with `-c`:

```bash
-profile singularity -c /path/to/slurm.config
```

See `conf/slurm.config.example` for a template.

## Examples

```bash
# FastQC only
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --skip_fastq_screen \
  --skip_kraken2

# FastQC + FastQ Screen
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --fastq_screen_conf /path/to/fastq_screen.conf \
  --skip_kraken2

# Full pipeline on SLURM
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --fastq_screen_conf /path/to/fastq_screen.conf \
  --kraken2_db /path/to/kraken2_db \
  --project_name "MyProject" \
  -profile singularity \
  -c /path/to/slurm.config

# Resume interrupted run
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  -resume
```

## Resource Requirements

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| FastQC | 4 | 4 GB | 4h |
| FastQ Screen | 8 | 16 GB | 8h |
| Kraken2 (mtDNA) | 8 | 32 GB | 4h |
| Kraken2 (rRNA/SILVA) | 8 | 32 GB | 4h |
| SortMeRNA | 8 | 16 GB | 4h |
| RiboDetector | 8 | 32 GB | 4h |
| MultiQC | 2 | 8 GB | 2h |

## Running on SLURM

### Production Runs

Use `submit_pipeline.sh` for production runs. The fourth argument (application type) automatically selects the matching experiment preset:

```bash
# Usage
sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name] [application]

# Examples
sbatch submit_pipeline.sh inputs/CGLZOO_18.csv results/CGLZOO_18 CGLZOO_18 "WGS"
sbatch submit_pipeline.sh inputs/CGLZOO_01.csv results/CGLZOO_01 CGLZOO_01 "RNA-seq"
sbatch submit_pipeline.sh inputs/CGLZOO_02.csv results/CGLZOO_02 CGLZOO_02 "Amplicon"
```

Recognised application strings → profiles: `WGS` → `wgs`, `RNA-seq`/`RNAseq` → `rnaseq`, `Amplicon` → `amplicon`. Unrecognised values fall back to the 1 M default.

### Testing

Use `test/submit_tests.sh` to run test pipelines:

```bash
sbatch test/submit_tests.sh --full          # Full pipeline
sbatch test/submit_tests.sh --fastqc_only   # FastQC only
sbatch test/submit_tests.sh --kraken-only   # Kraken2 mtDNA only
sbatch test/submit_tests.sh --kraken-fresh  # Kraken2 no resume
sbatch test/submit_tests.sh --sortmerna     # SortMeRNA rRNA
sbatch test/submit_tests.sh --sex           # Sex determination
sbatch test/submit_tests.sh --rrna_kraken2  # rRNA Kraken2 (SILVA SSU)
```

## License

This project is licensed under the MIT License.

## Author

CryoZoo Project
