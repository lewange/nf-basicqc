# nf-basicqc

A Nextflow pipeline for basic quality control analysis of Illumina FASTQ sequencing files.

## Overview

This pipeline performs quality control, contamination screening, and taxonomic classification on raw sequencing data. It aggregates results into an interactive MultiQC report.

**Key analyses:**
- **FastQC** - Sequence quality metrics (per-base quality, GC content, adapter detection)
- **FastQ Screen** - Multi-genome contamination screening
- **Kraken2** - Taxonomic classification
- **MultiQC** - Aggregated interactive report

## Requirements

- Nextflow ≥23.04.0
- Singularity or Docker
- FastQ Screen configuration file and genome database (if using FastQ Screen)
- Kraken2 database (if using Kraken2)

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
| `--kraken2_db` | - | Kraken2 database path |
| `--kraken2_subsample` | 5000000 | Number of reads to subsample for Kraken2 |
| `--project_name` | - | Project name for MultiQC header |
| `--application` | - | Application type for MultiQC header |

### Skip Options

| Parameter | Description |
|-----------|-------------|
| `--skip_fastqc` | Skip FastQC analysis |
| `--skip_fastq_screen` | Skip FastQ Screen |
| `--skip_kraken2` | Skip Kraken2 classification |

## Output

```
results/
├── fastqc/           # FastQC reports (HTML + ZIP)
├── fastq_screen/     # FastQ Screen reports
├── kraken2/          # Kraken2 taxonomy reports
├── multiqc/
│   └── basicqc_multiqc_report.html  # Main report
└── pipeline_info/    # Execution reports
```

## Profiles

```bash
-profile singularity   # Use Singularity containers
-profile docker        # Use Docker containers
-profile conda         # Use Conda environments
-profile test          # Test with reduced resources
```

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
| Kraken2 | 8 | 64 GB | 12h |
| MultiQC | 2 | 8 GB | 2h |

## Running on SLURM

### Production Runs

Use `submit_pipeline.sh` for production runs:

```bash
# Usage: sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name]

# Example
sbatch submit_pipeline.sh inputs/CGLZOO_01.csv results/CGLZOO_01 CGLZOO_01
```

### Testing

Use `test/submit_tests.sh` to run test pipelines:

```bash
# Full pipeline test
sbatch test/submit_tests.sh --full

# FastQC only test
sbatch test/submit_tests.sh --fastqc_only

# Kraken2 only test
sbatch test/submit_tests.sh --kraken-only

# Fresh Kraken2 test (no resume)
sbatch test/submit_tests.sh --kraken-fresh
```

## License

This project is licensed under the MIT License.

## Author

CryoZoo Project
