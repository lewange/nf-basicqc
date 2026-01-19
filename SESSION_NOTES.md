# Session Notes - nf-basicqc

## Project Overview
Nextflow pipeline for basic quality control analysis of Illumina FASTQ sequencing files.

**Key components:**
- FastQC - Sequence quality metrics
- FastQ Screen - Multi-genome contamination screening
- Kraken2 - Taxonomic classification
- MultiQC - Aggregated interactive report

## Current State

**Last updated:** 2026-01-19

**Git status:** On branch `main`, modified `test/submit_tests.sh`

**Recent commits:**
- `4c01769` modified submit_tests.sh
- `46cddb8` updated README
- `9af7b1e` added slurm submission scripts and updated README
- `f794e25` Adding README
- `3e765ca` Initial commit: BasicQC Nextflow pipeline

## Key Files

| File | Purpose |
|------|---------|
| `main.nf` | Main Nextflow workflow |
| `nextflow.config` | Pipeline configuration |
| `submit_pipeline.sh` | SLURM submission script for production |
| `test/submit_tests.sh` | SLURM submission script for testing |
| `conf/` | Configuration profiles |
| `modules/` | Nextflow process modules |

## Session History

### Session 1 - 2026-01-19
- Initial session notes file created
- Project structure reviewed

---

## Open Tasks
_None currently_

## Notes for Next Session
_Add notes here before ending each session_

---

## Quick Reference

```bash
# Run minimal pipeline (FastQC + MultiQC)
nextflow run main.nf --input samplesheet.csv --outdir results --skip_fastq_screen --skip_kraken2

# Run full pipeline
nextflow run main.nf --input samplesheet.csv --outdir results --fastq_screen_conf /path/to/conf --kraken2_db /path/to/db -profile singularity

# Resume interrupted run
nextflow run main.nf --input samplesheet.csv --outdir results -resume
```
