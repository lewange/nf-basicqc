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
- **Fixed MultiQC sample naming bugs** in `modules/prepare_multiqc_config.nf`:
  - **Bug 1**: MultiQC only found 1 report per module despite multiple files present
    - Root cause: Default `fn_clean_trim` patterns were stripping sample names at underscores, reducing `HFYMJDSXC_1_8bp-UDP0032` to just `HFYMJDSXC`
    - Fix: Added `fn_clean_trim: []` to disable aggressive default trimming
  - **Bug 2**: Sample names not showing as "SampleID (Species)" format
    - Root cause: `extra_fn_clean_exts` was removing `_1`/`_2` suffixes, causing R1/R2 deduplication; also patterns didn't account for R1/R2
    - Fix: Changed to `fn_clean_exts` (replaces defaults) and updated `sample_names_replace` to generate separate patterns for `_1`, `_2`, and base name
  - Expected result: Samples now display as `BB1523_R1 (Callithrix geoffroyi)`, `BB1523_R2 (...)`, etc.

- **Optimized Kraken2 batch processing** in `modules/kraken2_batch.nf`:
  - **Problem**: Low CPU utilization, sequential sample processing taking 25+ min for 2 samples
  - **Solution**: Parallel sample processing using bash background jobs
    - Added `--memory-mapping` flag so multiple kraken2 processes share the DB
    - Samples now run concurrently instead of sequentially
    - Proper error handling tracks failed samples and reports them
  - **Dynamic SLURM resources** based on sample count:
    - CPUs: 4 per sample (min 8, capped by `max_cpus`)
    - Memory: 64GB base + 2GB per sample
    - Time: 1h base + 15min per sample
  - Updated `nextflow.config` to use dynamic resources from module

---

## Open Tasks
- Verify MultiQC fix works after pipeline re-run
- Verify Kraken2 parallel processing speedup

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
