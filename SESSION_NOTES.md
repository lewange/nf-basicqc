# Session Notes - nf-basicqc

## Project Overview
Nextflow pipeline for basic quality control analysis of Illumina FASTQ sequencing files.

**Key components:**
- FastQC - Sequence quality metrics
- FastQ Screen - Multi-genome contamination screening
- Kraken2 - Taxonomic classification
- MultiQC - Aggregated interactive report

## Current State

**Last updated:** 2026-01-22

**Git status:** On branch `main`

**Recent commits:**
- `9980064` Switch Kraken2 to per-sample processing with one FASTQ per sample
- `64fdd3c` changed test script
- `90fb542` Optimize Kraken2 batch processing with parallel execution
- `54b8015` Fix MultiQC sample naming and report detection
- `4c01769` modified submit_tests.sh

## Key Files

| File | Purpose |
|------|---------|
| `main.nf` | Main Nextflow workflow |
| `nextflow.config` | Pipeline configuration |
| `submit_pipeline.sh` | SLURM submission script for production |
| `test/submit_tests.sh` | SLURM submission script for testing |
| `conf/` | Configuration profiles |
| `modules/` | Nextflow process modules |
| `modules/sex_determination.nf` | Sex determination from genetic markers |
| `scripts/compare_kraken2_dbs.sh` | Compare full vs mtDNA Kraken2 databases |

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

### Session 2 - 2026-01-20
- **Added custom Kraken2 summary columns to MultiQC** in `modules/summarize_kraken2.nf`:
  - **Problem**: MultiQC's built-in Kraken2 module creates one column per species detected (e.g., "Callithrix_jacchus"), which only makes sense if all samples have the same top species
  - **Solution**: New `SUMMARIZE_KRAKEN2` module that:
    - Parses each `.kraken2.report.txt` file
    - Finds the species (rank 'S') with highest percent of reads
    - Outputs a MultiQC custom content file with two columns:
      - `% Top Species` - Percent of reads assigned to top species
      - `Top Species` - Name of the most abundant species
  - **Files modified**:
    - `modules/summarize_kraken2.nf` (new) - Parses reports and generates custom MultiQC content
    - `modules/prepare_multiqc_config.nf` - Hides default Kraken module columns
    - `main.nf` - Integrates new module into workflow
  - Expected result: Each sample shows its own top species with percentage in a clean format

---

### Session 3 - 2026-01-21
- **Confirmed** Kraken2 summary columns working correctly in MultiQC report
- **Reverted Kraken2 from batch to per-sample processing**:
  - **Problem**: Batch processing with parallel execution (Session 1) was slower than sequential due to 350GB database overhead
  - **Solution**: Switch back to per-sample KRAKEN2 jobs, but only process ONE FASTQ per sample_name
  - **Changes made**:
    - `main.nf`:
      - Switched from `KRAKEN2_BATCH` to per-sample `KRAKEN2` module
      - Added `parse_samplesheet_kraken2()` function that groups by `sample_name` and takes only the first FASTQ pair
      - Kraken2 now runs as separate jobs, one per unique sample_name
    - `modules/summarize_kraken2.nf`:
      - Updated metadata structure from `fli/sample_name/species` to `sample_name/species`
  - **Result**: Fewer Kraken2 jobs (one per sample_name instead of one per FASTQ), each loading DB independently

- **Started building filtered Kraken2 database** (step 2 of optimization):
  - **Goal**: Reduce database from 233GB to ~30-50GB for faster loading
  - **Approach**: Build new database containing only:
    - Tetrapods (taxid 32523) - mammals, birds, reptiles, amphibians (full species-level)
    - Bacteria (taxid 2) - for contamination detection
    - Fungi (taxid 4751) - for contamination detection
  - **Script**: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/kraken/build_filtered_db.sh`
  - **How it works**:
    1. Downloads taxonomy + nt library from NCBI (using `--use-ftp` to avoid rsync issues)
    2. Extracts all descendant taxids for target groups from `nodes.dmp`
    3. Filters `seqid2taxid.map` to only include sequences from target taxa
    4. Builds database from filtered sequences
  - **Status**: SLURM job submitted, nt library downloading (will take many hours)
  - **Output database**: `k2_tetrapod_bac_fungi` in kraken directory

---

### Session 4 - 2026-01-22
- **Pivoted Kraken2 database strategy** to mitochondrial-only approach:
  - **Previous approach** (Session 3): Download full nt database and filter for tetrapods + bacteria/fungi
    - Problem: nt download is 300+ GB, filtering is complex and slow
  - **New approach**: Build small mtDNA-only database (~1GB) for species identification
  - **Rationale**:
    - mtDNA is naturally enriched in all sequencing types (WGS, RNA-seq, ATAC-seq, methylation)
    - mtDNA has excellent species resolution (faster evolution than nuclear)
    - Much smaller database = faster loading and searches
    - Works for sample verification use case (detecting sample swaps)

- **Created mtDNA Kraken2 database build script**:
  - **Location**: `/Users/lucas/scratch/data/Illumina_CryoZoo/genomes/kraken/build_mtdna_db.sh`
  - Downloads all RefSeq mitochondrial genomes
  - Builds Kraken2 database (no filtering needed - non-tetrapod mtDNA won't match tetrapod samples)
  - **Output database**: `k2_mtdna`
  - **Status**: Building on cluster

- **Updated SUMMARIZE_KRAKEN2 module** to report % mtDNA:
  - Now reports three columns in MultiQC:
    - `% mtDNA` - Percent of reads classified (= mitochondrial content)
    - `% Top Species` - Percent assigned to top species
    - `Top Species` - Name of most abundant species
  - % mtDNA is useful QC metric on its own (indicates sample quality, cell lysis, etc.)

- **Created database comparison script**:
  - **Location**: `scripts/compare_kraken2_dbs.sh`
  - Runs Kraken2 with both full database and mtDNA database on same samples
  - Compares species assignments and reports agreement
  - Outputs CSV and detailed text report

- **Added sex determination feature**:
  - **Goal**: Detect genetic sex from sequencing reads for QC/sample verification
  - **Approach**: K-mer matching against sex-specific markers
    - Mammals: Y chromosome markers (SRY, AMELY, DDX3Y, ZFY, USP9Y)
    - Birds: W chromosome markers (CHD1W, NIPBLW, SPINW)
  - **Files created**:
    - `modules/sex_determination.nf` - SEX_DETERMINATION and SUMMARIZE_SEX processes
    - `/Users/lucas/scratch/data/Illumina_CryoZoo/genomes/sex_markers/build_sex_marker_db.sh` - Database builder
  - **Pipeline integration**:
    - New parameters: `--sex_markers_db`, `--skip_sex_determination`
    - Uses same subsampled reads as Kraken2 (efficient)
    - MultiQC columns: Sex, Sex Conf., Y Hits (hidden), W Hits (hidden)
  - **How it works**:
    1. Builds k-mer index (k=31) from sex marker sequences
    2. Scans subsampled reads for matching k-mers
    3. Counts Y (male mammal) vs W (female bird) marker hits
    4. Infers sex based on marker counts and reports confidence
  - **Status**: Sex marker database building on cluster

---

## Open Tasks
- Test mtDNA Kraken2 database when build completes
- Run comparison script to validate mtDNA vs full database results
- Test sex determination with known samples
- Consider adding species class inference from Kraken2 to improve sex determination

## Notes for Next Session
- mtDNA database will be at: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/kraken/k2_mtdna`
- Sex markers database will be at: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/sex_markers/all_sex_markers.fasta`
- To use new databases:
  ```bash
  nextflow run main.nf \
      --input samplesheet.csv \
      --outdir results \
      --kraken2_db /path/to/k2_mtdna \
      --sex_markers_db /path/to/sex_markers/all_sex_markers.fasta \
      -profile singularity
  ```

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
