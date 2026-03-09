#!/usr/bin/env nextflow

/*
========================================================================================
    BasicQC Pipeline
========================================================================================
    A Nextflow pipeline for basic QC of Illumina FASTQ files
    - FastQC: Quality control metrics
    - FastQ Screen: Species/contamination detection
    - Kraken2: Taxonomic classification
    - MultiQC: Aggregate reports
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Import modules
include { FASTQC                                          } from './modules/fastqc'
include { FASTQ_SCREEN                                    } from './modules/fastq_screen'
include { SEQTK_SUBSAMPLE                                 } from './modules/seqtk_subsample'
include { SEQTK_SUBSAMPLE as SEQTK_SUBSAMPLE_RRNA         } from './modules/seqtk_subsample'
include { KRAKEN2                                         } from './modules/kraken2'
include { KRAKEN2 as KRAKEN2_RRNA                         } from './modules/kraken2'
include { SUMMARIZE_KRAKEN2                               } from './modules/summarize_kraken2'
include { SUMMARIZE_RRNA_KRAKEN2                          } from './modules/summarize_kraken2'
include { SEX_DETERMINATION                               } from './modules/sex_determination'
include { SUMMARIZE_SEX                                   } from './modules/sex_determination'
include { SORTMERNA_INDEX                                 } from './modules/sortmerna'
include { SORTMERNA                                       } from './modules/sortmerna'
include { RIBODETECTOR                                    } from './modules/ribodetector'
include { SUMMARIZE_RRNA                                  } from './modules/summarize_rrna'
include { MULTIQC                                         } from './modules/multiqc'
include { PREPARE_MULTIQC_CONFIG                          } from './modules/prepare_multiqc_config'
include { SUMMARIZE_RESULTS                               } from './modules/summarize_results'

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================
     BasicQC Pipeline v1.0
    =========================================
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results

    Mandatory arguments:
      --input           Path to input samplesheet (CSV format)
      --outdir          Output directory for results

    Optional arguments:
      --fastq_screen_conf   Path to fastq_screen configuration file
      --kraken2_db          Path to Kraken2 database
      --subsample_reads     Reads to subsample per sample for all tools (default: 1000000)
      --kraken2_subsample   Override subsample depth for Kraken2 + sex determination only
      --rrna_subsample      Override subsample depth for SortMeRNA + RiboDetector only
      --sex_markers_db      Path to sex marker FASTA for sex determination
      --sortmerna_db        Path to directory containing rRNA FASTA database files
      --rrna_kraken2_db     Path to Kraken2 database for rRNA-based species ID
      --read_length         Sequenced read length in bp for RiboDetector (default: 150)
      --skip_fastqc         Skip FastQC step
      --skip_fastq_screen   Skip FastQ Screen step
      --skip_kraken2        Skip Kraken2 step
      --skip_sex_determination  Skip sex determination step
      --skip_sortmerna      Skip SortMeRNA rRNA quantification step
      --project_name        Project name for MultiQC report header (e.g., 'CGLZOO_01')
      --application         Application type for MultiQC header (e.g., 'RNA-seq')
      -profile              Configuration profile (singularity, docker, conda)

    Samplesheet format (CSV):
      sample,fastq_1,fastq_2,sample_name,species
      HFYMJDSXC_1_8bp-UDP0032,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,BB1523,Callithrix geoffroyi
      HFYMJDSXC_1_8bp-UDP0034,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,BB1525,Gorilla gorilla

    Note: sample_name and species columns are optional but enable better MultiQC grouping
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check mandatory parameters
if (!params.input) {
    error "Please provide an input samplesheet with --input"
}

if (!params.outdir) {
    error "Please provide an output directory with --outdir"
}

// Check input file exists
input_file = file(params.input)
if (!input_file.exists()) {
    error "Input samplesheet not found: ${params.input}"
}

/*
========================================================================================
    INPUT CHANNEL
========================================================================================
*/

def parse_samplesheet(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample = row.sample
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = row.fastq_2 ? file(row.fastq_2) : null

            if (!fastq_1.exists()) {
                error "FASTQ file not found: ${row.fastq_1}"
            }
            if (fastq_2 && !fastq_2.exists()) {
                error "FASTQ file not found: ${row.fastq_2}"
            }

            return fastq_2 ? tuple(sample, [fastq_1, fastq_2]) : tuple(sample, [fastq_1])
        }
}

// Parse samplesheet for metadata (sample_name, species)
def parse_samplesheet_metadata(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            [
                fli: row.sample,
                sample_name: row.sample_name ?: row.sample,
                species: row.species ?: ''
            ]
        }
        .collect()
}

// Parse samplesheet per-sample: groups by sample_name, takes the first FLI's reads.
// Used for all tools that run once per sample (Kraken2, sex det, SortMeRNA, RiboDetector).
def parse_samplesheet_per_sample(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_name = row.sample_name ?: row.sample
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = row.fastq_2 ? file(row.fastq_2) : null
            def species = row.species ?: ''

            return fastq_2 ?
                tuple(sample_name, species, [fastq_1, fastq_2]) :
                tuple(sample_name, species, [fastq_1])
        }
        .groupTuple(by: [0, 1])  // Group by sample_name and species
        .map { sample_name, species, reads_list ->
            // Take only the first FASTQ pair for this sample
            tuple(sample_name, species, reads_list[0])
        }
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Parse samplesheet and create input channel
    ch_reads = parse_samplesheet(params.input)

    // Parse sample metadata for MultiQC config
    ch_sample_metadata = parse_samplesheet_metadata(params.input)

    // Initialize empty channels for MultiQC
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { it[1] })
    }

    //
    // MODULE: FastQ Screen
    //
    if (!params.skip_fastq_screen) {
        // Check for fastq_screen config
        if (!params.fastq_screen_conf) {
            log.warn "No FastQ Screen config provided (--fastq_screen_conf). Skipping FastQ Screen."
        } else {
            ch_fastq_screen_conf = file(params.fastq_screen_conf)
            FASTQ_SCREEN(ch_reads, ch_fastq_screen_conf)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_SCREEN.out.txt.map { it[1] })
        }
    }

    //
    // Per-sample subsampling — shared by all downstream tools
    // (Kraken2, sex determination, SortMeRNA, RiboDetector)
    //
    // All tools run once per sample_name (first FLI when multiple exist).
    // FastQC and FastQ Screen run per-FLI above; everything else uses this subsample.
    //
    // Default: one subsample at subsample_reads (1M) shared by all tools.
    // Override: set --kraken2_subsample and --rrna_subsample to different values
    //           to get separate subsampled sets for each group.
    //
    def run_sortmerna    = !params.skip_sortmerna && params.sortmerna_db
    def run_ribodetector = !params.skip_ribodetector
    def run_rrna_kraken2 = !params.skip_rrna_kraken2 && params.rrna_kraken2_db && run_ribodetector

    def k2_n   = (params.kraken2_subsample != null) ? params.kraken2_subsample : params.subsample_reads
    def rrna_n = (params.rrna_subsample    != null) ? params.rrna_subsample    : params.subsample_reads

    def needs_subsample = (!params.skip_kraken2 && params.kraken2_db) ||
                          (!params.skip_sex_determination && params.sex_markers_db) ||
                          run_sortmerna || run_ribodetector

    if (needs_subsample) {
        // One FASTQ per sample_name (first FLI): same approach for all tools
        ch_per_sample_reads = parse_samplesheet_per_sample(params.input)
            .map { sample_name, species, reads -> tuple(sample_name, reads) }

        SEQTK_SUBSAMPLE(ch_per_sample_reads, k2_n)
        ch_subsampled = SEQTK_SUBSAMPLE.out.reads

        // rRNA reads: same subsample unless explicitly configured otherwise
        ch_rrna_reads = ch_subsampled
        if (k2_n != rrna_n) {
            SEQTK_SUBSAMPLE_RRNA(ch_per_sample_reads, rrna_n)
            ch_rrna_reads = SEQTK_SUBSAMPLE_RRNA.out.reads
        }
    }

    //
    // Shared per-sample metadata (sample_name + species) — used by Kraken2 and sex det
    //
    ch_per_sample_meta = parse_samplesheet_per_sample(params.input)
        .map { sample_name, species, reads -> [sample_name: sample_name, species: species] }
        .collect()

    //
    // MODULE: Kraken2 mtDNA — one sample per sample_name, subsampled reads
    //
    if (!params.skip_kraken2) {
        if (!params.kraken2_db) {
            log.warn "No Kraken2 database provided (--kraken2_db). Skipping Kraken2."
        } else {
            ch_kraken2_db = file(params.kraken2_db)
            KRAKEN2(ch_subsampled, ch_kraken2_db)

            ch_kraken2_reports = KRAKEN2.out.report.map { it[1] }.collect()
            SUMMARIZE_KRAKEN2(ch_kraken2_reports, ch_per_sample_meta)
            ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_KRAKEN2.out.summary)
            ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_KRAKEN2.out.classified_reports.flatten())
        }
    }

    //
    // MODULE: Sex Determination — reuses ch_subsampled (same reads as Kraken2)
    //
    if (!params.skip_sex_determination) {
        if (!params.sex_markers_db) {
            log.warn "No sex markers database provided (--sex_markers_db). Skipping sex determination."
        } else {
            ch_sex_markers_db = file(params.sex_markers_db)
            SEX_DETERMINATION(ch_subsampled, ch_sex_markers_db, 'unknown')

            ch_sex_results = SEX_DETERMINATION.out.results.map { it[1] }.collect()
            SUMMARIZE_SEX(ch_sex_results, ch_per_sample_meta)
            ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_SEX.out.summary)
        }
    }

    //
    // MODULE: SortMeRNA — per-sample subsampled reads (ch_rrna_reads)
    // Index is built once per run by SORTMERNA_INDEX and shared across all samples.
    //
    if (run_sortmerna) {
        ch_sortmerna_fastas = file(params.sortmerna_db).isDirectory()
            ? Channel.fromPath("${params.sortmerna_db}/*.{fasta,fa,fna}").collect()
            : Channel.of(file(params.sortmerna_db)).collect()

        SORTMERNA_INDEX(ch_sortmerna_fastas)
        ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()

        // SortMeRNA provides % rRNA metric only; reads no longer passed to Kraken2
        SORTMERNA(ch_rrna_reads, ch_sortmerna_fastas, ch_sortmerna_index, 'false')
    } else if (!params.skip_sortmerna) {
        log.warn "No SortMeRNA database provided (--sortmerna_db). Skipping SortMeRNA."
    }

    //
    // MODULE: RiboDetector — per-sample subsampled reads (ch_rrna_reads)
    // Provides % rRNA metric and, when rrna_kraken2_db is set, saves rRNA reads
    // for downstream Kraken2 species ID (fewer mRNA false positives than SortMeRNA).
    //
    if (run_ribodetector) {
        RIBODETECTOR(ch_rrna_reads, params.read_length, run_rrna_kraken2.toString())

        if (run_rrna_kraken2) {
            ch_rrna_kraken2_db = file(params.rrna_kraken2_db)
            KRAKEN2_RRNA(
                RIBODETECTOR.out.rrna_reads.map { sample, reads -> tuple("${sample}_rrna", reads) },
                ch_rrna_kraken2_db
            )
            SUMMARIZE_RRNA_KRAKEN2(
                KRAKEN2_RRNA.out.report.map { it[1] }.collect(),
                ch_sample_metadata
            )
            ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_RRNA_KRAKEN2.out.summary)
        }
    }

    //
    // MODULE: Summarize rRNA quantification for MultiQC (proper per-sample naming)
    //
    if (run_sortmerna || run_ribodetector) {
        ch_smr_logs  = run_sortmerna    ? SORTMERNA.out.log.map    { it[1] }.collect() : Channel.of(file("NO_SORTMERNA"))
        ch_ribo_logs = run_ribodetector ? RIBODETECTOR.out.log.map { it[1] }.collect() : Channel.of(file("NO_RIBODETECTOR"))
        SUMMARIZE_RRNA(ch_smr_logs, ch_ribo_logs, ch_sample_metadata)
        ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_RRNA.out.summary)
    }

    //
    // Generate MultiQC config with sample metadata
    //
    PREPARE_MULTIQC_CONFIG(
        ch_sample_metadata,
        params.project_name ?: '',
        params.application ?: ''
    )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files
        .flatten()
        .collect()
        .filter { it.size() > 0 }
        .set { ch_multiqc_input }

    MULTIQC(
        ch_multiqc_input,
        PREPARE_MULTIQC_CONFIG.out.config
    )

    //
    // MODULE: Generate consolidated summary table
    //
    // Collect FastQC zip files
    ch_fastqc_for_summary = params.skip_fastqc
        ? Channel.of(file("NO_FASTQC"))
        : FASTQC.out.zip.map { it[1] }.collect()

    // Get Kraken2 summary (or placeholder)
    ch_kraken2_for_summary = (params.skip_kraken2 || !params.kraken2_db)
        ? Channel.of(file("NO_KRAKEN2"))
        : SUMMARIZE_KRAKEN2.out.summary

    // Get sex determination summary (or placeholder)
    ch_sex_for_summary = (params.skip_sex_determination || !params.sex_markers_db)
        ? Channel.of(file("NO_SEX"))
        : SUMMARIZE_SEX.out.summary

    // Get SortMeRNA logs (or placeholder)
    ch_sortmerna_for_summary = (params.skip_sortmerna || !params.sortmerna_db)
        ? Channel.of(file("NO_SORTMERNA"))
        : SORTMERNA.out.log.map { it[1] }.collect()

    // Get RiboDetector logs (or placeholder)
    ch_ribodetector_for_summary = params.skip_ribodetector
        ? Channel.of(file("NO_RIBODETECTOR"))
        : RIBODETECTOR.out.log.map { it[1] }.collect()

    // Get rRNA Kraken2 reports (or placeholder)
    ch_rrna_kraken2_for_summary = run_rrna_kraken2
        ? KRAKEN2_RRNA.out.report.map { it[1] }.collect()
        : Channel.of(file("NO_RRNA_KRAKEN2"))

    // Parse sample info for summary
    ch_summary_sample_info = ch_sample_metadata.collect()

    SUMMARIZE_RESULTS(
        ch_fastqc_for_summary,
        ch_kraken2_for_summary,
        ch_sex_for_summary,
        ch_sortmerna_for_summary,
        ch_ribodetector_for_summary,
        ch_rrna_kraken2_for_summary,
        ch_summary_sample_info
    )
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'failed'}"
    log.info "Results saved to: ${params.outdir}"
    log.info ""
    log.info "Key outputs:"
    log.info "  - Summary table: ${params.outdir}/summary/qc_summary.tsv"
    log.info "  - MultiQC report: ${params.outdir}/multiqc/"
    log.info ""
}
