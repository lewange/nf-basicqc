#!/bin/bash
#SBATCH --job-name=nf_basicqc
#SBATCH --output=nf_basicqc_%j.out
#SBATCH --error=nf_basicqc_%j.err
#SBATCH --partition=genD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --qos=marathon
#SBATCH --time=48:00:00

# BasicQC Nextflow Pipeline - Production Run
#
# Usage:
#   sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name]
#
# Example:
#   sbatch submit_pipeline.sh inputs/CGLZOO_01.csv results/CGLZOO_01 CGLZOO_01

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name]"
    exit 1
fi

SAMPLESHEET="$1"
OUTDIR="$2"
PROJECT_NAME="${3:-basicqc}"

# Set paths
PIPELINE_DIR="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/BasicQC/nf-basicqc"
FASTQ_SCREEN_CONF="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/FastQ_Screen_Genomes/fastq_screen.conf"
KRAKEN2_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/kraken"

cd $PIPELINE_DIR

echo "$(date) Starting BasicQC pipeline"
echo "=================================="
echo "Samplesheet: $SAMPLESHEET"
echo "Output dir:  $OUTDIR"
echo "Project:     $PROJECT_NAME"
echo ""

# Load modules if needed (uncomment/modify as needed)
# module load nextflow
# module load singularity

# Run the pipeline
nextflow run main.nf \
    --input $SAMPLESHEET \
    --outdir $OUTDIR \
    --fastq_screen_conf $FASTQ_SCREEN_CONF \
    --kraken2_db $KRAKEN2_DB \
    --kraken2_subsample 5000000 \
    --multiqc_title "$PROJECT_NAME" \
    -profile singularity,slurm \
    -resume

echo "$(date) Pipeline complete"
echo "Results in: $OUTDIR"
