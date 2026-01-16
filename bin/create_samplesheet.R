#!/usr/bin/env Rscript
# Script to create a samplesheet from CNAG XLS file
# Usage: Rscript create_samplesheet.R <project.xls> <fastq_path> <output.csv>

library(readxl)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript create_samplesheet.R <project.xls> <fastq_path> <output.csv>")
}

xls_path <- args[1]
fastq_path <- args[2]
output_csv <- args[3]

# Read the Excel file (skip first 2 rows as in original script)
xls_data <- read_excel(xls_path, skip = 2)

# Create samplesheet dataframe
samplesheet <- data.frame(
  sample = xls_data$FLI,
  fastq_1 = paste0(fastq_path, "/", xls_data$FLI, "_1.fastq.gz"),
  fastq_2 = paste0(fastq_path, "/", xls_data$FLI, "_2.fastq.gz")
)

# Write output CSV
write.csv(samplesheet, output_csv, row.names = FALSE, quote = FALSE)

cat("Samplesheet created:", output_csv, "\n")
cat("Number of samples:", nrow(samplesheet), "\n")
