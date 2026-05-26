# ==============================================================================
# GWAS summary stats pre-QC & filtering
# Description: Aggregates summary statistics, filters SNPs based on MAF and HWE,
# and generates a consolidated pre-QC dataset alongside an exclusion list.
# Usage: Rscript script_name.R <COHORT_NAME>
# ==============================================================================

# Load libraries
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: Please provide a cohort name. Example: Rscript pre_qc.R COHORTX")
}
cohort <- args[1]

# define paths & params
DIR_STATS   <- "stats"
OUT_EXCLUDE <- paste0(cohort, "_exclude_SNPs.txt")
OUT_RDS     <- paste0(cohort, "_all-stats_preQC.RDS")

THRESH_HWE  <- 1e-6
THRESH_MAF  <- 0.01
# THRESH_INFO <- 0.8

# import -------------------------------------------------
file_list <- list.files(path = DIR_STATS, pattern = "\\.snp-stats\\.gz$", full.names = TRUE)

if (length(file_list) == 0) {
  stop(paste("Error: No .snp-stats.gz files found in", DIR_STATS))
}

cat(sprintf("[%s] Reading and binding %d stats files...\n", cohort, length(file_list)))

stats_list <- lapply(file_list, fread)
stats <- rbindlist(stats_list)

cat("Total SNPs loaded:", nrow(stats), "\n")

# QC filtering -------------------------------------------------
# SNPs failing HWE
exclude_hwe <- stats[HWE < THRESH_HWE]
cat(sprintf("SNPs failing HWE (< %g): %d\n", THRESH_HWE, nrow(exclude_hwe)))

# SNPs failing MAF
exclude_maf <- stats[MAF < THRESH_MAF]
cat(sprintf("SNPs failing MAF (< %g): %d\n", THRESH_MAF, nrow(exclude_maf)))

# (Optional) identify SNPs failing imputation quality
# exclude_info <- stats[impute_info < THRESH_INFO]
# cat(sprintf("SNPs failing INFO (< %g): %d\n", THRESH_INFO, nrow(exclude_info)))

# combine & export
excluded_rsids <- unique(c(exclude_hwe$rsid, exclude_maf$rsid))
cat("Total UNIQUE SNPs flagged for exclusion:", length(excluded_rsids), "\n")

cat("Saving exclusion list to:", OUT_EXCLUDE, "\n")
fwrite(data.table(rsid = excluded_rsids), OUT_EXCLUDE, col.names = FALSE, quote = FALSE)

cat("Saving aggregated pre-QC stats to:", OUT_RDS, "\n")
saveRDS(stats, OUT_RDS)

cat("QC aggregation complete.\n")
