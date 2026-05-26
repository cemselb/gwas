# ==============================================================================
# GWAS survival results aggregation & visualisation
# Description: aggregates chromosome-level Cox proportional hazards results (.coxph),
# merges with pre-QC summary statistics, checks against excluded SNPs, and
# generates Manhattan plots for various filtering thresholds.
# ==============================================================================

# load libraries ---------------------------------------------------------
# setwd("path/to/working_directory/")

library(data.table)
library(qqman)

# paths and params
DIR_COX_RESULTS  <- "results/coxph/"
PATH_EXCLUDE     <- "data/cohort_exclude_SNPs.txt"
PATH_PREQC_STATS <- "data/cohort_all-stats_preQC.RDS"
OUT_RESULTS      <- "output/Cox_TTFR_all_chrs_results.txt"

MANHATTAN_COLORS <- c("#8060E0", "#C69ADD", "#80E372", "#D1DAA2", "#DC52D3", 
                      "#D5DA54", "#D7D1D7", "#76AECE", "#D59A67", "#DC6788", "#7DE0C8")

# aggregate (chr 1-22) -----------------------------------------------
cat("Aggregating chromosome-level results...\n")

cox_files <- paste0(DIR_COX_RESULTS, "Cox_gen_chr", 1:22, ".coxph")
results_list <- lapply(cox_files, function(f) fread(f, stringsAsFactors = FALSE, header = TRUE))
results <- rbindlist(results_list)

cat("Total SNPs loaded:", nrow(results), "\n")

# QC ----------------------------------------------------
# check for missing p-values
cat("SNPs with missing P-values:", nrow(results[is.na(results$PVALUE), ]), "\n")

# check overlap with previously excluded SNPs
excl <- fread(PATH_EXCLUDE, header = FALSE)
cat("Excluded SNPs present in results:", nrow(results[results$RSID %in% excl$V1, ]), "\n")

# Optional: quick initial Manhattan plot for P < 0.01
# manhattan(results[results$PVALUE < 1e-2, ], bp = "POS", chr = "CHR", snp = "RSID", p = "PVALUE")

# merge with preQC stats ----------------------------------------------
cat("Loading pre-QC stats and merging...\n")
stats <- readRDS(PATH_PREQC_STATS)

stats_subset <- stats[, c("rsid", "minor_allele_frequency", "SAMP_MAF", "HW_exact_p_value", "info")]

results_merged <- merge(as.data.frame(results), stats_subset, all.x = TRUE, by.x = "RSID", by.y = "rsid")

cat("SNPs with global MAF < 0.01:", nrow(results_merged[which(results_merged$minor_allele_frequency < 0.01), ]), "\n")
cat("SNPs with sample MAF < 0.01:", nrow(results_merged[which(results_merged$SAMP_MAF < 0.01), ]), "\n")
cat("SNPs failing HWE (p < 1e-6):", nrow(results_merged[which(results_merged$HW_exact_p_value < 1e-6), ]), "\n")

# export
cat("Saving aggregated results to:", OUT_RESULTS, "\n")
fwrite(results_merged, OUT_RESULTS, sep = "\t", quote = FALSE, row.names = FALSE)

#plots -----------------------------------------------------------
cat("Generating Manhattan plots...\n")

results_plot <- results_merged[-which(results_merged$PVALUE > 1.5e-2), ]

# plot 1: unfiltered (but restricted to p <= 0.015)
manhattan(results_plot, 
          chr = "CHR", snp = "RSID", bp = "POS", p = "PVALUE",
          main = "GWAS Survival Analysis (TTFR) - Unfiltered", 
          col = MANHATTAN_COLORS, ylim = c(1, 10), cex = 0.5, 
          chrlabs = as.character(c(1:22)), cex.axis = 0.8)

# plot 2: filtered by INFO > 0.8 and sample MAF > 0.05
manhattan(results_plot[which(results_plot$info > 0.8 & results_plot$SAMP_MAF > 0.05), ], 
          chr = "CHR", snp = "RSID", bp = "POS", p = "PVALUE",
          main = "GWAS Survival Analysis (TTFR) - INFO > 0.8 & MAF > 0.05", 
          col = MANHATTAN_COLORS, ylim = c(1, 10), cex = 0.5, 
          chrlabs = as.character(c(1:22)), cex.axis = 0.8)

# plot 3: filtered by sample MAF > 0.05 only
manhattan(results_plot[which(results_plot$SAMP_MAF > 0.05), ], 
          chr = "CHR", snp = "RSID", bp = "POS", p = "PVALUE",
          main = "GWAS Survival Analysis (TTFR) - MAF > 0.05", 
          col = MANHATTAN_COLORS, ylim = c(1, 10), cex = 0.5, 
          chrlabs = as.character(c(1:22)), cex.axis = 0.8)
