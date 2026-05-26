# ==============================================================================
# GWAS survival analysis
# Description: Merges imputed sample data with covariates, recodes variables, 
# handles missing data, and runs Cox proportional hazards models per chromosome.
# Usage: Rscript run_gwasurvivr.R <CHROMOSOME_NUMBER>
# ==============================================================================

library(gwasurvivr)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: Please provide a chromosome number. Example: Rscript run_gwasurvivr.R 22")
}
chr_num <- args[1]

# paths
PATH_SAMPLE    <- "data/cohort_imputed.sample"
PATH_COVARIATE <- "data/cohort_covariates_pheno.txt"
DIR_IMPUTED    <- "data/imputed_genotypes/"
OUT_PREFIX     <- "output/Cox_surv_chr"

# import & merge -----------------------------------------------------
# read files from plink
sample_data <- read.table(PATH_SAMPLE, header = TRUE)
sample_subset <- as.data.frame(sample_data[-1, c(1:2)]) 

covariate_data <- read.table(PATH_COVARIATE, header = TRUE)
covfile <- merge(sample_subset, covariate_data, by = c("ID_1", "ID_2"), all = TRUE)

covfile <- covfile[, c(1, 3:ncol(covfile))]

# recode and filter ------------------------------------------------------
# recode Ssex: Assuming original is 1=Male, 2=Female -> recode to 0=Male, 1=Female
covfile$Sex <- ifelse(covfile$Sex == 1, 0, 
               ifelse(covfile$Sex == 2, 1, -9))

# recode event: Assuming original is 1=Censored, 2=Event -> recode to 0=Censored, 1=Event
covfile$Event <- ifelse(covfile$Event == 1, 0, 
                 ifelse(covfile$Event == 2, 1, NA))

# identify sample IDs with valid TTFR entry
sample.ids <- covfile[!is.na(covfile$TTFR), ]$ID_1

columns_to_fill <- c("TTFR", "AAO", paste0("C", 1:8), "Sex", "Event")

for (col in columns_to_fill) {
  if (col %in% colnames(covfile)) {
    covfile[[col]][is.na(covfile[[col]])] <- -9
  }
}


# gwasurvivr --------------------------------------------------------
cat(sprintf("Running gwasurvivr Cox model for Chromosome %s...\n", chr_num))

impute2CoxSurv(
  impute.file    = paste0(DIR_IMPUTED, "cohort_imputed_chr", chr_num, ".gen.gz"),
  sample.file    = PATH_SAMPLE, 
  covariate.file = covfile,
  id.column      = "ID_1",
  chr            = chr_num,
  sample.ids     = sample.ids,
  time.to.event  = "TTFR",
  event          = "Event",
  covariates     = c(paste0("C", 1:8), "AAO", "Sex"),
  inter.term     = NULL,
  chunk.size     = 10000,
  print.covs     = "only",
  flip.dosage    = TRUE, 
  verbose        = TRUE,
  out.file       = paste0(OUT_PREFIX, chr_num),
  maf.filter     = 0.01,
  clusterObj     = NULL
)

cat(sprintf("Chromosome %s analysis complete.\n", chr_num))
