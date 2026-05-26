# ==============================================================================
# Data preparation & survival analysis formatting
# Description: Merges clinical and genetic cohort data, calculates Time to 
# First Relapse (TTFR), processes PCA covariates, and formats files for gwasurvivr.
# ==============================================================================

# Load libraries ---------------------------------------------------------
# set working dir
# setwd("path/to/your/working_directory")

library(readxl)
library(lubridate)
library(ggplot2)
library(gwasurvivr)

# define paths
PATH_SUMMARY_DATA <- "data/summary_data.RDS"
PATH_PATIENT_INFO <- "data/patient_selection_info.xlsx"
PATH_RELAPSE_DATA <- "data/patient_relapse_info.xlsx"
PATH_GENOTYPE_SAMPLE <- "data/genotype_1000G.sample"
PATH_AAO_DATA <- "data/basisdok_aao.txt" ##clinical/insurance information file for age at onset (AAO)
PATH_MDS_COVAR <- "data/qc_cluster.mds"
PATH_HRC_SAMPLE <- "data/imputed_HRC.sample"
OUT_DIR <- "output/"

# import data and merge -----------------------------------------------------
summary_data <- readRDS(PATH_SUMMARY_DATA)
ptlist <- read_xlsx(PATH_PATIENT_INFO)

# map cohort ID to the patient list
ptlist$cohort_id <- summary_data$Cohort_ID[match(ptlist$PatID, summary_data$Master_Pat_ID)]

# QC check: count missing IDs
cat("Missing Cohort IDs:\n")
print(table(is.na(ptlist$cohort_id)))

# export only valid patients
valid_ptlist <- ptlist[!is.na(ptlist$cohort_id), ]
write.table(valid_ptlist[, 6], paste0(OUT_DIR, "patientlist_valid_IDs.txt"), 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# merge relapse data
relapse <- read_xlsx(PATH_RELAPSE_DATA)
ptlist_merged <- merge(ptlist, relapse, by.x = "PatID", by.y = "PATNR")

# filter for patients present in the genotype sample file
geno_sample <- read.table(PATH_GENOTYPE_SAMPLE, header = TRUE)
ptlist_filtered <- ptlist_merged[paste0(ptlist_merged$cohort_id, "_", ptlist_merged$cohort_id) %in% geno_sample$ID_2, ]

# Phenotype processing (TTFR) ------------------------------
ptlist_filtered$TherapyBegin <- as.Date(ptlist_filtered$TherapyBegin)
ptlist_filtered$DateRelapse <- as.Date(ptlist_filtered$DateRelapse)

# add 2-year follow-up threshold
ptlist_filtered$FollowupEnd <- ptlist_filtered$TherapyBegin %m+% years(2)

# merge Age at Onset (AAO)
aao <- read.table(PATH_AAO_DATA, header = TRUE)
ptlist_filtered$AAO <- aao$aao2[match(ptlist_filtered$PatID, aao$PATNR_EXP_ADD)]

# subset relevant clinical columns
clinical_subset <- ptlist_filtered[, c(1, 2, 7, 16, 18, 24, 25)]

# calculate TTFR in months
clinical_subset$TTFR <- ifelse(
  (clinical_subset$DateRelapse >= clinical_subset$TherapyBegin) & (clinical_subset$DateRelapse <= clinical_subset$FollowupEnd),
  interval(ymd(clinical_subset$TherapyBegin), ymd(clinical_subset$DateRelapse)) %/% months(1),
  ifelse(clinical_subset$DateRelapse > clinical_subset$FollowupEnd, 25, NA)
)

# handle cases where TTFR = 0 (relapse within the first month)
ttfr_zero <- clinical_subset[which(clinical_subset$TTFR == 0), ]
ttfr_zero$intervalDate <- interval(ymd(ttfr_zero$TherapyBegin), ymd(ttfr_zero$DateRelapse)) %/% days(1)
ttfr_zero$TTFR2 <- ifelse(ttfr_zero$intervalDate <= 14, 0.5, 1)

clinical_subset[which(clinical_subset$TTFR == 0), ]$TTFR <- 
  ttfr_zero$TTFR2[match(clinical_subset[which(clinical_subset$TTFR == 0), ]$PatID, ttfr_zero$PatID)]

# prep plink .fam parameters
clinical_subset$PID <- 0
clinical_subset$MID <- 0
clinical_subset$sex <- ifelse(clinical_subset$Gender.x == "m", 1, 2)

fam_file <- clinical_subset[, c(3, 3, 9, 10, 11, 8)]
write.table(fam_file[, 1:6], paste0(OUT_DIR, "famfile_formatted.txt"), 
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

fam_file$AAO <- ptlist_filtered$AAO[match(fam_file$cohort_id, ptlist_filtered$cohort_id)]

# process PCA/covariates -----------------------------------------------
mds_data <- read.table(PATH_MDS_COVAR, header = TRUE)
covariates <- mds_data[, c(1, 2, 4:ncol(mds_data))]

# plot relative variance of MDS components
vData <- covariates[, 3:ncol(covariates)]
vData <- data.frame(sapply(vData, var))
colnames(vData) <- "variance"
vData$relVariance <- vData$variance / sum(vData$variance)
vData$PC <- factor(rownames(vData), levels = rownames(vData))

pca_plot <- ggplot(vData, aes(x = PC, y = relVariance, group = 1)) +
  geom_line(colour = "#6baed6", size = 2) +
  geom_point(colour = "#08519c", size = I(5)) +
  xlab("MDS Components") + ylab("Relative variance") + theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey60")) +
  scale_y_continuous(breaks = seq(0, 1, 0.001))

print(pca_plot)

# retain top n components (for example up to C8)
covariates <- covariates[, 1:which(colnames(covariates) == "C8")]

# merge phenotypes with covariates
phenotype <- merge(fam_file, covariates, by.x = c("cohort_id", "cohort_id.1"), by.y = c("FID", "IID"))
phenotype$cohort_id <- paste0(phenotype$cohort_id, "_", phenotype$cohort_id)
phenotype$cohort_id.1 <- paste0(phenotype$cohort_id.1, "_", phenotype$cohort_id.1)

# Format files for gwasurvivr -------------------------------------------
hrc_sample <- read.table(PATH_HRC_SAMPLE, header = TRUE)
hrc_sample <- as.data.frame(hrc_sample[-1, ])

sample_merged <- merge(hrc_sample, phenotype, by.x = c("ID_1", "ID_2"), by.y = c("cohort_id", "cohort_id.1"), all.x = TRUE)
sample_subset <- sample_merged[, c(1, 2, 6:16)] 

# work on pheno file
sample_codes <- as.data.frame(t(c(rep(0, 2), 2, 1, 2, paste0(rep(2, 8)))))
colnames(sample_codes) <- colnames(sample_subset)
pheno_final <- rbind(sample_codes, sample_subset)

# export
write.table(pheno_final, paste0(OUT_DIR, "cohort_phenotype.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

pheno_subset <- pheno_final[which(pheno_final$ID_1 %in% paste0(fam_file$cohort_id, "_", fam_file$cohort_id)), ]
write.table(pheno_subset[, c(1:2, 5:13)], paste0(OUT_DIR, "cohort_covar.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

#export event files
fam_file$cohort_id <- paste0(fam_file$cohort_id, "_", fam_file$cohort_id)
fam_file$cohort_id.1 <- paste0(fam_file$cohort_id.1, "_", fam_file$cohort_id.1)

#binary phenotype: 2 eEvent occurred within 24 months), 1 (censored/no event)
fam_file$phen <- ifelse(fam_file$TTFR <= 24, 2, 1)

fam_merged <- merge(hrc_sample, fam_file, by.x = c("ID_1", "ID_2"), by.y = c("cohort_id", "cohort_id.1"), all.x = TRUE)
colnames(fam_merged)[1:2] <- c("FID", "IID")

write.table(fam_merged[, c(1:2, 4:6, 9)], paste0(OUT_DIR, "cohort_TTFR_phen.fam"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(fam_merged[, c(1:2, 9)], paste0(OUT_DIR, "cohort_TTFR_phen.phen"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(pheno_subset[, c(1:2, 4:13)], paste0(OUT_DIR, "cohort_covar_phen.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# gwasurvivr files ----------------------------------------------------
cov2 <- read.table(paste0(OUT_DIR, "cohort_covar_phen.txt"), header = TRUE)
phen2 <- read.table(paste0(OUT_DIR, "cohort_TTFR_phen.phen"), header = TRUE)
fam2 <- read.table(paste0(OUT_DIR, "cohort_TTFR_phen.fam"))

cov2$Event <- phen2$phen[match(cov2$ID_1, phen2$FID)]
cov2$Sex <- fam2$V5[match(cov2$ID_1, fam2$V1)]

write.table(cov2, paste0(OUT_DIR, "cohort_covar_phen_event.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Example gwasurvivr implementation 
# impute2CoxSurv(
#   impute.file = "path/to/impute", 
#   sample.file = "path/to/sample", 
#   chr = 1, 
#   covariate.file = cov2, 
#   id.column = "ID_1",
#   time.to.event = "TTFR", 
#   event = "Event", 
#   covariates = c("Sex", "C1", "C2", "C3"),
#   out.file = paste0(OUT_DIR, "survival_results")
# )
