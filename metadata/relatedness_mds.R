# ==============================================================================
# GWAS relatedness QC
# Description: Parses PLINK IBD (.genome) and missingness (.imiss) files to 
# identify related individuals (cryptic relatedness) and generates a removal 
# list prioritising the retention of samples with lower missingness.
# ==============================================================================

# set working dir -----------------------------------------------------
# setwd("path/to/working_directory/mds/")

# define paths
PATH_IBD_GENOME    <- "data/plink.genome"
PATH_MISSINGNESS   <- "data/plink.imiss"
PATH_MDS_POST_QC   <- "data/pca_norel_cluster.mds"
OUT_REMOVE_LIST    <- "output/remove_relatives.txt"

# params
# 0.125 = 3rd degree relatives (1st cousins). Can adjust to 0.0625 for 4th degree.
PI_HAT_THRESHOLD   <- 0.125 

# load IBD file ----------------------------------------------
genome <- read.table(PATH_IBD_GENOME, header = TRUE)

# subset relatives based on PI_HAT threshold (extracting FID1, IID1, FID2, IID2, PI_HAT)
relatives <- genome[which(genome$PI_HAT >= PI_HAT_THRESHOLD), c(1, 2, 3, 4, 10)]
relatives <- relatives[order(-relatives$PI_HAT), ]

cat("Total related pairs found (PI_HAT >=", PI_HAT_THRESHOLD, "):", nrow(relatives), "\n")

# Optional: check how many related individuals belong to your specific target cohort
# Assuming 'target_sample_ids' exists in your environment from previous filtering
# cat("Relatives in target cohort (FID1):", 
#     nrow(relatives[which(paste0(relatives$FID1, "_", relatives$FID1) %in% target_sample_ids), ]), "\n")
# cat("Relatives in target cohort (FID2):", 
#     nrow(relatives[which(paste0(relatives$FID2, "_", relatives$FID2) %in% target_sample_ids), ]), "\n")

# missingness and generate removal list -------------------------------
missingness <- read.table(PATH_MISSINGNESS, header = TRUE)
rel_miss <- merge(relatives, missingness[, c(1:2, 6)], 
                  by.x = c("FID1", "IID1"), by.y = c("FID", "IID"))

rel_miss <- merge(rel_miss, missingness[, c(1:2, 6)], 
                  by.x = c("FID2", "IID2"), by.y = c("FID", "IID"), 
                  suffix = c("1", "2"), all.x = TRUE)

# identify which individual in each pair has the higher missingness rate
drop_indiv1 <- as.character(rel_miss[which(rel_miss$F_MISS1 >= rel_miss$F_MISS2), ]$IID1)
drop_indiv2 <- as.character(rel_miss[which(rel_miss$F_MISS2 > rel_miss$F_MISS1), ]$IID2)

# combine and format for PLINK removal list (FID and IID)
samples_to_remove <- unique(c(drop_indiv1, drop_indiv2))
samples_to_remove_df <- data.frame(FID = samples_to_remove, IID = samples_to_remove)

cat("Total unique individuals marked for removal:", nrow(samples_to_remove_df), "\n")

# export list for plink (--remove)
write.table(samples_to_remove_df, OUT_REMOVE_LIST, 
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

# After PLINK removal - load postQC covars -----------------------------
# assuming you ran PLINK here with the generated 'remove_relatives.txt' to 
# generate a new MDS/PCA file without related individuals.

mds_post_qc <- read.table(PATH_MDS_POST_QC, header = TRUE)
covar <- mds_post_qc[, c(1, 2, 4:ncol(mds_post_qc))]

head(covar)
