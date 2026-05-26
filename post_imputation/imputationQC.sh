#!/bin/bash
# ==============================================================================
# GWAS imputation QC & summary stats calculation
# Description: Iterates through chromosome-level imputed data (.gen.gz) and 
# calculates SNP statistics using QCTOOL. Includes an optional PLINK dosage run.
# Usage: ./run_qctool_stats.sh
# ==============================================================================

DATA_DIR="data/imputed"
OUT_DIR="output/stats"
COHORT_PREFIX="cohort_imputed"

SAMPLE_FILE="data/cohort.sample"
KEEP_SAMPLES="data/keep_samples.txt" ##list of QC-passing individuals to retain

mkdir -p "${OUT_DIR}"

echo "Starting QCTOOL summary statistics calculation..."

for chr in {1..22}; do
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Processing Chromosome ${chr}..."
    
    # Generate SNP statistics
    qctool -g "${DATA_DIR}/${COHORT_PREFIX}_chr${chr}.gen.gz" \
           -s "${SAMPLE_FILE}" \
           -incl-samples "${KEEP_SAMPLES}" \
           -snp-stats "${OUT_DIR}/${COHORT_PREFIX}_chr${chr}.snp-stats"
           
    # --------------------------------------------------------------------------
    # OPTIONAL DOWNSTREAM: PLINK dosage association testing
    # --------------------------------------------------------------------------
    # FAM_FILE="data/cohort_survival.fam"
    # PHEN_FILE="data/cohort_survival.phen"
    # COVAR_FILE="data/cohort_covar.txt"
    # PHENO_NAME="TTFR"
    # 
    # plink --dosage "${DATA_DIR}/${COHORT_PREFIX}_chr${chr}.gen.gz" noheader skip0=1 skip1=2 format=3 sex \
    #       --fam "${FAM_FILE}" \
    #       --pheno "${PHEN_FILE}" \
    #       --pheno-name "${PHENO_NAME}" \
    #       --covar "${COVAR_FILE}" \
    #       --out "${OUT_DIR}/${COHORT_PREFIX}_chr${chr}"
           
done

echo "All chromosomes processed successfully!"
