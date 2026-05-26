#!/bin/bash
# ==============================================================================
# Imputation to dosage conversion
# Description: Iterates through chromosome-level IMPUTE2/genotype files (.gen) 
# and converts them to standard dosage text files using a custom Perl script.
# Usage: ./convert_to_dosage.sh
# ==============================================================================

DATA_DIR="data/imputed"
OUT_DIR="output/dosages"
COHORT_PREFIX="cohort_imputed"

PERL_SCRIPT="scripts/to_dosage_from_impute2.pl"

mkdir -p "${OUT_DIR}"

echo "Starting IMPUTE2 to Dosage conversion..."

for chr in {1..22}; do
    IN_FILE="${DATA_DIR}/${COHORT_PREFIX}_chr${chr}.gen"
    OUT_FILE="${OUT_DIR}/${COHORT_PREFIX}_chr${chr}_dosage.txt"

    # Safety check: ensure the input file exists before running
    if [ -f "${IN_FILE}" ]; then
        echo "[$(date +'%H:%M:%S')] Converting Chromosome ${chr} to dosage format..."
        perl "${PERL_SCRIPT}" "${IN_FILE}" "${OUT_FILE}"
    else
        echo "Warning: Input file ${IN_FILE} not found. Skipping Chromosome ${chr}."
    fi
done

echo "Dosage conversion complete!"
