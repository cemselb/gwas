#!/bin/bash

#usage: ./prepare_vcf.sh <input.vcf> <output_prefix>

set -euo pipefail # Fail fast on errors

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.vcf> <output_prefix>"
    exit 1
fi

INPUT_VCF=$1
OUTPUT_PREFIX=$2

echo "[INFO] Compressing ${INPUT_VCF}..."
bcftools view "${INPUT_VCF}" -Oz -o "${OUTPUT_PREFIX}.vcf.gz"

echo "[INFO] Indexing ${OUTPUT_PREFIX}.vcf.gz..."
tabix -p vcf "${OUTPUT_PREFIX}.vcf.gz"

echo "[INFO] Filling tags (AN, AC) into ${OUTPUT_PREFIX}_tagged.vcf..."
bcftools +fill-tags "${OUTPUT_PREFIX}.vcf.gz" -Ov -o "${OUTPUT_PREFIX}_tagged.vcf" -- -t AN,AC

echo "[INFO] Processing complete."
