```markdown
# 🧬 GWAS QC & Downstream Analysis Pipeline

![Bash](https://img.shields.io/badge/Bash-4EAA25?style=for-the-badge&logo=gnu-bash&logoColor=white)
![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Perl](https://img.shields.io/badge/Perl-39457E?style=for-the-badge&logo=perl&logoColor=white)

A computational pipeline for genotype data processing and GWAS analysis. This repository contains the scripts and utilities required to take raw genotype data through quality control, post-imputation processing, and downstream statistical/survival analysis.

## 📂 Repository Structure

The workflow is divided into logical, independent modules:

*   **`utils/`**: Core data-wrangling scripts (e.g., `prep_vcf.sh` for VCF indexing, compression, and tagging, if working with vcf files).
*   **`post_imputation/`**: Scripts for handling imputed data chunks (for example: automated logistic regression for specific chromosomes).
*   **`analysis/`**: Downstream statistical modeling, including GWAS time-to-event/survival analysis.
*   **`outputs/`**: R scripts for generating publication-ready plots (e.g. `enhanced_manhattan.R` for highlighted Manhattan plots).
*   **`metadata/`**: Reference tracking and configuration.


### 1. Environment Setup
To ensure reproducibility, this pipeline relies on Conda. 
```bash
conda env create -f environment.yml
conda activate gwas_env
