# GWAS Analysis of Inflammatory Bowel Disease (IBD)

A genome-wide association study analyzing genetic variants associated with Inflammatory Bowel Disease, developed for the "Statistical Models in Human Genetics" module in the Genomic Sciences program at UNAM Juriquilla.

## Overview

This project implements a complete GWAS pipeline including:
- Data preprocessing and quality control
- Sample and SNP filtering
- Population stratification correction via PCA
- Association testing using logistic regression
- Results visualization (Manhattan and QQ plots)

## Data

The analysis uses the PennCath cohort dataset:
- **Genotype data**: ~500,000 SNPs across 1,401 individuals (PLINK format: `.bed`, `.bim`, `.fam`)
- **Clinical data**: `Clinical_IBD.csv` containing phenotype (IBD case/control status), age, sex, and lipid measurements

## Pipeline

| Script | Description |
|--------|-------------|
| `MyGWAS_Lab1.R` | Data import, initial SNP QC (call rate ≥95%, MAF >1%) |
| `MyGWAS_Lab2.R` | Sample QC: heterozygosity, call rate, HWE, IBD/relatedness filtering, PCA |
| `MyGWAS_Lab3.R` | GWAS execution using parallelized GLM |
| `MyGWAS_LabFunctions.R` | Visualization functions (Manhattan plot, QQ plot, λGC calculation) |
| `GWAS_IBD_RodriguezAxel.R` | Complete integrated pipeline |

## Quality Control Thresholds

- SNP call rate: ≥95%
- Minor allele frequency: >1%
- Sample call rate: ≥90%
- Inbreeding coefficient \|F\|: ≤0.055
- Hardy-Weinberg equilibrium: p > 10⁻⁶ (controls only)
- Kinship coefficient: <0.05

## Dependencies
```r
# Bioconductor
BiocManager::install(c("snpStats", "SNPRelate"))

# CRAN
install.packages(c("dplyr", "doParallel", "qqman", "gdsfmt"))
```

## Usage

Run scripts sequentially (Lab1 → Lab2 → Lab3) or use the integrated script:
```r
source("GWAS_IBD_RodriguezAxel.R")
```

**Note**: Update file paths in scripts to match your local directory structure.

## Output

- `GWAA.txt`: Association results (SNP, estimate, SE, t-value, p-value)
- Manhattan and QQ plots for visualization
- Genomic inflation factor (λGC) for assessing population stratification

## Author

Axel Rodriguez Perez  
Genomic Sciences, UNAM · November 2023
