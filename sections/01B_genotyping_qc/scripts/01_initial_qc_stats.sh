#!/bin/bash

################################################################################
# Section 1B: Genotyping QC — Step 01: Initial Quality Assessment
# 
# PURPOSE:
#   Compute initial quality statistics on the raw genotype data to understand
#   the baseline data quality and identify which variants/samples need filtering.
#   This is a purely observational step — we don't filter yet, just measure.
#
# INPUT:
#   - pdac_demo.bed/bim/fam (raw genotypes)
#
# OUTPUT:
#   - pdac_demo.afreq        (allele frequencies)
#   - pdac_demo.het          (individual heterozygosity)
#   - pdac_demo.missingxy    (missing data by chr, sex)
#   - pdac_demo.nobs         (number of observations per variant)
#
# NOTES:
#   - GRCh38 genome build (assumed in the bim file)
#   - PLINK 2 syntax (plink2)
#   - No filtering applied yet
#
################################################################################

set -e  # Exit immediately on error

# Configuration
SEED="${1:-2026}"
DATASET_DIR="../../demo_dataset/data"
DATASET_NAME="pdac_demo"

# ============================================================================
# STEP 1: Compute allele frequency spectrum
# ============================================================================
echo ""
echo "=== Computing allele frequencies ==="
echo ""

# Why --freq?
#   Allele frequency is essential for:
#   - Understanding the allele distribution (rare vs common variants)
#   - Computing Hardy-Weinberg p-values later (if needed)
#   - Deciding minimum allele frequency filters
#
# Why --ref-allele?
#   By default, PLINK2 assumes the last allele in .bim is the reference.
#   This is already set correctly in our prepared dataset.

plink2 \
  --bfile ${DATASET_DIR}/${DATASET_NAME} \
  --freq counts \
  --out ${DATASET_NAME}_qc

echo "✓ Allele frequencies saved to ${DATASET_NAME}_qc.afreq"

# ============================================================================
# STEP 2: Compute individual heterozygosity and inbreeding coefficient
# ============================================================================
echo ""
echo "=== Computing heterozygosity per individual ==="
echo ""

# Why --het?
#   Heterozygosity is used to identify:
#   - Potential sample contamination (unusually high heterozygosity)
#   - Potential inbreeding (unusually low heterozygosity)
#   - Sample swaps or quality issues
#
# The F coefficient (inbreeding) ranges from -1 (all heterozygous) to +1 (all homozygous).
# We'll use this in Step 4 to filter outliers.

plink2 \
  --bfile ${DATASET_DIR}/${DATASET_NAME} \
  --het \
  --out ${DATASET_NAME}_qc

echo "✓ Heterozygosity computed, saved to ${DATASET_NAME}_qc.het"

# ============================================================================
# STEP 3: Compute missing data rate by individual and variant
# ============================================================================
echo ""
echo "=== Computing missing data rates ==="
echo ""

# Why missing data rates?
#   High missingness can indicate:
#   - Poor genotyping quality
#   - PCR failures or technical issues
#   - Low-quality samples
#
# We'll use the per-sample and per-variant missingness to apply filters in Step 2.

plink2 \
  --bfile ${DATASET_DIR}/${DATASET_NAME} \
  --missing \
  --out ${DATASET_NAME}_qc

echo "✓ Missing data rates saved to:"
echo "  - ${DATASET_NAME}_qc.imiss (per individual)"
echo "  - ${DATASET_NAME}_qc.lmiss (per variant/locus)"

# ============================================================================
# STEP 4: Inspect basic summary statistics
# ============================================================================
echo ""
echo "=== Computing summary statistics ==="
echo ""

# Variant count
NVAR=$(awk 'NR==2 {print $5}' ${DATASET_NAME}_qc.afreq)

# Sample count
NSAMP=$(tail -n +2 ${DATASET_NAME}_qc.imiss | wc -l)

# Median heterozygosity
MEDIAN_HET=$(tail -n +2 ${DATASET_NAME}_qc.het | awk '{print $(NF-1)}' | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)]}')

echo "Baseline statistics:"
echo "  - Variants:         ${NVAR}"
echo "  - Individuals:      ${NSAMP}"
echo "  - Median F (inbreeding): ${MEDIAN_HET}"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the variant filtering step:"
echo ""
echo "  bash 02_variant_filtering.sh"
echo ""
echo "This will apply filters for:"
echo "  - Variant call rate (--geno)"
echo "  - Hardy-Weinberg equilibrium (--hwe)"
echo "  - Minor allele frequency (--maf) [optional]"
echo ""
