#!/bin/bash

################################################################################
# Section 1B: Genotyping QC — Step 07: Relatedness & Sample Pruning (IBD/PI_HAT)
# 
# PURPOSE:
#   Identify and remove related samples. Relatedness violates the assumption
#   of independent samples in GWAS, leading to inflated test statistics and
#   false positive associations.
#
# INPUT:
#   - pdac_demo_06_filt.bed/bim/fam (from Step 06, HWE-filtered)
#
# OUTPUT:
#   - pdac_demo_07_filt.bed/bim/fam — pruned to unrelated samples
#   - pdac_demo_07_kinship.king.gz — kinship coefficients (binary)
#   - pdac_demo_07_removal_list.txt — related pairs identified
#
# THRESHOLD:
#   PI_HAT > 0.1875  (2nd-degree relatives: parent-child, sibling)
#   Note: PI_HAT is the proportion of genome identical-by-descent (IBD)
#
# METHOD:
#   1. LD-prune variants (keep ~30K for kinship inference; too many cause noise)
#   2. Compute kinship coefficients (IBD) using KING algorithm
#   3. Identify related pairs (PI_HAT > 0.1875)
#   4. Remove lower-quality sample from each related pair
#
# NOTES:
#   - PDAC may have related samples (family cohorts)
#   - Removing one per related pair preserves power while maintaining independence
#   - KING algorithm is fast and accurate for unrelated to 3rd-degree relatives
#
################################################################################

set -e

# Configuration
SEED="${1:-2026}"
DATASET_INPUT="${2:-.}"
DATASET_NAME="pdac_demo"

# ============================================================================
# STEP 1: LD-based SNP pruning
# ============================================================================
echo ""
echo "=== LD-pruning variants for kinship inference ==="
echo ""

# Why prune SNPs for kinship?
#   Too many variants make kinship estimates noisy (redundant information).
#   Standard practice: keep ~30K-50K independent variants.
#   Parameters:
#   --indep-pairwise <window_size> <window_step> <r2_threshold>
#   - 50 bp window (50 SNPs at a time)
#   - 5 SNP step (move by 5 SNPs)
#   - r2 > 0.2 (remove if more correlated than this)

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_06_filt \
  --indep-pairwise 50 5 0.2 \
  --out ${DATASET_NAME}_07_prune

echo "✓ LD-pruning complete. Kept $(wc -l < ${DATASET_NAME}_07_prune.prune.in) independent variants"

# ============================================================================
# STEP 2: Compute kinship coefficients
# ============================================================================
echo ""
echo "=== Computing kinship coefficients (KING algorithm) ==="
echo ""

# Why KING algorithm?
#   KING (Kinship-based INference for Gwas) is robust to population stratification.
#   Computes PI_HAT (proportion of genome IBD):
#   - Unrelated: PI_HAT ≈ 0
#   - 2nd degree (siblings, parent-child): PI_HAT ≈ 0.25
#   - 1st degree (parent-child, sibling): PI_HAT ≈ 0.5
#   - Threshold 0.1875 catches 2nd-degree relatives
#
# Binary KING output (.king.gz) is memory-efficient for large cohorts.

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_06_filt \
  --extract ${DATASET_NAME}_07_prune.prune.in \
  --king-cutoff 0.1875 \
  --out ${DATASET_NAME}_07_kinship

echo "✓ Kinship analysis complete"

# ============================================================================
# STEP 3: Extract removal list
# ============================================================================
echo ""
echo "=== Identifying samples to remove ==="
echo ""

# PLINK2 --king-cutoff automatically creates a list of samples to keep
# We need to invert this to get a removal list

# Get samples to KEEP
if [ -f ${DATASET_NAME}_07_kinship.king.cutoff.in.id ]; then
  awk '{print $1, $2}' ${DATASET_NAME}_07_kinship.king.cutoff.in.id > ${DATASET_NAME}_07_keep_samples.txt
else
  echo "No king.cutoff file found; using all samples"
  awk '{print $1, $2}' ${DATASET_INPUT}/${DATASET_NAME}_06_filt.fam > ${DATASET_NAME}_07_keep_samples.txt
fi

# Get all samples and find those NOT in the keep list
awk '{print $1, $2}' ${DATASET_INPUT}/${DATASET_NAME}_06_filt.fam > ${DATASET_NAME}_07_all_samples.txt
comm -23 \
  <(sort ${DATASET_NAME}_07_all_samples.txt) \
  <(sort ${DATASET_NAME}_07_keep_samples.txt) \
  > ${DATASET_NAME}_07_removal_list.txt

NREMOVE=$(wc -l < ${DATASET_NAME}_07_removal_list.txt)
echo "Related samples to remove: ${NREMOVE}"

# ============================================================================
# STEP 4: Create kinship-pruned dataset
# ============================================================================
echo ""
echo "=== Creating pruned (unrelated) dataset ==="
echo ""

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_06_filt \
  --remove ${DATASET_NAME}_07_removal_list.txt \
  --make-bed \
  --out ${DATASET_NAME}_07_filt

echo "✓ Relatedness-pruned dataset: ${DATASET_NAME}_07_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NSAMP_BEFORE=$(tail -n +2 ${DATASET_INPUT}/${DATASET_NAME}_06_filt.fam | wc -l)
NSAMP_AFTER=$(tail -n +2 ${DATASET_NAME}_07_filt.fam | wc -l)
NSAMP_REMOVED=$((NSAMP_BEFORE - NSAMP_AFTER))

echo "Samples before relatedness filter: ${NSAMP_BEFORE}"
echo "Samples after relatedness filter:  ${NSAMP_AFTER}"
echo "Related samples removed:          ${NSAMP_REMOVED}"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the MAF (minor allele frequency) filter:"
echo ""
echo "  bash 08_maf_filter.sh"
echo ""
echo "This will:"
echo "  - Apply MAF threshold (--maf, context-dependent)"
echo "  - For PDAC: use --maf 0.01 to retain rare variants"
echo "  - Create final clean dataset for association testing"
echo ""
