#!/bin/bash

################################################################################
# Section 1B: Genotyping QC — Step 03: Sex Check
# 
# PURPOSE:
#   Identify and flag samples with sex discordance (genetic sex does not match
#   reported sex). This catches:
#   - Sample swaps
#   - Data entry errors
#   - Potential contamination (XY in supposed females, XX in males)
#
# INPUT:
#   - pdac_demo_02_filt.bed/bim/fam (from Step 02, sample call rate filtered)
#
# OUTPUT:
#   - pdac_demo.03_sexcheck.sexcheck — sex concordance report
#   - pdac_demo.03_sexcheck_discordant.txt — list of discordant samples to remove
#
# METHOD:
#   Compute X chromosome heterozygosity (F statistic on chrX).
#   - Females (XX): low heterozygosity (F near 0)
#   - Males (XY): high heterozygosity (F near 1)
#   Samples with F inconsistent with reported sex are flagged.
#
# NOTES:
#   - Requires sex information in .fam file (column 4: 1=M, 2=F)
#   - Our demo dataset has simulated chrX; real data uses genotyped chrX
#   - Standard thresholds: F_LOW < 0.2 (females), F_HIGH > 0.8 (males)
#
################################################################################

set -e

# Configuration
SEED="${1:-2026}"
DATASET_INPUT="${2:-.}"
DATASET_NAME="pdac_demo"

# ============================================================================
# STEP 1: Compute X-chromosome genotyping frequency (F statistic)
# ============================================================================
echo ""
echo "=== Computing X chromosome heterozygosity (sex check) ==="
echo ""

# Why --check-sex?
#   PLINK2 computes the inbreeding coefficient F on the X chromosome only.
#   - Females (XX) genotypes: F near 0 (mostly heterozygous relative to population)
#   - Males (XY) genotypes: F near 1 (all X alleles homozygous, haploid)
#
# The rationale:
#   Males have one X chromosome, so all X markers are "homozygous" (haploid).
#   Females have two X chromosomes, so they can be heterozygous.
#   If a sample labeled "female" shows high F (close to 1), it's likely mislabeled
#   or swapped. Conversely, a sample labeled "male" with low F is suspicious.

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_02_filt \
  --check-sex \
  --out ${DATASET_NAME}_03_sexcheck

echo "✓ Sex check complete: ${DATASET_NAME}_03_sexcheck.sexcheck"

# ============================================================================
# STEP 2: Identify discordant samples
# ============================================================================
echo ""
echo "=== Identifying sex-discordant samples ==="
echo ""

# Extract discordant samples (STATUS = "PROBLEM")
# Standard thresholds:
#   - Females (reported sex 2): F should be < 0.2 (mostly heterozygous on X)
#   - Males (reported sex 1): F should be > 0.8 (mostly homozygous on X, haploid)

echo "#FID IID" > ${DATASET_NAME}_03_sexcheck_discordant.txt
awk '$5 == "PROBLEM" {print $1, $2}' ${DATASET_NAME}_03_sexcheck.sexcheck >> ${DATASET_NAME}_03_sexcheck_discordant.txt

NDISCORDANT=$(tail -n +2 ${DATASET_NAME}_03_sexcheck_discordant.txt | wc -l)

echo "Discordant samples identified: ${NDISCORDANT}"

if [ $NDISCORDANT -gt 0 ]; then
  echo ""
  echo "Samples to review:"
  head -20 ${DATASET_NAME}_03_sexcheck_discordant.txt
fi

# ============================================================================
# STEP 3: Create dataset without discordant samples
# ============================================================================
echo ""
echo "=== Removing discordant samples ==="
echo ""

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_02_filt \
  --remove ${DATASET_NAME}_03_sexcheck_discordant.txt \
  --make-bed \
  --out ${DATASET_NAME}_03_filt

echo "✓ Sex-filtered dataset: ${DATASET_NAME}_03_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NSAMP_BEFORE=$(tail -n +2 ${DATASET_INPUT}/${DATASET_NAME}_02_filt.fam | wc -l)
NSAMP_AFTER=$(tail -n +2 ${DATASET_NAME}_03_filt.fam | wc -l)
NSAMP_REMOVED=$((NSAMP_BEFORE - NSAMP_AFTER))

echo "Samples before sex check:      ${NSAMP_BEFORE}"
echo "Samples after sex check:       ${NSAMP_AFTER}"
echo "Sex-discordant samples removed: ${NSAMP_REMOVED}"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the heterozygosity check:"
echo ""
echo "  bash 04_heterozygosity.sh"
echo ""
echo "This will:"
echo "  - Compute per-sample heterozygosity (F coefficient on autosomes)"
echo "  - Identify outliers (contamination, inbreeding, population outliers)"
echo "  - Generate visualization plots"
echo ""
