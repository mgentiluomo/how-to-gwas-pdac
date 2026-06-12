#!/bin/bash

################################################################################
# Section 1B: Genotyping QC — Step 06: Hardy-Weinberg Equilibrium (HWE) Test
# 
# PURPOSE:
#   Remove variants that deviate significantly from Hardy-Weinberg equilibrium.
#   HWE violations can indicate:
#   - Genotyping errors (false heterozygotes or homozygotes)
#   - Population stratification
#   - Selection pressure
#
# INPUT:
#   - pdac_demo_05_filt.bed/bim/fam (from Step 05, geno-filtered)
#   - phenotype.txt (case/control status)
#
# OUTPUT:
#   - pdac_demo_06_filt.bed/bim/fam — HWE-filtered genotypes
#
# THRESHOLD:
#   --hwe 1e-6 (p > 1e-6 in controls; more lenient in cases)
#
# STRATEGY:
#   1. Test controls ONLY for HWE (p > 1e-6)
#   2. Caution: cases may deviate due to disease association
#   3. Remove variants failing HWE in controls
#   4. Keep all variants passing in controls (even if cases deviate)
#
# NOTES:
#   - PDAC phenotype is in ../../demo_dataset/data/phenotype.txt
#   - Disease-associated variants may show HWE deviation in cases (expected)
#   - Standard threshold: 1e-6 is stringent; 1e-4 is more lenient
#
################################################################################

set -e

# Configuration
SEED="${1:-2026}"
DATASET_INPUT="${2:-.}"
DATASET_NAME="pdac_demo"
PHENOTYPE_FILE="../../demo_dataset/data/phenotype.txt"

# ============================================================================
# STEP 1: Extract control IDs
# ============================================================================
echo ""
echo "=== Testing Hardy-Weinberg equilibrium in controls ==="
echo ""

# Why test controls only?
#   HWE is based on the assumption of no selection, no mutation, no migration.
#   In PDAC GWAS, *disease variants* may violate HWE in CASES (expected).
#   Testing in controls only:
#   - Isolates genotyping errors (which violate HWE in controls)
#   - Avoids removing disease-associated variants (which deviate in cases)
#   - Follows standard GWAS practice
#
# PHENOTYPE file format: FID IID PHENO (1=control, 2=case)

echo "Extracting control samples..."

awk '$3 == 1 {print $1, $2}' ${PHENOTYPE_FILE} > ${DATASET_NAME}_controls.txt

NCTRL=$(wc -l < ${DATASET_NAME}_controls.txt)
echo "Control samples: ${NCTRL}"

# ============================================================================
# STEP 2: Test HWE in controls
# ============================================================================
echo ""
echo "=== Computing HWE p-values in controls ==="
echo ""

# Why --hwe 1e-6?
#   HWE is a statistical test; with many variants, some will fail by chance.
#   Bonferroni correction: α = 0.05 / number_of_variants
#   With 420K variants: α ≈ 1e-7 (very stringent)
#   Using 1e-6 is conservative but allows small p-values from genotyping errors.
#
# If you want a more lenient threshold (to retain more variants):
#   Use --hwe 1e-4 (or even 1e-3 for discovery)

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_05_filt \
  --keep ${DATASET_NAME}_controls.txt \
  --hardy \
  --out ${DATASET_NAME}_06_hwe

echo "✓ HWE test results: ${DATASET_NAME}_06_hwe.hardy"

# ============================================================================
# STEP 3: Identify HWE-failing variants
# ============================================================================
echo ""
echo "=== Identifying HWE violations ==="
echo ""

# Extract variants failing HWE (p < 1e-6)
# hardy file format: CHR SNP TEST NOBS OBS_CT EXP_CT P
# We want: TEST == "ALL" (or "UNAFF" for unaffected/controls) and P < 1e-6

awk 'NR > 1 && $7 < 1e-6 {print $2}' ${DATASET_NAME}_06_hwe.hardy > ${DATASET_NAME}_06_hwe_exclude.txt

NHWE_FAIL=$(wc -l < ${DATASET_NAME}_06_hwe_exclude.txt)
echo "Variants failing HWE (p < 1e-6): ${NHWE_FAIL}"

# ============================================================================
# STEP 4: Create HWE-filtered dataset
# ============================================================================
echo ""
echo "=== Removing HWE-failing variants ==="
echo ""

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_05_filt \
  --exclude ${DATASET_NAME}_06_hwe_exclude.txt \
  --make-bed \
  --out ${DATASET_NAME}_06_filt

echo "✓ HWE-filtered dataset: ${DATASET_NAME}_06_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NVAR_BEFORE=$(tail -n +2 ${DATASET_INPUT}/${DATASET_NAME}_05_filt.bim | wc -l)
NVAR_AFTER=$(tail -n +2 ${DATASET_NAME}_06_filt.bim | wc -l)
NVAR_REMOVED=$((NVAR_BEFORE - NVAR_AFTER))

echo "Variants before HWE filter: ${NVAR_BEFORE}"
echo "Variants after HWE filter:  ${NVAR_AFTER}"
echo "HWE-failing variants:       ${NVAR_REMOVED}"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the relatedness check:"
echo ""
echo "  bash 07_relatedness.sh"
echo ""
echo "This will:"
echo "  - Compute kinship coefficients (IBD, PI_HAT)"
echo "  - Identify related pairs (2nd degree or closer)"
echo "  - Prune to one sample per related pair"
echo ""
