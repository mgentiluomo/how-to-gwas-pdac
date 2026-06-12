#!/bin/bash

################################################################################
# Section 1B: Genotyping QC — Step 05: Variant Call Rate Filtering
# 
# PURPOSE:
#   Remove variants with excessive missing genotypes. Variants with high
#   missingness indicate:
#   - Assay design problems (SNP not well-genotyped on platform)
#   - Systematic technical failures
#   - Unreliable allele frequency estimates
#
# INPUT:
#   - pdac_demo_04_filt.bed/bim/fam (from Step 04, het-filtered samples)
#
# OUTPUT:
#   - pdac_demo_05_filt.bed/bim/fam — genotypes with variant filter applied
#
# THRESHOLD:
#   --geno 0.05 (remove variants with >5% missing)
#
# NOTES:
#   - Applied AFTER sample filters to get accurate missing-data patterns
#   - Standard threshold, often lowered to 0.02 for more stringent QC
#   - Variants with poor genotyping across all samples are unreliable
#
################################################################################

set -e

# Configuration
SEED="${1:-2026}"
DATASET_INPUT="${2:-.}"
DATASET_NAME="pdac_demo"

# ============================================================================
# STEP 1: Filter variants by call rate
# ============================================================================
echo ""
echo "=== Filtering variants by call rate (--geno 0.05) ==="
echo ""

# Why --geno 0.05?
#   "geno" = "genotype missingness"
#   --geno 0.05 removes any variant with >5% missing genotypes.
#   This threshold:
#   - Catches problematic SNPs (assay failures, polymorphism issues)
#   - Retains 95% of genotype calls per variant (high quality)
#   - Balances stringency with variant retention
#
# Why 5% and not more stringent (e.g. 1%)?
#   Variants with 1-5% missingness are usually fine.
#   Very stringent thresholds can remove true variants with rare technical issues.
#   Post-QC, you can always be more stringent if power is not a concern.

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_04_filt \
  --geno 0.05 \
  --make-bed \
  --out ${DATASET_NAME}_05_filt

echo "✓ Variant-filtered dataset: ${DATASET_NAME}_05_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NVAR_BEFORE=$(tail -n +2 ${DATASET_INPUT}/${DATASET_NAME}_04_filt.bim | wc -l)
NVAR_AFTER=$(tail -n +2 ${DATASET_NAME}_05_filt.bim | wc -l)
NVAR_REMOVED=$((NVAR_BEFORE - NVAR_AFTER))
PERCENT_RETAINED=$(echo "scale=1; $NVAR_AFTER * 100 / $NVAR_BEFORE" | bc)

echo "Variants before geno filter: ${NVAR_BEFORE}"
echo "Variants after geno filter:  ${NVAR_AFTER}"
echo "Variants removed:            ${NVAR_REMOVED}"
echo "% retained:                  ${PERCENT_RETAINED}%"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the Hardy-Weinberg equilibrium test:"
echo ""
echo "  bash 06_hardy_weinberg.sh"
echo ""
echo "This will:"
echo "  - Test variants for deviation from Hardy-Weinberg equilibrium"
echo "  - Apply separate filters for cases and controls"
echo "  - (Disease variants may deviate in cases)"
echo ""
