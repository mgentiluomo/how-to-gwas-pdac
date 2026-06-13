#!/usr/bin/env bash

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
#   - PDAC phenotype is in demo_data/phenotype.txt
#   - Disease-associated variants may show HWE deviation in cases (expected)
#   - Standard threshold: 1e-6 is stringent; 1e-4 is more lenient
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "$SCRIPT_DIR/../../scripts/dev" ]; then
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
elif [ -d "$SCRIPT_DIR/../../../scripts/dev" ]; then
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
else
  PROJECT_ROOT="$(pwd)"
fi
cd "$PROJECT_ROOT"

# Configuration
SEED="${1:-2026}"
DATASET_INPUT="${2:-results/qc}"
DATASET_NAME="pdac_demo"
PHENOTYPE_FILE="${3:-demo_data/phenotype.txt}"
OUT_DIR="${4:-results/qc}"

mkdir -p "$OUT_DIR"

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

awk '($3 == 1 || $3 == "1") {print $1, $2}' "$PHENOTYPE_FILE" > "${OUT_DIR}/${DATASET_NAME}_06_controls.txt"

NCTRL=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_06_controls.txt")
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
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_05_filt" \
  --keep "${OUT_DIR}/${DATASET_NAME}_06_controls.txt" \
  --hardy \
  --out "${OUT_DIR}/${DATASET_NAME}_06_hwe"

echo "✓ HWE test results: ${OUT_DIR}/${DATASET_NAME}_06_hwe.hardy"

# ============================================================================
# STEP 3: Identify HWE-failing variants
# ============================================================================
echo ""
echo "=== Identifying HWE violations ==="
echo ""

# Extract variants failing HWE (p < 1e-6)
# hardy file format: CHR SNP TEST NOBS OBS_CT EXP_CT P
# We want: TEST == "ALL" (or "UNAFF" for unaffected/controls) and P < 1e-6

awk '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      gsub(/^#/, "", $i)
      if ($i == "ID") id = i
      if ($i == "P") p = i
    }
    next
  }
  id && p && $p != "NA" && $p < 1e-6 { print $id }
' "${OUT_DIR}/${DATASET_NAME}_06_hwe.hardy" > "${OUT_DIR}/${DATASET_NAME}_06_hwe_exclude.txt"

NHWE_FAIL=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_06_hwe_exclude.txt")
echo "Variants failing HWE (p < 1e-6): ${NHWE_FAIL}"

# ============================================================================
# STEP 4: Create HWE-filtered dataset
# ============================================================================
echo ""
echo "=== Removing HWE-failing variants ==="
echo ""

if [ -s "${OUT_DIR}/${DATASET_NAME}_06_hwe_exclude.txt" ]; then
  plink2 \
    --bfile "${DATASET_INPUT}/${DATASET_NAME}_05_filt" \
    --exclude "${OUT_DIR}/${DATASET_NAME}_06_hwe_exclude.txt" \
    --make-bed \
    --out "${OUT_DIR}/${DATASET_NAME}_06_filt"
else
  plink2 \
    --bfile "${DATASET_INPUT}/${DATASET_NAME}_05_filt" \
    --make-bed \
    --out "${OUT_DIR}/${DATASET_NAME}_06_filt"
fi

echo "✓ HWE-filtered dataset: ${OUT_DIR}/${DATASET_NAME}_06_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NVAR_BEFORE=$(wc -l < "${DATASET_INPUT}/${DATASET_NAME}_05_filt.bim")
NVAR_AFTER=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_06_filt.bim")
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
echo "  bash scripts/01B_genotyping_qc/07_relatedness.sh"
echo ""
echo "This will:"
echo "  - Compute kinship coefficients (IBD, PI_HAT)"
echo "  - Identify related pairs (2nd degree or closer)"
echo "  - Prune to one sample per related pair"
echo ""
