#!/usr/bin/env bash

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
DATASET_DIR="${2:-demo_data}"
DATASET_NAME="pdac_demo"
OUT_DIR="${3:-results/qc}"
OUT_PREFIX="${OUT_DIR}/${DATASET_NAME}_01_qc"

mkdir -p "$OUT_DIR"

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
  --bfile "${DATASET_DIR}/${DATASET_NAME}" \
  --freq \
  --out "$OUT_PREFIX"

echo "✓ Allele frequencies saved to ${OUT_PREFIX}.afreq"

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
  --bfile "${DATASET_DIR}/${DATASET_NAME}" \
  --het \
  --out "$OUT_PREFIX"

echo "✓ Heterozygosity computed, saved to ${OUT_PREFIX}.het"

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
  --bfile "${DATASET_DIR}/${DATASET_NAME}" \
  --missing \
  --out "$OUT_PREFIX"

echo "✓ Missing data rates saved to:"
echo "  - ${OUT_PREFIX}.smiss (per sample)"
echo "  - ${OUT_PREFIX}.vmiss (per variant)"

# ============================================================================
# STEP 4: Inspect basic summary statistics
# ============================================================================
echo ""
echo "=== Computing summary statistics ==="
echo ""

# Variant count
NVAR=$(wc -l < "${DATASET_DIR}/${DATASET_NAME}.bim")

# Sample count
NSAMP=$(wc -l < "${DATASET_DIR}/${DATASET_NAME}.fam")

# Median heterozygosity
MEDIAN_HET=$(awk '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      gsub(/^#/, "", $i)
      if ($i == "F") f_col = i
    }
    next
  }
  f_col && $f_col != "nan" { print $f_col }
' "${OUT_PREFIX}.het" | sort -n | awk '{a[NR]=$1} END {if (NR) print a[int((NR + 1) / 2)]; else print "NA"}')

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
echo "  bash scripts/01B_genotyping_qc/02_sample_callrate.sh"
echo ""
echo "This will apply filters for:"
echo "  - Variant call rate (--geno)"
echo "  - Hardy-Weinberg equilibrium (--hwe)"
echo "  - Minor allele frequency (--maf) [optional]"
echo ""
