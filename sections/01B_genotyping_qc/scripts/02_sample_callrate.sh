#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 02: Sample Call Rate Filtering
# 
# PURPOSE:
#   Remove samples with excessive missing genotypes. High missingness can
#   indicate platform failures, lab issues, or poor DNA quality. This is the
#   first sample-level filter and is applied BEFORE variant filters because
#   low-quality samples skew allele frequency estimates.
#
# INPUT:
#   - pdac_demo.bed/bim/fam (raw genotypes)
#
# OUTPUT:
#   - pdac_demo.02_samples_callrate_pass.txt — list of samples passing filter
#
# THRESHOLD:
#   --mind 0.02  (allow max 2% missing genotypes per sample)
#   For 430K variants: 0.02 * 430K = 8,600 missing calls per sample is OK
#
# NOTES:
#   - Applied before any variant filters (avoids bias from rare/poor-quality variants)
#   - Context: PDAC has limited samples; we're stringent to retain power
#   - Individual samples are evaluated; variants are kept
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
DATASET_INPUT="${2:-demo_data}"
DATASET_NAME="pdac_demo"
OUT_DIR="${3:-results/qc}"

mkdir -p "$OUT_DIR"

# ============================================================================
# STEP 1: Filter samples by call rate
# ============================================================================
echo ""
echo "=== Filtering samples by call rate (--mind 0.02) ==="
echo ""

# Why --mind 0.02?
#   "mind" = "missing individual"
#   --mind 0.02 removes any sample with >2% missing genotypes.
#   This is a standard threshold that:
#   - Catches platform failures and DNA quality issues
#   - Retains power by not being too stringent
#   - Applied FIRST (before variant filters) to avoid skewed allele counts
#
# Why 2% and not stricter (e.g. 1%)?
#   In PDAC with ~250 cases, each sample is valuable for power.
#   2% balances data quality with sample retention.
#   If post-hoc QC reveals more outliers, case/control-specific thresholds can be used.

plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}" \
  --mind 0.02 \
  --write-samples \
  --out "${OUT_DIR}/${DATASET_NAME}_02_samples_callrate"

awk 'NR > 1 {print $1, $2}' "${OUT_DIR}/${DATASET_NAME}_02_samples_callrate.psam" > "${OUT_DIR}/${DATASET_NAME}_02_samples_keep.txt"

echo "✓ Samples passing call rate filter: $(wc -l < "${OUT_DIR}/${DATASET_NAME}_02_samples_keep.txt")"

# ============================================================================
# STEP 2: Create filtered dataset
# ============================================================================
echo ""
echo "=== Creating sample-filtered dataset ==="
echo ""

# Keep only samples that passed the filter
plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}" \
  --keep "${OUT_DIR}/${DATASET_NAME}_02_samples_keep.txt" \
  --make-bed \
  --out "${OUT_DIR}/${DATASET_NAME}_02_filt"

echo "✓ Filtered dataset: ${OUT_DIR}/${DATASET_NAME}_02_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NSAMP_BEFORE=$(wc -l < "${DATASET_INPUT}/${DATASET_NAME}.fam")
NSAMP_AFTER=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_02_filt.fam")
NSAMP_REMOVED=$((NSAMP_BEFORE - NSAMP_AFTER))
PERCENT_RETAINED=$(awk -v after="$NSAMP_AFTER" -v before="$NSAMP_BEFORE" 'BEGIN {printf "%.1f", after * 100 / before}')

echo "Samples before filter:  ${NSAMP_BEFORE}"
echo "Samples after filter:   ${NSAMP_AFTER}"
echo "Samples removed:        ${NSAMP_REMOVED}"
echo "% retained:             ${PERCENT_RETAINED}%"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the sex check:"
echo ""
echo "  bash scripts/01B_genotyping_qc/03_sex_check.sh"
echo ""
echo "This will:"
echo "  - Compute sex-specific allele frequencies (chrX)"
echo "  - Identify sex-discordant samples"
echo "  - Flag samples for removal if sex is miscoded or discordant with genotypes"
echo ""
