#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 08: Minor Allele Frequency (MAF) Filter
# 
# PURPOSE:
#   Remove or keep variants based on minor allele frequency. MAF filtering is
#   context-dependent and depends on:
#   - Study design (discovery vs replication)
#   - Sample size and power
#   - Disease rarity and heritability architecture
#
# INPUT:
#   - pdac_demo_07_filt.bed/bim/fam (from Step 07, relatedness-pruned)
#
# OUTPUT:
#   - pdac_demo_08_filt.bed/bim/fam — MAF-filtered genotypes (ready for association!)
#
# THRESHOLD FOR PDAC:
#   --maf 0.01  (keep variants with MAF ≥ 1%)
#   NOT --maf 0.05 (which would discard rare variants with large effects)
#
# RATIONALE FOR LOW-MAF THRESHOLD IN PDAC:
#   1. Rare variants carry larger effect sizes in complex diseases
#   2. Limited sample size: we need all signals (can't afford to discard rare variants)
#   3. Pancreatic cancer is rare (~2% lifetime risk): causal variants may be rare
#   4. Common variant + rare variant mixed model explains most heritability
#
# NOTES:
#   - Alternative thresholds:
#     - --maf 0.05: classic GWAS (assumes common-disease common-variant)
#     - --maf 0.001 (or --mac 3): discovery in very small cohorts
#   - If computational power is limited, can increase to 0.05 later
#   - Statistical power decreases with lower MAF; adjust sample size expectations
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
OUT_DIR="${3:-results/qc}"
MAF_THRESHOLD="${4:-0.01}"  # Default: 1% (PDAC context)

mkdir -p "$OUT_DIR"

# ============================================================================
# STEP 1: Apply MAF filter
# ============================================================================
echo ""
echo "=== Applying MAF filter (--maf ${MAF_THRESHOLD}) ==="
echo ""

# Why --maf 0.01 for PDAC?
#   Standard GWAS uses --maf 0.05 (common-disease common-variant model).
#   PDAC context requires a different approach:
#
#   - Pancreatic cancer is RARE (2% lifetime risk)
#   - Genetic architecture likely includes rare-to-intermediate variants
#   - Removing variants with 1-5% MAF discards potential signal
#   - Sample size is limited (~250 EUR cases): power is already low
#     → Retaining rare variants is necessary to detect effects
#
# MAF = minor allele frequency = frequency of less-common allele
# --maf 0.01 keeps variants where the rarer allele appears in ≥1% of the sample

plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_07_filt" \
  --maf "${MAF_THRESHOLD}" \
  --make-bed \
  --out "${OUT_DIR}/${DATASET_NAME}_08_filt"

echo "✓ MAF-filtered dataset: ${OUT_DIR}/${DATASET_NAME}_08_filt.bed/bim/fam"

# ============================================================================
# STEP 2: Display filtered allele spectrum
# ============================================================================
echo ""
echo "=== Allele frequency spectrum (post-QC) ==="
echo ""

# Compute final allele frequencies
plink2 \
  --bfile "${OUT_DIR}/${DATASET_NAME}_08_filt" \
  --freq \
  --out "${OUT_DIR}/${DATASET_NAME}_08_afreq"

# Summary statistics on final dataset
echo "Final dataset summary:"
NVAR_FINAL=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_08_filt.bim")
NSAMP_FINAL=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_08_filt.fam")

# Count variants by MAF bins
echo ""
echo "Variants by MAF bin (post-QC):"
awk 'NR > 1 {
  maf = $5 < 0.5 ? $5 : 1 - $5
  if (maf < 0.001) bin = "< 0.1%"
  else if (maf < 0.01) bin = "0.1%-1%"
  else if (maf < 0.05) bin = "1%-5%"
  else bin = "> 5%"
  count[bin]++
}
END {
  for (b in count) print "  " b ": " count[b]
}' "${OUT_DIR}/${DATASET_NAME}_08_afreq.afreq" | sort

echo ""
echo "Total variants: ${NVAR_FINAL}"
echo "Total samples:  ${NSAMP_FINAL}"
echo ""

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NVAR_BEFORE=$(wc -l < "${DATASET_INPUT}/${DATASET_NAME}_07_filt.bim")
NVAR_AFTER=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_08_filt.bim")
NVAR_REMOVED=$((NVAR_BEFORE - NVAR_AFTER))

echo "Variants before MAF filter: ${NVAR_BEFORE}"
echo "Variants after MAF filter:  ${NVAR_AFTER}"
echo "Variants removed:           ${NVAR_REMOVED}"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the final QC summary:"
echo ""
echo "  bash scripts/01B_genotyping_qc/09_qc_summary.sh"
echo ""
echo "This will:"
echo "  - Summarize total samples and variants filtered at each step"
echo "  - Show filtering impact on sample/variant counts"
echo "  - Create a decision-tree figure for reproducibility and documentation"
echo ""
echo "After that, your dataset is ready for:"
echo "  - Population stratification analysis (Section 2: PCA)"
echo "  - Imputation (Section 3)"
echo "  - Association testing (Section 4)"
echo ""
