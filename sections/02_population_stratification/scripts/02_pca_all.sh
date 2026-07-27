#!/usr/bin/env bash

################################################################################
# Section 2: Population stratification — Step 02: PCA on the full QC-passed set
#
# PURPOSE:
#   Compute principal components on all QC-passed individuals, across all three
#   ancestry groups. This is the picture that shows what the cohort actually
#   contains, and it is the evidence on which the decision to restrict the
#   primary analysis to a single ancestry rests.
#
#   PCA makes no assumption about population structure. It reduces the genotype
#   matrix to a set of uncorrelated axes ordered by the variance they explain;
#   ancestry simply happens to be the dominant biological signal in that
#   variance when a cohort spans continents.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.bed/bim/fam
#   - results/pca/pdac_demo_02_prune.prune.in   (from Step 01)
#
# OUTPUT:
#   - results/pca/pdac_demo_02_pca_all.eigenvec  one row per individual
#   - results/pca/pdac_demo_02_pca_all.eigenval  variance per component
#
# NOTE:
#   Chromosome X is excluded from the genetic relationship matrix automatically.
#   In this dataset chrX is simulated for the sex check only, so it must not
#   contribute to ancestry inference.
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

IN="${1:-results/qc/pdac_demo_08_filt}"
OUT_DIR="${2:-results/pca}"
PRUNE_IN="${OUT_DIR}/pdac_demo_02_prune.prune.in"
OUT_PREFIX="${OUT_DIR}/pdac_demo_02_pca_all"

mkdir -p "$OUT_DIR"

echo ""
echo "=== Principal component analysis, all ancestry groups ==="
echo ""

# Why ten components?
#   Ten is enough to see where the signal stops. How many to *use* as covariates
#   is decided later, from the scree plot and from testing each component
#   against case status, not by adopting a default.
plink2 \
  --bfile "$IN" \
  --extract "$PRUNE_IN" \
  --pca 10 \
  --out "$OUT_PREFIX"

echo ""
echo "Outputs:"
echo "  ${OUT_PREFIX}.eigenvec   PC scores, one row per individual"
echo "  ${OUT_PREFIX}.eigenval   eigenvalues, one per component"
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Rscript scripts/02_population_stratification/03_pca_plots.R"
echo ""
