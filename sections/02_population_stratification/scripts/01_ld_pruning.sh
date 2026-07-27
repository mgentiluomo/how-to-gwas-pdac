#!/usr/bin/env bash

################################################################################
# Section 2: Population stratification — Step 01: LD pruning
#
# PURPOSE:
#   Principal component analysis must be computed on variants that are
#   approximately independent of one another. Variants in linkage disequilibrium
#   (LD) carry overlapping information, and dense regions such as the major
#   histocompatibility complex on chromosome 6 would otherwise dominate the
#   components, so that the leading axes describe local genomic architecture
#   rather than ancestry.
#
#   This step keeps roughly one representative variant per LD block.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.bed/bim/fam   (QC-passed data, Section 1B)
#
# OUTPUT:
#   - results/pca/pdac_demo_02_prune.prune.in    (variants to keep)
#   - results/pca/pdac_demo_02_prune.prune.out   (variants pruned away)
#
# PARAMETERS (--indep-pairwise 50 5 0.2):
#   50   window size, in variants
#   5    step size: the window advances five variants at a time
#   0.2  r-squared threshold; within a window, one of any pair above this is
#        dropped
#
#   The r-squared threshold is a genuine decision, not a standard. 0.1 is more
#   stringent (|r| about 0.32), 0.2 more permissive (|r| about 0.45). What is
#   required is to report the value used and to confirm that the leading
#   components do not change materially when it is varied.
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
OUT_PREFIX="${OUT_DIR}/pdac_demo_02_prune"

mkdir -p "$OUT_DIR"

echo ""
echo "=== LD pruning ==="
echo ""

plink2 \
  --bfile "$IN" \
  --indep-pairwise 50 5 0.2 \
  --out "$OUT_PREFIX"

KEPT=$(wc -l < "${OUT_PREFIX}.prune.in")
DROPPED=$(wc -l < "${OUT_PREFIX}.prune.out")

echo ""
echo "Variants retained for PCA:  $KEPT"
echo "Variants pruned away:       $DROPPED"
echo ""
echo "A substantial fraction is expected to be pruned: this is redundancy being"
echo "removed, not data being lost. The pruned set is used only for PCA and"
echo "relatedness, never for association testing."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  bash scripts/02_population_stratification/02_pca_all.sh"
echo ""
