#!/usr/bin/env bash

################################################################################
# Section 2: Population stratification — Step 05: PCA within the analysis set
#
# PURPOSE:
#   Recompute principal components using only the individuals who will actually
#   be analysed.
#
#   This step is frequently skipped, and skipping it is a mistake. The
#   components from Step 02 describe continental separation: they are dominated
#   by the differences between Europe, Africa and East Asia. Once the non-target
#   groups are removed, those axes no longer describe anything present in the
#   data. The structure that remains, and that can still confound the
#   association test, is the fine-scale variation within the retained group, and
#   only a PCA computed within that group can see it.
#
#   LD pruning is also repeated, because allele frequencies, and therefore which
#   variants are correlated, differ between populations.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.bed/bim/fam
#   - results/pca/pdac_demo_02_eur_keep.txt   (from Step 04)
#
# OUTPUT:
#   - results/pca/pdac_demo_02_pca_eur.eigenvec   covariates for Section 4A
#   - results/pca/pdac_demo_02_pca_eur.eigenval
#   - results/pca/pdac_demo_02_pca_eur_scree.png
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
KEEP="${OUT_DIR}/pdac_demo_02_eur_keep.txt"
PRUNE="${OUT_DIR}/pdac_demo_02_prune_eur"
OUT_PREFIX="${OUT_DIR}/pdac_demo_02_pca_eur"

echo ""
echo "=== LD pruning within the analysis set ==="
echo ""

plink2 \
  --bfile "$IN" \
  --keep "$KEEP" \
  --indep-pairwise 50 5 0.2 \
  --out "$PRUNE"

echo ""
echo "=== PCA within the analysis set ==="
echo ""

plink2 \
  --bfile "$IN" \
  --keep "$KEEP" \
  --extract "${PRUNE}.prune.in" \
  --pca 10 \
  --out "$OUT_PREFIX"

# A scree plot for the within-group PCA. Compare it with the one from Step 03:
# the continental structure is gone, and what remains is a much flatter curve.
Rscript --vanilla -e "
val <- scan('${OUT_PREFIX}.eigenval', quiet = TRUE)
pct <- 100 * val / sum(val)
png('${OUT_PREFIX}_scree.png', width = 1600, height = 1200, res = 200)
plot(seq_along(pct), pct, type = 'b', pch = 19, col = '#0072B2',
     xlab = 'Principal component', ylab = 'Variance explained (%)',
     main = 'Scree plot within the analysis set')
dev.off()
cat('Variance explained, PC1 to PC4:', sprintf('%.2f%%', pct[1:4]), '\n')
"

echo ""
echo "Outputs:"
echo "  ${OUT_PREFIX}.eigenvec       covariates for association testing"
echo "  ${OUT_PREFIX}.eigenval"
echo "  ${OUT_PREFIX}_scree.png"
echo ""
echo "How many components to carry forward is decided in Section 4A, from this"
echo "scree plot, from testing each component against case status, and from"
echo "checking that the results are stable across a range of counts."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Section 4A: Rscript scripts/04A_association_binary/01_make_covariates.R"
echo ""
