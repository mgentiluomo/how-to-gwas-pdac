#!/usr/bin/env bash

################################################################################
# Section 6: Fine mapping — Step 01: define the region and measure LD
#
# PURPOSE:
#   An association signal is not one variant. Because nearby variants are
#   inherited together, a single causal variant drags its neighbours up with it,
#   and the association test cannot tell them apart. Fine mapping asks which of
#   the correlated variants is most likely to be the cause, and returns a
#   credible set: a shortlist that, under the assumptions of the method, contains
#   the causal variant with a stated probability.
#
#   This step prepares the two things every fine mapping method needs: the
#   association statistics for a region, and a matrix of correlations between the
#   variants in it.
#
# CHOOSING THE REGION:
#   Wide enough to contain the whole signal including its flanks, narrow enough
#   that the single causal variant assumption is plausible. A window of a few
#   hundred kilobases either side of the lead variant is conventional. Too narrow
#   and the true causal variant may fall outside it; too wide and unrelated
#   signals are drawn in.
#
# WHICH LD:
#   In sample LD, computed from the genotypes actually analysed, is the correct
#   choice whenever individual level data are available, and it is what is used
#   here. When only summary statistics exist, an external reference panel is
#   substituted, and it must match the ancestry of the study. A mismatched LD
#   reference is the commonest cause of a confidently wrong credible set.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt              QC-passed genotypes
#   - results/pca/pdac_demo_02_eur_keep.txt     the analysis set
#
# OUTPUT:
#   - results/finemap/abo_eur.bed/bim/fam       genotypes for the region
#   - results/finemap/abo_ld.vcor               pairwise r2
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

QC="${1:-results/qc/pdac_demo_08_filt}"
KEEP="${2:-results/pca/pdac_demo_02_eur_keep.txt}"
OUT_DIR="${3:-results/finemap}"
CHR=9
FROM=132900000
TO=133600000

mkdir -p "$OUT_DIR"

echo ""
echo "=== Extracting the region: chr${CHR}:${FROM}-${TO} ==="
echo ""
plink2 --bfile "$QC" --keep "$KEEP" \
       --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
       --make-bed --out "${OUT_DIR}/abo_eur"

N=$(wc -l < "${OUT_DIR}/abo_eur.bim")
echo ""
echo "Variants in the region: ${N}"
echo ""
echo "Note how few they are. On a genotyping array a 700 kb window carries a"
echo "couple of hundred markers, and the resolution of any fine mapping is"
echo "bounded by that. Imputation increases the count but does not by itself"
echo "increase the information, as Section 3 explains."
echo ""

echo "=== Pairwise linkage disequilibrium ==="
echo ""
plink2 --bfile "${OUT_DIR}/abo_eur" \
       --r2-unphased --ld-window-r2 0 \
       --out "${OUT_DIR}/abo_ld"

echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Rscript scripts/06_finemapping_annotation/02_finemap_abf.R"
echo ""
