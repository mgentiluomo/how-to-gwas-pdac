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
#   variants are correlated, differ between populations. --keep is applied
#   before both --maf and --indep-pairwise, so the frequency floor and the LD
#   estimates both use the analysis set's own genotypes rather than the pooled
#   multi-ancestry ones.
#
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.bed/bim/fam
#   - results/pca/pdac_demo_02_eur_keep.txt   (from Step 04)
#   - data_processed/highLD_b38.bed           (long-range LD regions, GRCh38)
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
LRLD_BED="${3:-data_processed/highLD_b38.bed}"
PRUNE_MAF="${4:-0.05}"
KEEP="${OUT_DIR}/pdac_demo_02_eur_keep.txt"
PRUNE="${OUT_DIR}/pdac_demo_02_prune_eur"
OUT_PREFIX="${OUT_DIR}/pdac_demo_02_pca_eur"

# ---------------------------------------------------------------------------
# Long-range LD guard. Fails loudly rather than proceeding without the
# exclusion, because a silent no-op produces output indistinguishable from a
# correct run and would put the Methods text out of step with the code.
# ---------------------------------------------------------------------------
require_lrld() {
  local bed="$1" bim="$2" n regions
  if [ ! -s "$bed" ]; then
    echo "✗ Long-range LD region file not found: $bed" >&2
    echo "  Create it before running this step, or pass a path as argument 3." >&2
    exit 1
  fi
  n=$(awk '
    FNR == NR {
      if ($0 ~ /^#/ || NF < 3) next
      k++; c[k] = $1; s[k] = $2 + 0; e[k] = $3 + 0
      next
    }
    { for (i = 1; i <= k; i++) if ($1 == c[i] && $4 > s[i] && $4 <= e[i]) { hit++; break } }
    END { print hit + 0 }
  ' "$bed" "$bim")
  regions=$(awk '!/^#/ && NF >= 3' "$bed" | wc -l)
  if [ "$n" -eq 0 ]; then
    echo "✗ No variants fall inside the long-range LD regions listed in $bed" >&2
    echo "  The exclusion would be a silent no-op. Usual causes:" >&2
    echo "    - chromosome codes differ ('1' in the .bim vs 'chr1' in the BED)" >&2
    echo "    - the BED is on a different build than the data (GRCh38 expected)" >&2
    echo "    - bed0 vs bed1 coordinate convention mismatch" >&2
    exit 1
  fi
  echo "Long-range LD exclusion: ${regions} regions, ${n} variants removed before pruning"
}

if [ ! -s "$KEEP" ]; then
  echo "✗ Analysis set keep list not found: $KEEP" >&2
  echo "  Run Step 04 first:" >&2
  echo "    bash scripts/02_population_stratification/04_define_analysis_set.sh" >&2
  exit 1
fi

echo ""
echo "=== LD pruning within the analysis set ==="
echo ""

require_lrld "$LRLD_BED" "${IN}.bim"
echo ""

plink2 \
  --bfile "$IN" \
  --keep "$KEEP" \
  --autosome \
  --maf "$PRUNE_MAF" \
  --exclude bed0 "$LRLD_BED" \
  --indep-pairwise 50 5 0.2 \
  --out "$PRUNE"

NVAR=$(wc -l < "${PRUNE}.prune.in")
echo ""
echo "Variants retained for the within-group PCA: ${NVAR}"

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
echo "Report ${NVAR} as the number of variants the components were computed on."
echo "If PC3 or PC4 has moved appreciably against a run without the long-range"
echo "LD exclusion, that component was tracking a haplotype block rather than"
echo "ancestry, and the comparison is worth recording."
echo ""
echo "How many components to carry forward is decided in Section 4A, from this"
echo "scree plot, from testing each component against case status, and from"
echo "checking that the results are stable across a range of counts. Note that"
echo "--pca 10 above caps extraction at ten components: raise it if the"
echo "sensitivity check calls for comparing 5, 10 and 20."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Section 4A: Rscript scripts/04A_association_binary/01_make_covariates.R"
echo ""
