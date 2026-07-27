#!/usr/bin/env bash

################################################################################
# Section 4A: Association testing — Step 02: logistic regression, genome-wide
#
# PURPOSE:
#   Test every variant for association with case status under an additive model,
#   adjusted for sex, age and the principal components from Section 2.
#
# WHY --glm firth-fallback:
#   Standard logistic regression becomes unreliable when a variant is rare, when
#   cases are few, or when a genotype almost perfectly separates cases from
#   controls: the estimate diverges and the P value cannot be trusted. Firth's
#   penalised likelihood corrects this. The fallback form uses ordinary logistic
#   regression by default and switches to Firth only where it fails, which is
#   both faster and, for the variants where it matters, more accurate.
#
#   In a rare disease this is not an optional refinement. With a few hundred
#   cases, separation is common at exactly the low-frequency variants that are
#   of most interest.
#
# WHY --chr 1-22:
#   Association testing here is restricted to the autosomes. The chromosome X in
#   this demonstration dataset is simulated so that the sex check in Section 1B
#   has data to work on; it carries no phenotype signal. PLINK 2 also models sex
#   separately on chrX, which collides with sex supplied as a user covariate and
#   will stop the run with a CORR_TOO_HIGH error. Analysing the X chromosome
#   properly is a real task with its own rules; it is out of scope here.
#
# WHY --covar-variance-standardize:
#   Puts covariates on a common scale, which improves convergence. It does not
#   change the genotype effect estimates.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.bed/bim/fam
#   - results/pca/pdac_demo_02_eur_keep.txt
#   - demo_data/phenotype.txt
#   - results/assoc/pdac_demo_04A_covar.txt
#
# OUTPUT:
#   - results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid
#   - results/assoc/pdac_demo_04A_gwas.log      keep this: it records the
#                                               case and control counts analysed
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
KEEP="${2:-results/pca/pdac_demo_02_eur_keep.txt}"
PHENO="${3:-demo_data/phenotype.txt}"
OUT_DIR="${4:-results/assoc}"
COVAR="${OUT_DIR}/pdac_demo_04A_covar.txt"
OUT_PREFIX="${OUT_DIR}/pdac_demo_04A_gwas"

mkdir -p "$OUT_DIR"

echo ""
echo "=== Genome-wide association test ==="
echo ""

plink2 \
  --bfile "$IN" \
  --chr 1-22 \
  --keep "$KEEP" \
  --pheno "$PHENO" \
  --covar "$COVAR" \
  --covar-name SEX,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --covar-variance-standardize \
  --glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err \
  --out "$OUT_PREFIX"

echo ""
echo "Check the log for the line reporting how many cases and controls were"
echo "analysed. That number, not the number of individuals in the file, is what"
echo "the Methods section must report:"
echo ""
grep -E "cases and .* controls remaining" "${OUT_PREFIX}.log" || true
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Rscript scripts/04A_association_binary/03_qq_lambda.R"
echo ""
