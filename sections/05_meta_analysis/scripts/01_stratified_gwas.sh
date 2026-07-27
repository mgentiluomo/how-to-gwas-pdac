#!/usr/bin/env bash

################################################################################
# Section 5: Meta-analysis — Step 01: ancestry-stratified association testing
#
# PURPOSE:
#   Run the Section 4A analysis separately in each ancestry group, producing one
#   set of summary statistics per stratum. These are the inputs to the
#   meta-analysis.
#
# WHY STRATIFY RATHER THAN POOL:
#   Pooling individuals of different ancestries into one regression assumes that
#   the effect of each variant is the same in every group, and that principal
#   components alone can absorb the frequency differences between them. Neither
#   assumption is safe. Causal-allele frequencies and linkage disequilibrium
#   structure differ between populations, so the same marker can tag a causal
#   variant well in one group and poorly in another.
#
#   Analysing each group separately and combining the results afterwards makes
#   no such assumption. It also makes any difference between groups visible and
#   testable, through the heterogeneity statistics computed in Step 03, rather
#   than hidden inside a single pooled estimate.
#
# EACH STRATUM GETS ITS OWN:
#   - LD pruning, because which variants are correlated depends on the
#     population;
#   - principal components, because the axes of variation within Africa are not
#     the axes of variation within Europe;
#   - association test.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.bed/bim/fam
#   - demo_data/sample_ancestry.tsv, phenotype.txt, covariates.txt
#
# OUTPUT, for each of eur, afr, eas:
#   - results/meta/<group>/pdac_demo_05_<group>_gwas.PHENO.glm.logistic.hybrid
#   - results/meta/<group>/pdac_demo_05_<group>_gwas.log
#   - results/meta/pdac_demo_05_strata_counts.tsv
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
ANC="${2:-demo_data/sample_ancestry.tsv}"
PHENO="${3:-demo_data/phenotype.txt}"
COVAR_IN="${4:-demo_data/covariates.txt}"
OUT_DIR="${5:-results/meta}"

# NOTE: do not name this variable GROUPS. In bash, GROUPS is a reserved array
# holding the current user's Unix group IDs; assigning to it fails silently and
# the loop below would iterate over a group ID instead of an ancestry label.
STRATA="eur afr eas"
NPC=5   # see the note below

mkdir -p "$OUT_DIR"
COUNTS="${OUT_DIR}/pdac_demo_05_strata_counts.tsv"
echo -e "stratum\tindividuals\tcases\tcontrols\teffective_n" > "$COUNTS"

# A note on the number of components.
#   Section 4A used ten in the European stratum, which has 625 individuals. The
#   African and East Asian strata have fewer than 320 each. Ten components
#   estimated in 300 people is a large number of parameters for the information
#   available, and over-adjustment costs power. Five is used here for every
#   stratum, so that the strata are treated identically and the meta-analysis
#   combines like with like. Step 03 reports the inflation factor per stratum,
#   which is the check that five was enough.

for G in $STRATA; do
  echo ""
  echo "############################################################"
  echo "###  Stratum: ${G}"
  echo "############################################################"

  GDIR="${OUT_DIR}/${G}"
  mkdir -p "$GDIR"
  KEEP="${GDIR}/keep.txt"
  PRUNE="${GDIR}/prune"
  PCA="${GDIR}/pca"
  COVAR="${GDIR}/covar.txt"
  GWAS="${GDIR}/pdac_demo_05_${G}_gwas"

  # --- who is in this stratum ------------------------------------------------
  awk -v target="$G" '
    NR==FNR { if (tolower($2) == target) k[$1] = 1; next }
    ($2 in k) { print $1"\t"$2 }
  ' "$ANC" "${IN}.fam" > "$KEEP"

  echo "Individuals in stratum: $(wc -l < "$KEEP")"

  # --- LD pruning and PCA, within this stratum -------------------------------
  plink2 --bfile "$IN" --keep "$KEEP" \
         --indep-pairwise 50 5 0.2 --out "$PRUNE" > /dev/null

  plink2 --bfile "$IN" --keep "$KEEP" \
         --extract "${PRUNE}.prune.in" \
         --pca $NPC --out "$PCA" > /dev/null

  # --- covariates ------------------------------------------------------------
  Rscript --vanilla -e "
    cov <- read.table('${COVAR_IN}', header = TRUE, comment.char = '', check.names = FALSE)
    pcs <- read.table('${PCA}.eigenvec', header = TRUE, comment.char = '', check.names = FALSE)
    names(cov)[1] <- sub('^#', '', names(cov)[1])
    names(pcs)[1] <- sub('^#', '', names(pcs)[1])
    m <- merge(cov, pcs, by = intersect(c('FID','IID'), intersect(names(cov), names(pcs))))
    names(m)[1] <- paste0('#', names(m)[1])
    write.table(m, '${COVAR}', sep = '\t', quote = FALSE, row.names = FALSE)
    cat('covariate rows:', nrow(m), '\n')
  "

  # --- association -----------------------------------------------------------
  PC_LIST=$(seq -s, -f 'PC%g' 1 $NPC)
  plink2 --bfile "$IN" \
         --chr 1-22 \
         --keep "$KEEP" \
         --pheno "$PHENO" \
         --covar "$COVAR" \
         --covar-name "SEX,AGE,${PC_LIST}" \
         --covar-variance-standardize \
         --glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err \
         --out "$GWAS"

  # --- record the composition ------------------------------------------------
  LINE=$(grep -E "cases and .* controls remaining" "${GWAS}.log" | head -1)
  NCA=$(echo "$LINE" | grep -oE "^[0-9]+")
  NCO=$(echo "$LINE" | grep -oE "and [0-9]+" | grep -oE "[0-9]+")
  NTOT=$(wc -l < "$KEEP")
  NEFF=$(awk -v a="$NCA" -v b="$NCO" 'BEGIN{ if (a+b>0) printf "%.1f", 4*a*b/(a+b); else print "NA" }')
  echo -e "${G}\t${NTOT}\t${NCA}\t${NCO}\t${NEFF}" >> "$COUNTS"
done

echo ""
echo "=== Stratum composition ==="
cat "$COUNTS"

echo ""
echo "Note how the effective sample size behaves. A stratum with an almost"
echo "balanced case-control ratio contributes nearly its full size, while an"
echo "unbalanced one contributes much less than its headline count suggests."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Rscript scripts/05_meta_analysis/02_harmonise.R"
echo ""
