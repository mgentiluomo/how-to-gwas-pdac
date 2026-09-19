#!/usr/bin/env bash

################################################################################
# Section 5: Meta-analysis — Step 01: ancestry-stratified association testing
#
# PURPOSE:
#   Run the Section 4A analysis separately in each ancestry group, producing one
#   set of summary statistics per stratum. These are the inputs to the
#   meta-analysis.
#
# WHAT CHANGED IN THIS REVISION
#
#   1. INPUT DATASET. Reads results/pca/pdac_demo_02_hwe_filt, the dataset
#      Section 2 declares as carried into association testing, rather than
#      pdac_demo_08_filt, which precedes the within-ancestry Hardy-Weinberg
#      filter.
#
#   2. PRUNING. The per-stratum pruning now matches every other pruning call in
#      the pipeline: autosomes only, a minor allele frequency floor evaluated
#      inside the stratum, and long-range linkage disequilibrium regions
#      excluded. Without the region exclusion a handful of very long haplotype
#      blocks survive pruning and can produce one of the leading components,
#      which is precisely the failure the components are meant to correct.
#
#   3. VISIBILITY. The pruning and PCA output is no longer sent to /dev/null.
#      The number of variants each stratum's components were computed on is a
#      number the Methods section has to report, and suppressing it hides both
#      that and any warning PLINK issues.
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
#   - results/pca/pdac_demo_02_hwe_filt.bed/bim/fam
#   - demo_data/sample_ancestry.tsv, phenotype.txt, covariates.txt
#   - data_processed/highLD_b38.bed
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

IN="${1:-results/pca/pdac_demo_02_hwe_filt}"
ANC="${2:-demo_data/sample_ancestry.tsv}"
PHENO="${3:-demo_data/phenotype.txt}"
COVAR_IN="${4:-demo_data/covariates.txt}"
OUT_DIR="${5:-results/meta}"
LRLD="${6:-data_processed/highLD_b38.bed}"
PRUNE_MAF="${7:-0.05}"

# NOTE: do not name this variable GROUPS. In bash, GROUPS is a reserved array
# holding the current user's Unix group IDs; assigning to it fails silently and
# the loop below would iterate over a group ID instead of an ancestry label.
STRATA="eur afr eas"
NPC=5   # see the note below

mkdir -p "$OUT_DIR"

# --- long-range LD guard -----------------------------------------------------
# Fails loudly rather than proceeding without the exclusion: a silent no-op
# produces output indistinguishable from a correct run.
require_lrld() {
  local bed="$1" bim="$2" n regions
  if [ ! -s "$bed" ]; then
    echo "X Long-range LD region file not found: $bed" >&2
    echo "  Create it before running this step, or pass a path as argument 6." >&2
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
    echo "X No variants fall inside the long-range LD regions listed in $bed" >&2
    echo "  The exclusion would be a silent no-op. Usual causes:" >&2
    echo "    - chromosome codes differ ('1' in the .bim vs 'chr1' in the BED)" >&2
    echo "    - the BED is on a different build than the data (GRCh38 expected)" >&2
    echo "    - bed0 vs bed1 coordinate convention mismatch" >&2
    exit 1
  fi
  echo "Long-range LD exclusion: ${regions} regions, ${n} variants in scope"
}

[ -s "${IN}.bed" ] || {
  echo "X ${IN}.bed not found." >&2
  echo "  Section 2 must be run first; 06_hwe_within_ancestry.sh writes it." >&2
  exit 1
}
require_lrld "$LRLD" "${IN}.bim"

COUNTS="${OUT_DIR}/pdac_demo_05_strata_counts.tsv"
echo -e "stratum\tindividuals\tcases\tcontrols\teffective_n\tpca_variants" > "$COUNTS"

# A note on the number of components.
#   Section 4A fits the European stratum, of about 600 individuals, with ten
#   components and again with the subset that its screen retained. The African
#   and East Asian strata here have fewer than 320 each. Ten components
#   estimated in 300 people is a large number of parameters for the information
#   available, and over-adjustment costs power. Five is used for every stratum,
#   so that the strata are treated identically and the meta-analysis combines
#   like with like. Step 03 reports the inflation factor per stratum, which is
#   the check that five was enough.

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
  # --keep is applied before --maf and before --indep-pairwise, so both the
  # frequency floor and the LD estimates use this stratum's own genotypes.
  plink2 --bfile "$IN" --keep "$KEEP" \
         --autosome \
         --maf "$PRUNE_MAF" \
         --exclude bed0 "$LRLD" \
         --indep-pairwise 50 5 0.2 \
         --out "$PRUNE"

  NPRUNE=$(wc -l < "${PRUNE}.prune.in")
  echo "  ${G}: ${NPRUNE} pruned variants for the within-stratum PCA"

  plink2 --bfile "$IN" --keep "$KEEP" \
         --extract "${PRUNE}.prune.in" \
         --pca $NPC --out "$PCA"

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
  echo -e "${G}\t${NTOT}\t${NCA}\t${NCO}\t${NEFF}\t${NPRUNE}" >> "$COUNTS"
done

echo ""
echo "=== Stratum composition ==="
cat "$COUNTS"

echo ""
echo "Note how the effective sample size behaves. A stratum with an almost"
echo "balanced case-control ratio contributes nearly its full size, while an"
echo "unbalanced one contributes much less than its headline count suggests."
echo ""
echo "The pca_variants column is the number of pruned variants each stratum's"
echo "components were computed on. These differ between strata by design,"
echo "because linkage disequilibrium structure differs between populations."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Rscript scripts/05_meta_analysis/02_harmonise.R"
echo ""
