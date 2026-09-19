#!/usr/bin/env bash

################################################################################
# Section 4A: Association testing — Step 02: logistic regression, genome-wide
#
# PURPOSE:
#   Test every variant for association with case status under an additive model,
#   adjusted for sex, age and the principal components from Section 2.
#
# WHAT CHANGED IN THIS REVISION
#
#   1. INPUT DATASET. The previous version read results/qc/pdac_demo_08_filt,
#      which is the dataset before the within-ancestry Hardy-Weinberg filter.
#      Section 2 ends by writing results/pca/pdac_demo_02_hwe_filt and declaring
#      it the dataset carried into association testing, so that is what is read
#      here. The difference is small on the demonstration data, but a quality
#      control step whose output nothing consumes has not been applied.
#
#   2. TWO COVARIATE MODELS. Step 01 screens each component for association with
#      case status and prints the result. Passing all ten components regardless
#      makes that screen decorative. This script now fits both models:
#
#        primary     SEX, AGE, PC1-PC10        the conventional choice
#        screened    SEX, AGE, associated PCs  the components the screen kept
#
#      Neither is automatically correct. Ten components in a few hundred
#      individuals is a large number of parameters for the information
#      available, and over-adjustment costs power; but a component that is not
#      nominally associated with status may still be absorbing structure that
#      matters, and dropping components on a P < 0.05 screen is itself a
#      selection procedure. Fitting both and reporting the comparison is the
#      honest response, and it supplies part of the stability check the guide
#      asks for.
#
#      The primary model keeps the original output prefix, so the downstream
#      scripts (03_qq_lambda.R, 04_manhattan.R, 05_summary.R) are unaffected.
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
#   Note that this makes the analysed variant count smaller than the post-QC
#   count reported in Table 2. Both numbers should be stated.
#
# WHY --covar-variance-standardize:
#   Puts covariates on a common scale, which improves convergence. It does not
#   change the genotype effect estimates.
#
# INPUT:
#   - results/pca/pdac_demo_02_hwe_filt.bed/bim/fam
#   - results/pca/pdac_demo_02_eur_keep.txt
#   - demo_data/phenotype.txt
#   - results/assoc/pdac_demo_04A_covar.txt
#   - results/assoc/pdac_demo_04A_pc_screen.tsv
#
# OUTPUT:
#   - results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid        primary
#   - results/assoc/pdac_demo_04A_gwas_pcscreened.PHENO.glm.*           sensitivity
#   - results/assoc/pdac_demo_04A_pc_model_comparison.tsv
#   - the .log files: keep these, they record the case and control counts
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
KEEP="${2:-results/pca/pdac_demo_02_eur_keep.txt}"
PHENO="${3:-demo_data/phenotype.txt}"
OUT_DIR="${4:-results/assoc}"
NPC_PRIMARY="${5:-10}"
COVAR="${OUT_DIR}/pdac_demo_04A_covar.txt"
SCREEN="${OUT_DIR}/pdac_demo_04A_pc_screen.tsv"
OUT_PRIMARY="${OUT_DIR}/pdac_demo_04A_gwas"
OUT_SCREENED="${OUT_DIR}/pdac_demo_04A_gwas_pcscreened"
COMPARISON="${OUT_DIR}/pdac_demo_04A_pc_model_comparison.tsv"

mkdir -p "$OUT_DIR"

for f in "${IN}.bed" "$KEEP" "$PHENO" "$COVAR"; do
  [ -s "$f" ] || { echo "X required input missing: $f" >&2; exit 1; }
done

if [ ! -s "$IN.bed" ]; then
  echo "X $IN.bed not found." >&2
  echo "  Section 2 must be run first; 06_hwe_within_ancestry.sh writes it." >&2
  exit 1
fi

# --- build the two covariate lists -------------------------------------------
PC_PRIMARY=$(seq -s, -f 'PC%g' 1 "$NPC_PRIMARY")
COVAR_PRIMARY="SEX,AGE,${PC_PRIMARY}"

if [ -s "$SCREEN" ]; then
  # Column 1 is the PC name, the last column is yes/no. Read by header name so
  # the parse survives a change in column order.
  PC_SCREENED=$(awk -F'\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) { if ($i == "PC") pc = i; if ($i == "nominally_associated") a = i }
      next
    }
    $a == "yes" { printf "%s%s", (n++ ? "," : ""), $pc }
  ' "$SCREEN")
else
  echo "!! $SCREEN not found; skipping the screened model." >&2
  PC_SCREENED=""
fi

if [ -n "$PC_SCREENED" ]; then
  COVAR_SCREENED="SEX,AGE,${PC_SCREENED}"
else
  COVAR_SCREENED="SEX,AGE"
fi

echo ""
echo "=== Covariate models ==="
echo ""
echo "  primary   ${COVAR_PRIMARY}"
echo "  screened  ${COVAR_SCREENED}"
if [ -z "$PC_SCREENED" ]; then
  echo ""
  echo "  No component reached P < 0.05 in the Step 01 screen, so the screened"
  echo "  model adjusts for sex and age alone. That is a legitimate result on a"
  echo "  single-ancestry analysis set, not a failure."
fi

# --- primary model -----------------------------------------------------------
echo ""
echo "=== Genome-wide association test: primary model ==="
echo ""

plink2 \
  --bfile "$IN" \
  --chr 1-22 \
  --keep "$KEEP" \
  --pheno "$PHENO" \
  --covar "$COVAR" \
  --covar-name "$COVAR_PRIMARY" \
  --covar-variance-standardize \
  --glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err \
  --out "$OUT_PRIMARY"

echo ""
echo "Check the log for the line reporting how many cases and controls were"
echo "analysed. That number, not the number of individuals in the file, is what"
echo "the Methods section must report:"
echo ""
grep -E "cases and .* controls remaining" "${OUT_PRIMARY}.log" || true

# --- screened model ----------------------------------------------------------
echo ""
echo "=== Genome-wide association test: screened model ==="
echo ""

plink2 \
  --bfile "$IN" \
  --chr 1-22 \
  --keep "$KEEP" \
  --pheno "$PHENO" \
  --covar "$COVAR" \
  --covar-name "$COVAR_SCREENED" \
  --covar-variance-standardize \
  --glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err \
  --out "$OUT_SCREENED"

# --- compare the two ---------------------------------------------------------
echo ""
echo "=== Comparing the two covariate models ==="
echo ""

Rscript --vanilla -e "
suppressWarnings({
  find_glm <- function(prefix) {
    f <- Sys.glob(paste0(prefix, '*.glm.logistic.hybrid'))
    if (!length(f)) f <- Sys.glob(paste0(prefix, '*.glm.logistic'))
    if (!length(f)) stop('no glm output found for ', prefix)
    f[1]
  }
  load_res <- function(prefix) {
    d <- read.table(find_glm(prefix), header = TRUE, sep = '\t',
                    comment.char = '', check.names = FALSE)
    names(d)[1] <- sub('^#', '', names(d)[1])
    d
  }
  lambda <- function(p) {
    p <- p[is.finite(p) & p > 0 & p <= 1]
    stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) /
      stats::qchisq(0.5, 1, lower.tail = FALSE)
  }
  summarise <- function(d, label, covars) {
    p <- suppressWarnings(as.numeric(d\$P))
    ok <- is.finite(p)
    id <- if ('ID' %in% names(d)) d\$ID else rep(NA_character_, nrow(d))
    lead <- if (any(ok)) id[which.min(replace(p, !ok, Inf))] else NA_character_
    leadp <- if (any(ok)) min(p[ok]) else NA_real_
    data.frame(
      model            = label,
      covariates       = covars,
      variants_tested  = nrow(d),
      valid_P          = sum(ok),
      no_valid_P       = sum(!ok),
      lambda           = round(lambda(p), 4),
      genome_wide_sig  = sum(ok & p < 5e-8),
      suggestive       = sum(ok & p < 1e-5),
      lead_variant     = lead,
      lead_P           = signif(leadp, 3),
      stringsAsFactors = FALSE
    )
  }

  a <- load_res('${OUT_PRIMARY}')
  b <- load_res('${OUT_SCREENED}')
  cmp <- rbind(summarise(a, 'primary',  '${COVAR_PRIMARY}'),
               summarise(b, 'screened', '${COVAR_SCREENED}'))
  write.table(cmp, '${COMPARISON}', sep = '\t', quote = FALSE, row.names = FALSE)
  print(cmp[, c('model','variants_tested','no_valid_P','lambda',
                'genome_wide_sig','suggestive','lead_variant','lead_P')],
        row.names = FALSE)

  # How far do the two models move the individual estimates?
  key <- intersect(names(a), 'ID')
  if (length(key) && 'BETA' %in% names(a) && 'BETA' %in% names(b)) {
    m <- merge(a[, c('ID','BETA','P')], b[, c('ID','BETA','P')],
               by = 'ID', suffixes = c('.primary','.screened'))
    m <- m[is.finite(m\$BETA.primary) & is.finite(m\$BETA.screened), ]
    if (nrow(m) > 100) {
      cat('\n  Pearson r between the two sets of effect estimates: ',
          round(stats::cor(m\$BETA.primary, m\$BETA.screened), 5), '\n', sep = '')
      cat('  Largest absolute difference in log odds: ',
          signif(max(abs(m\$BETA.primary - m\$BETA.screened)), 3), '\n', sep = '')
      cat('  Genome-wide significant in one model only: ',
          sum(xor(m\$P.primary < 5e-8, m\$P.screened < 5e-8)), '\n', sep = '')
    }
  }
})
"

echo ""
echo "Written: ${COMPARISON}"
echo ""
echo "Report both models. If the inflation factor and the leading associations"
echo "are stable across them, the choice of component count did not drive the"
echo "result, which is the claim the guide asks to be substantiated rather than"
echo "assumed. If they differ, that difference is itself the finding and the"
echo "primary model must be justified rather than defaulted to."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Rscript scripts/04A_association_binary/03_qq_lambda.R"
echo ""
