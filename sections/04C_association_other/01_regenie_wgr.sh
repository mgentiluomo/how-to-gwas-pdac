#!/usr/bin/env bash

################################################################################
# Whole-genome regression: modelling relatedness instead of pruning
#
# THE QUESTION
#   Relatedness pruning removed 121 individuals from the European cohort. In a
#   rare-trait study, people removed are information lost, so the removal has to
#   be justified rather than assumed. The alternative is to keep the related
#   individuals and model the relatedness, which is what REGENIE's whole-genome
#   regression does.
#
#   Running REGENIE on the post-QC dataset would answer nothing: the relatives
#   have already been removed, so there is no relatedness left to model. The
#   informative comparison goes back to the dataset before relatedness pruning
#   and asks what the pruning bought.
#
# FOUR ANALYSES, THE SAME VARIANTS THROUGHOUT
#
#   A  pruned set,   unrelated only,      logistic regression (plink2)
#   B  unpruned set, relatedness ignored,  logistic regression (plink2)
#   C  unpruned set, relatedness modelled, REGENIE with the step 1 offset
#   D  unpruned set, relatedness ignored,  REGENIE without the offset
#
#   A against B is the total cost of the pruning decision. Note that it is a
#   compound comparison: the two sets differ in size AND in case-control
#   balance, because the phenotype-aware rule sacrificed controls
#   preferentially. Effective sample size is reported so that the balance shift
#   is visible rather than buried.
#
#   B against C is the comparison that actually tests the mixed model: the same
#   individuals, the same balance, relatedness ignored against modelled.
#
#   D separates the software from the method. If C differs from B, D tells you
#   whether that is the leave-one-chromosome-out offset doing work or simply
#   REGENIE and plink2 computing slightly different things. Without D the two
#   are confounded. Set RUN_D=0 to skip it.
#
# A CAVEAT SPECIFIC TO THIS DEMONSTRATION
#   The phenotype was simulated from two causal variants with no polygenic
#   background. REGENIE's step 1 fits a whole-genome ridge regression to capture
#   polygenic signal, so with none present the fitted predictor has little to
#   fit and will largely absorb the two causal variants instead. The step 1
#   cross-validated R-squared is printed below for exactly this reason: it is
#   the number that tells you how much the offset is carrying, and whether C is
#   a fair test of the method or an artefact of the simulation.
#
# INPUT
#   results/qc/pdac_demo_06_filt.bed/bim/fam    pre-pruning dataset
#   results/pca/pdac_demo_02_hwe_exclude.txt    within-ancestry HWE failures
#   results/pca/pdac_demo_02_eur_keep.txt       the pruned Europeans
#   results/assoc/pdac_demo_04A_covar.txt       covariates for the pruned set
#   demo_data/sample_ancestry.tsv, phenotype.txt, covariates.txt
#   data_processed/highLD_b38.bed
#
# OUTPUT (default results/regenie, override with GWAS_OUT_DIR)
#   pdac_demo_regenie_A_pruned.*      pdac_demo_regenie_B_unpruned.*
#   pdac_demo_regenie_C_step2_*       pdac_demo_regenie_D_nopred_*
#   pdac_demo_regenie_comparison.tsv  pdac_demo_regenie_summary.txt
#
# USAGE
#   bash 01_regenie_regenie.sh [n_pcs] [maf_floor] [prune_maf]
#   Resources are detected at run time; override with GWAS_THREADS, GWAS_MEM_MB.
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

N_PC="${1:-10}"
MAF_MIN="${2:-0.01}"
PRUNE_MAF="${3:-0.05}"
RUN_D="${RUN_D:-1}"

IN_UNPRUNED="results/qc/pdac_demo_06_filt"
HWE_EXCLUDE="results/pca/pdac_demo_02_hwe_exclude.txt"
KEEP_PRUNED="results/pca/pdac_demo_02_eur_keep.txt"
COVAR_PRUNED="results/assoc/pdac_demo_04A_covar.txt"
ANC="demo_data/sample_ancestry.tsv"
PHENO="demo_data/phenotype.txt"
COVAR_IN="demo_data/covariates.txt"
LRLD="data_processed/highLD_b38.bed"
OUT_DIR="${GWAS_OUT_DIR:-results/regenie}"
TMP_DIR="${OUT_DIR}/tmp"

mkdir -p "$OUT_DIR" "$TMP_DIR"

# ============================================================================
# REGENIE writes its step 1 predictor list as a whitespace-delimited file with
# no quoting, so any output path containing a space makes step 2 fail with
# "step 1 list file is not in the right format". Project directories under
# OneDrive or "My Documents" hit this routinely. REGENIE therefore writes to a
# scratch directory guaranteed to be free of spaces, and the results are copied
# back afterwards. Everything else, including the input paths, is quoted and
# handles spaces normally.
# ============================================================================
SCRATCH="${GWAS_SCRATCH:-${HOME}/.cache/gwas_regenie}"
case "$SCRATCH" in
  *" "*) echo "X scratch path contains a space: $SCRATCH" >&2
         echo "  set GWAS_SCRATCH to a path without spaces" >&2; exit 1 ;;
esac
mkdir -p "$SCRATCH"
rm -f "${SCRATCH}"/pdac_demo_regenie_*

# ============================================================================
# Resource detection
#
# A scheduler's allocation always wins over what the hardware reports: on a
# shared node the two differ and only the first is yours to use.
# ============================================================================
detect_cores() {
  local v n q p c
  for v in GWAS_THREADS SLURM_CPUS_PER_TASK NSLOTS PBS_NUM_PPN OMP_NUM_THREADS; do
    n="${!v:-}"
    if [ -n "$n" ] && [ "$n" -gt 0 ] 2>/dev/null; then echo "$n"; return; fi
  done
  n="$(nproc 2>/dev/null || echo 1)"
  if [ -r /sys/fs/cgroup/cpu.max ]; then
    read -r q p < /sys/fs/cgroup/cpu.max || true
    if [ "${q:-max}" != "max" ] && [ "${p:-0}" -gt 0 ] 2>/dev/null; then
      c=$(( q / p )); [ "$c" -ge 1 ] && [ "$c" -lt "$n" ] && n="$c"
    fi
  fi
  [ "$n" -gt 1 ] && n=$(( n - 1 ))
  echo "$n"
}

detect_mem_mb() {
  local m lim
  if [ -n "${GWAS_MEM_MB:-}" ]; then echo "$GWAS_MEM_MB"; return; fi
  m="$(awk '/^MemAvailable:/{print int($2/1024)}' /proc/meminfo 2>/dev/null || true)"
  [ -z "$m" ] && m=2048
  if [ -r /sys/fs/cgroup/memory.max ]; then
    lim="$(cat /sys/fs/cgroup/memory.max)"
    if [ "$lim" != "max" ]; then
      lim=$(( lim / 1024 / 1024 ))
      [ "$lim" -gt 0 ] && [ "$lim" -lt "$m" ] && m="$lim"
    fi
  fi
  echo "$m"
}

NCORE="$(detect_cores)"
MEM_MB="$(detect_mem_mb)"
PLINK_MEM=$(( MEM_MB / 2 ))
# plink2 refuses --memory below 640; omit the flag rather than claim memory the
# machine does not have, and let plink2 detect for itself.
if [ "$PLINK_MEM" -ge 640 ]; then PLINK_MEM_FLAG="--memory $PLINK_MEM"; else PLINK_MEM_FLAG=""; fi
# REGENIE holds one block of genotypes at a time; the block size controls its
# footprint. Scale it, within sane bounds.
BSIZE1=1000; BSIZE2=400
if [ "$MEM_MB" -lt 2048 ];  then BSIZE1=400;  BSIZE2=200; fi
if [ "$MEM_MB" -gt 16384 ]; then BSIZE1=2000; BSIZE2=800; fi

echo ""
echo "=== Resources ==="
echo "  cores:  $NCORE"
echo "  memory: ${MEM_MB} MB   plink2 ${PLINK_MEM_FLAG:---memory auto}"
echo "  regenie block size: step 1 ${BSIZE1}, step 2 ${BSIZE2}"
echo "  regenie scratch: ${SCRATCH}"
echo "  outputs: ${OUT_DIR}"
echo "  override with GWAS_THREADS, GWAS_MEM_MB, GWAS_SCRATCH, GWAS_OUT_DIR"
echo ""

# ============================================================================
# Preflight
# ============================================================================
fail() { echo ""; echo "X $*" >&2; exit 1; }

for f in "${IN_UNPRUNED}.bed" "$KEEP_PRUNED" "$COVAR_PRUNED" "$ANC" "$PHENO" \
         "$COVAR_IN" "$LRLD"; do
  [ -s "$f" ] || fail "required input missing: $f"
done
command -v plink2 >/dev/null 2>&1 || fail "plink2 not on PATH"

REGENIE_CMD="regenie"
if ! command -v regenie >/dev/null 2>&1; then
  if command -v micromamba >/dev/null 2>&1; then
    REGENIE_CMD="micromamba run -n regenie_env regenie"
  else
    fail "regenie not on PATH. Re-run scripts/dev/tools_setup.sh"
  fi
fi
$REGENIE_CMD --version >/dev/null 2>&1 || fail "regenie will not run. If it was
  installed with micromamba, the symlink may not carry its shared libraries:
    micromamba run -n regenie_env regenie --version"
echo "regenie: $($REGENIE_CMD --version 2>&1 | head -1)"

# The flag for 1=control / 2=case encoding is --cc12 in REGENIE v3 and later,
# and was --1 in earlier releases. Detect rather than assume, because a
# miscoded binary trait produces plausible output that is entirely wrong.
CC_FLAG="${CC_FLAG:-}"
if [ -z "$CC_FLAG" ]; then
  HELP_TXT="$($REGENIE_CMD --help 2>&1 || true)"
  if printf '%s' "$HELP_TXT" | grep -q -- "--cc12"; then
    CC_FLAG="--cc12"
  elif printf '%s' "$HELP_TXT" | grep -q -- "--1 "; then
    CC_FLAG="--1"
  else
    CC_FLAG="--cc12"
    echo "  could not detect the encoding flag; defaulting to --cc12"
  fi
fi
echo "case/control encoding flag: ${CC_FLAG}"

# ============================================================================
# Analysis sets
# ============================================================================
echo ""
echo "=== Analysis sets ==="

KEEP_UNPRUNED="${TMP_DIR}/eur_unpruned.txt"
awk 'NR==FNR { if (tolower($2) == "eur") k[$1] = 1; next }
     ($2 in k) { print $1"\t"$2 }' "$ANC" "${IN_UNPRUNED}.fam" > "$KEEP_UNPRUNED"

N_UNPRUNED=$(wc -l < "$KEEP_UNPRUNED")
N_PRUNED=$(wc -l < "$KEEP_PRUNED")
N_REMOVED=$(( N_UNPRUNED - N_PRUNED ))
[ "$N_REMOVED" -gt 0 ] || fail "the two sets are the same size; check the inputs"

# What the pruning actually cost, by phenotype. On this dataset the
# phenotype-aware rule protected cases almost completely, so the loss falls on
# controls; that is a result about the rule and belongs in the output.
awk 'NR==FNR { p[$2] = $3; next }
     { seen[$2] = 1 }
     END { }' "$PHENO" "$KEEP_PRUNED" > /dev/null 2>&1 || true

REMOVED_BREAKDOWN=$(awk '
  FNR==NR { if (FNR > 1 || $1 != "FID") pheno[$2] = $3; next }
  FILENAME == kp { kept[$2] = 1; next }
  { if (!($2 in kept)) { if (pheno[$2] == 2) ca++; else if (pheno[$2] == 1) co++; else un++ } }
  END { printf "%d %d %d", ca+0, co+0, un+0 }
' "$PHENO" kp="$KEEP_PRUNED" "$KEEP_PRUNED" "$KEEP_UNPRUNED")
R_CASE=$(echo "$REMOVED_BREAKDOWN" | cut -d' ' -f1)
R_CTRL=$(echo "$REMOVED_BREAKDOWN" | cut -d' ' -f2)

echo "  pruned (unrelated) Europeans:   ${N_PRUNED}"
echo "  unpruned Europeans:             ${N_UNPRUNED}"
echo "  removed by relatedness pruning: ${N_REMOVED}  (${R_CASE} cases, ${R_CTRL} controls)"
echo ""
echo "  The two sets differ in case-control balance as well as in size, so the"
echo "  A against B comparison measures the whole pruning decision rather than"
echo "  relatedness alone. Effective sample size is reported for each."

# ============================================================================
# One variant set for all four analyses
# ============================================================================
echo ""
echo "=== Common variant set ==="

EXCLUDE_FLAG=""
[ -s "$HWE_EXCLUDE" ] && EXCLUDE_FLAG="--exclude $HWE_EXCLUDE"

plink2 --bfile "$IN_UNPRUNED" --keep "$KEEP_UNPRUNED" --chr 1-22 \
       --maf "$MAF_MIN" $EXCLUDE_FLAG \
       --threads "$NCORE" $PLINK_MEM_FLAG \
       --write-snplist --out "${TMP_DIR}/analysis_variants"

VARIANTS="${TMP_DIR}/analysis_variants.snplist"
N_VAR=$(wc -l < "$VARIANTS")
echo "  variants tested in all analyses: ${N_VAR}"

# ============================================================================
# Components within the unpruned set
#
# B, C and D cannot reuse the covariates built for the pruned set: those
# components were estimated in a different sample, and using them would
# confound the relatedness comparison with a covariate difference.
# ============================================================================
echo ""
echo "=== Components within the unpruned set ==="

plink2 --bfile "$IN_UNPRUNED" --keep "$KEEP_UNPRUNED" \
       --autosome --maf "$PRUNE_MAF" --exclude bed0 "$LRLD" \
       --indep-pairwise 50 5 0.2 \
       --threads "$NCORE" $PLINK_MEM_FLAG \
       --out "${TMP_DIR}/prune_unpruned"

N_PRUNE_VAR=$(wc -l < "${TMP_DIR}/prune_unpruned.prune.in")
echo "  pruned variants for the PCA and for REGENIE step 1: ${N_PRUNE_VAR}"

plink2 --bfile "$IN_UNPRUNED" --keep "$KEEP_UNPRUNED" \
       --extract "${TMP_DIR}/prune_unpruned.prune.in" --pca "$N_PC" \
       --threads "$NCORE" $PLINK_MEM_FLAG \
       --out "${TMP_DIR}/pca_unpruned"

COVAR_UNPRUNED="${TMP_DIR}/covar_unpruned.txt"
Rscript --vanilla -e "
  cov <- read.table('${COVAR_IN}', header = TRUE, comment.char = '', check.names = FALSE)
  pcs <- read.table('${TMP_DIR}/pca_unpruned.eigenvec', header = TRUE,
                    comment.char = '', check.names = FALSE)
  names(cov)[1] <- sub('^#', '', names(cov)[1])
  names(pcs)[1] <- sub('^#', '', names(pcs)[1])
  m <- merge(cov, pcs, by = intersect(c('FID','IID'), intersect(names(cov), names(pcs))))
  m <- m[complete.cases(m), ]
  write.table(m, '${COVAR_UNPRUNED}', sep = '\t', quote = FALSE, row.names = FALSE)
  cat('  covariate rows for the unpruned set:', nrow(m), '\n')
"

PC_LIST=$(seq -s, -f 'PC%g' 1 "$N_PC")

# ============================================================================
# A. Pruned set, logistic regression
# ============================================================================
echo ""
echo "############################################################"
echo "###  A: unrelated only, logistic regression"
echo "############################################################"

plink2 --bfile "$IN_UNPRUNED" --keep "$KEEP_PRUNED" --extract "$VARIANTS" \
       --pheno "$PHENO" --covar "$COVAR_PRUNED" \
       --covar-name "SEX,AGE,${PC_LIST}" --covar-variance-standardize \
       --glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err \
       --threads "$NCORE" $PLINK_MEM_FLAG \
       --out "${OUT_DIR}/pdac_demo_regenie_A_pruned"

# ============================================================================
# B. Unpruned set, relatedness ignored
# ============================================================================
echo ""
echo "############################################################"
echo "###  B: relatives retained, relatedness ignored"
echo "############################################################"

plink2 --bfile "$IN_UNPRUNED" --keep "$KEEP_UNPRUNED" --extract "$VARIANTS" \
       --pheno "$PHENO" --covar "$COVAR_UNPRUNED" \
       --covar-name "SEX,AGE,${PC_LIST}" --covar-variance-standardize \
       --glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err \
       --threads "$NCORE" $PLINK_MEM_FLAG \
       --out "${OUT_DIR}/pdac_demo_regenie_B_unpruned"

# ============================================================================
# C. Unpruned set, relatedness modelled
#
# Step 1 fits ridge regressions genome-wide on the pruned variant set and
# writes a leave-one-chromosome-out predictor per individual. Step 2 uses that
# predictor as an offset, so whatever polygenic and relatedness structure it
# captured is accounted for at every variant.
# ============================================================================
echo ""
echo "############################################################"
echo "###  C: relatives retained, relatedness modelled (REGENIE)"
echo "############################################################"
echo ""
echo "Step 1: whole-genome ridge regression"
echo ""

$REGENIE_CMD --step 1 \
  --bed "$IN_UNPRUNED" --keep "$KEEP_UNPRUNED" \
  --extract "${TMP_DIR}/prune_unpruned.prune.in" \
  --phenoFile "$PHENO" $CC_FLAG \
  --covarFile "$COVAR_UNPRUNED" --covarColList "SEX,AGE,${PC_LIST}" \
  --bt --bsize "$BSIZE1" \
  --lowmem --lowmem-prefix "${SCRATCH}/tmp" \
  --threads "$NCORE" \
  --out "${SCRATCH}/pdac_demo_regenie_C_step1" \
  || fail "REGENIE step 1 failed. At this sample size the cross-validated ridge
  can be unstable; read ${SCRATCH}/pdac_demo_regenie_C_step1.log. Consider --loocv."

PRED_LIST="${SCRATCH}/pdac_demo_regenie_C_step1_pred.list"
[ -s "$PRED_LIST" ] || fail "step 1 produced no predictor list: $PRED_LIST"

# The step 1 R-squared says how much the offset is carrying. On a phenotype
# with no polygenic background it is the number that tells you whether C is a
# fair test of the method.
STEP1_RSQ=$(grep -A8 "^phenotype 1" "${SCRATCH}/pdac_demo_regenie_C_step1.log" 2>/dev/null \
            | grep "min value" | sed 's/.*Rsq = \([0-9.]*\).*/\1/' | head -1)
[ -n "${STEP1_RSQ:-}" ] && echo "" && echo "  step 1 cross-validated Rsq: ${STEP1_RSQ}"

echo ""
echo "Step 2: association with the step 1 predictor as offset"
echo ""

$REGENIE_CMD --step 2 \
  --bed "$IN_UNPRUNED" --keep "$KEEP_UNPRUNED" --extract "$VARIANTS" \
  --phenoFile "$PHENO" $CC_FLAG \
  --covarFile "$COVAR_UNPRUNED" --covarColList "SEX,AGE,${PC_LIST}" \
  --bt --firth --approx --pThresh 0.01 \
  --pred "$PRED_LIST" --bsize "$BSIZE2" --threads "$NCORE" \
  --out "${SCRATCH}/pdac_demo_regenie_C_step2"

cp "${SCRATCH}"/pdac_demo_regenie_C_step1.log "$OUT_DIR/" 2>/dev/null || true
cp "${SCRATCH}"/pdac_demo_regenie_C_step2* "$OUT_DIR/" 2>/dev/null || true

# ============================================================================
# D. REGENIE without the offset
#
# Separates the software from the method. Any difference between B and D is
# plink2 against REGENIE; any difference between D and C is the offset.
# ============================================================================
if [ "$RUN_D" = "1" ]; then
  echo ""
  echo "############################################################"
  echo "###  D: REGENIE without the step 1 offset"
  echo "############################################################"
  echo ""

  $REGENIE_CMD --step 2 \
    --bed "$IN_UNPRUNED" --keep "$KEEP_UNPRUNED" --extract "$VARIANTS" \
    --phenoFile "$PHENO" $CC_FLAG \
    --covarFile "$COVAR_UNPRUNED" --covarColList "SEX,AGE,${PC_LIST}" \
    --bt --firth --approx --pThresh 0.01 \
    --ignore-pred --bsize "$BSIZE2" --threads "$NCORE" \
    --out "${SCRATCH}/pdac_demo_regenie_D_nopred"

  cp "${SCRATCH}"/pdac_demo_regenie_D_nopred* "$OUT_DIR/" 2>/dev/null || true
fi

# ============================================================================
# Comparison
#
# Case and control counts are computed from the keep lists and the phenotype
# file rather than scraped from logs, because plink2 and REGENIE report them in
# different formats and a parser that fails silently gives NA in a table that
# looks complete.
# ============================================================================
echo ""
echo "=== Comparison ==="
echo ""

Rscript --vanilla -e "
out_dir <- '${OUT_DIR}'
pheno   <- read.table('${PHENO}', header = TRUE, comment.char = '', check.names = FALSE)
names(pheno)[1] <- sub('^#', '', names(pheno)[1])
ph <- setNames(pheno[[3]], as.character(pheno[[2]]))

count_set <- function(keep_file) {
  k <- as.character(read.table(keep_file, stringsAsFactors = FALSE)[[2]])
  v <- ph[k]
  c(cases = sum(v == 2, na.rm = TRUE), controls = sum(v == 1, na.rm = TRUE))
}
neff <- function(ca, co) if (ca + co > 0) round(4 * ca * co / (ca + co), 1) else NA_real_

read_plink <- function(prefix) {
  f <- Sys.glob(paste0(prefix, '*.glm.logistic.hybrid'))
  if (!length(f)) f <- Sys.glob(paste0(prefix, '*.glm.logistic'))
  if (!length(f)) return(NULL)
  d <- read.table(f[1], header = TRUE, sep = '\t', comment.char = '', check.names = FALSE)
  names(d)[1] <- sub('^#', '', names(d)[1])
  data.frame(ID = d\$ID, P = suppressWarnings(as.numeric(d\$P)), stringsAsFactors = FALSE)
}
read_regenie <- function(prefix) {
  f <- Sys.glob(paste0(prefix, '*.regenie'))
  if (!length(f)) return(NULL)
  d <- read.table(f[1], header = TRUE, comment.char = '#', check.names = FALSE)
  # REGENIE reports -log10(P), not P
  data.frame(ID = d\$ID, P = 10^(-suppressWarnings(as.numeric(d\$LOG10P))),
             stringsAsFactors = FALSE)
}
lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1, lower.tail = FALSE)
}

CAUSAL <- '9:133273682:A:T'
n_pruned   <- count_set('${KEEP_PRUNED}')
n_unpruned <- count_set('${KEEP_UNPRUNED}')

sets <- list(
  list(l = 'A pruned, plink2',      d = read_plink(file.path(out_dir, 'pdac_demo_regenie_A_pruned')),   n = n_pruned),
  list(l = 'B unpruned, plink2',    d = read_plink(file.path(out_dir, 'pdac_demo_regenie_B_unpruned')), n = n_unpruned),
  list(l = 'C unpruned, REGENIE',   d = read_regenie(file.path(out_dir, 'pdac_demo_regenie_C_step2')),  n = n_unpruned))
if (${RUN_D} == 1)
  sets[[4]] <- list(l = 'D unpruned, REGENIE no offset',
                    d = read_regenie(file.path(out_dir, 'pdac_demo_regenie_D_nopred')), n = n_unpruned)

rows <- lapply(sets, function(s) {
  if (is.null(s\$d)) return(NULL)
  p <- s\$d\$P; ok <- is.finite(p)
  data.frame(analysis = s\$l, cases = s\$n[['cases']], controls = s\$n[['controls']],
             neff = neff(s\$n[['cases']], s\$n[['controls']]),
             variants = nrow(s\$d), lambda = round(lam(p), 4),
             sig_5e8 = sum(ok & p < 5e-8), sugg_1e5 = sum(ok & p < 1e-5),
             P_causal = signif(s\$d\$P[match(CAUSAL, s\$d\$ID)], 3),
             stringsAsFactors = FALSE)
})
cmp <- do.call(rbind, Filter(Negate(is.null), rows))
write.table(cmp, file.path(out_dir, 'pdac_demo_regenie_comparison.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)
print(cmp, row.names = FALSE)

cat('\n')
cat('How to read this table.\n')
cat('  A vs B  the whole pruning decision. The sets differ in size and in\n')
cat('          balance, so read neff, not the headline counts. If B is no more\n')
cat('          inflated than A, the retained relatedness cost nothing and the\n')
cat('          pruning removed information for no statistical gain.\n')
cat('  B vs C  the mixed model, on the same people. C lower means the offset\n')
cat('          absorbed real structure; C higher means it absorbed signal.\n')
cat('  D       REGENIE without the offset. B vs D isolates the software from\n')
cat('          the method; D vs C isolates the offset.\n')
cat('  P_causal is the simulated ABO variant. Watch it alongside lambda: an\n')
cat('          analysis that lowers inflation by discarding real signal is not\n')
cat('          an improvement.\n')
" | tee "${OUT_DIR}/pdac_demo_regenie_summary.txt"

echo ""
echo "Written to ${OUT_DIR}/"
echo ""
echo "The demonstration phenotype has no polygenic background, so the step 1"
echo "predictor has little to fit and may absorb the causal variants instead."
echo "Read the step 1 Rsq above before drawing conclusions about C."
echo ""