#!/usr/bin/env bash
################################################################################
# run_experiments.sh
#
# Runs the experiments that sit on top of the canonical run and are not part
# of it: Section 8 (replicates of the simulated phenotype, steps 01 to 10) and
# the Section 4B permutation null with ten random covariates (step 06).
#
# They read the outputs of the canonical run in results/, so run
#   bash scripts/dev/run_all.sh v0.3-data
# first. Nothing here changes a file the canonical run wrote, and
# RUN_MANIFEST.txt is left alone: this script writes its own manifest,
# EXPERIMENTS_MANIFEST.txt, with the commit, the canonical run it sat on, the
# seeds, the software versions and the MD5 of every result file.
#
# USAGE
#   bash scripts/dev/run_experiments.sh                  run everything (hours)
#   bash scripts/dev/run_experiments.sh --manifest-only  no step is run: the
#        manifest is written for the result files already in the section
#        folders, and says so
#
# A full run refuses a dirty working tree, as run_all.sh does.
#
# SEEDS (fixed in the scripts; listed in the manifest)
#   08 step 01  replicate i: 20261124 + i, 200 replicates (GWAS_08_N_REP)
#   08 step 03  bootstrap of the winner's curse: 20261124, 10,000 resamples
#   08 step 05  draw of the 2,000 null variants: 20261124
#   08 step 10  simulation: 20261124 + 10,000,000 x orientation + 100,000 x OR index + i
#   4B step 06  permutations 20261122 (those of the canonical 4B step 04),
#               random covariates 20261125
#
# SETTINGS (environment)
#   GWAS_08_PARALLEL, GWAS_08_THREADS   Section 8 parallelism (script defaults)
#   GWAS_THREADS                        cores of 4B step 06 (default: all)
#   GWAS_PRIVATE_ID_MAP, GWAS_SOURCE_META   only for 08 step 09 (source
#        populations of the analysis sets); they point to files that are not
#        distributed. Without them step 09 writes nothing and the committed
#        table is kept.
################################################################################
set -euo pipefail
cd "$(dirname "$0")/../.."
MODE="run"; [ "${1:-}" = "--manifest-only" ] && MODE="manifest-only"
die() { echo "ERROR: $*" >&2; exit 1; }
say() { echo "$*"; }

[ -f RUN_MANIFEST.txt ] || die "RUN_MANIFEST.txt not found: run scripts/dev/run_all.sh first"
grep -q "NOT A CANONICAL RUN" RUN_MANIFEST.txt && die "results/ holds a non-canonical run"
for f in results/pca/pdac_demo_02_hwe_filt.bed results/assoc/pdac_demo_04A_npc_chosen.txt \
         results/assoc/pdac_demo_04B_perm_lambda.tsv results/meta/eur/keep.txt; do
  [ -f "$f" ] || die "missing output of the canonical run: $f"
done
command -v plink2 >/dev/null || die "plink2 not on PATH"
command -v Rscript >/dev/null || die "Rscript not on PATH"

S8="sections/08_replicates/scripts"
STEPS=(
  "$S8/01_simulate_replicates.R"
  "$S8/02_analyse_replicates.R"
  "$S8/03_summarise.R"
  "$S8/04_expected_information.R"
  "$S8/05_test_statistics.R"
  "$S8/06_components.R"
  "$S8/07_causal_variants_and_structure.R"
  "$S8/08_expected_gain.R"
  "$S8/09_set_composition.R"
  "$S8/10_simulated_detectable_or.R"
  "sections/04B_association_survival/06_random_covariate_null.R"
)
STARTED="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
DIRTY="clean"; [ -n "$(git status --porcelain --untracked-files=no)" ] && DIRTY="dirty (uncommitted changes to tracked files)"

if [ "$MODE" = "run" ]; then
  [ "$DIRTY" = "clean" ] || die "working tree is dirty: commit or stash first"
  LOG_DIR="logs/experiments"; mkdir -p "$LOG_DIR"
  n=0
  for path in "${STEPS[@]}"; do
    n=$((n + 1)); [ -f "$path" ] || die "step missing: $path"
    log="$LOG_DIR/$(printf '%02d' "$n")_$(basename "$path" .R).log"
    printf '[%2d/%d] %s ... ' "$n" "${#STEPS[@]}" "$path"
    Rscript "$path" > "$log" 2>&1 || { echo "FAILED"; tail -20 "$log"; die "step $n failed, see $log"; }
    echo "ok"
  done
  # 4B step 06 writes beside the canonical 4B files in results/assoc
  cp -f results/assoc/pdac_demo_04B_randcov_* results/assoc/pdac_demo_04B_ties.tsv sections/04B_association_survival/results/
fi

FINISHED="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
can() { awk -v k="$1" '$1 == k { $1 = ""; sub(/^ +/, ""); print; exit }' RUN_MANIFEST.txt; }
{
  echo "Experiments manifest (Section 8 and 4B step 06; not the canonical run)"
  echo "======================================================================"
  echo
  if [ "$MODE" = "run" ]; then
    echo "mode           full run: every step below was run by this invocation"
    echo "started        $STARTED"
  else
    echo "mode           manifest only: no step was run by this invocation; the"
    echo "               checksums describe the result files found in the section folders"
  fi
  echo "written        $FINISHED"
  echo "host           $(hostname)"
  echo "commit         $(git rev-parse HEAD)"
  echo "working tree   $DIRTY"
  echo "branch         $(git rev-parse --abbrev-ref HEAD)"
  echo
  echo "canonical run read from results/"
  echo "  commit       $(can commit)"
  echo "  finished     $(can finished)"
  echo "  data release $(awk '/^data release/ { print $3 }' RUN_MANIFEST.txt)"
  echo "  4A components chosen  $(awk '$1 == "chosen_n_pcs" { print $2 }' results/assoc/pdac_demo_04A_npc_chosen.txt)"
  echo
  echo "steps"
  for path in "${STEPS[@]}"; do echo "  $path"; done
  echo
  echo "seeds and sizes"
  echo "  08 step 01   replicate i: 20261124 + i; replicates ${GWAS_08_N_REP:-200} (default 200)"
  echo "  08 step 03   bootstrap seed 20261124, 10000 resamples"
  echo "  08 step 05   null variants: seed 20261124, 2000 variants"
  echo "  08 step 10   seed 20261124, 2000 simulations per OR and orientation"
  echo "  4B step 06   permutations: seed ${GWAS_4B_PERM_SEED:-20261122}, n ${GWAS_4B_N_PERM:-100}; random covariates: seed ${GWAS_4B_RANDCOV_SEED:-20261125}"
  echo "  08 step 09   private inputs $( [ -n "${GWAS_PRIVATE_ID_MAP:-}" ] && [ -n "${GWAS_SOURCE_META:-}" ] && echo "set" || echo "not set: table as committed" )"
  echo
  plink2 --version
  R --version | head -1
  echo "logistf $(Rscript -e 'cat(tryCatch(as.character(packageVersion("logistf")), error = function(e) "not installed (only used for a check in 08 step 05)"))' 2>/dev/null)"
  echo "survival $(Rscript -e 'cat(as.character(packageVersion("survival")))' 2>/dev/null)"
  echo
  echo "result checksums (MD5)"
  for f in sections/08_replicates/results/pdac_demo_08_* sections/04B_association_survival/results/pdac_demo_04B_randcov_* \
           sections/04B_association_survival/results/pdac_demo_04B_ties.tsv; do
    [ -f "$f" ] && printf '  %s  %s\n' "$(md5sum "$f" | cut -d' ' -f1)" "$f"
  done
} > EXPERIMENTS_MANIFEST.txt
say "wrote EXPERIMENTS_MANIFEST.txt ($MODE)"
