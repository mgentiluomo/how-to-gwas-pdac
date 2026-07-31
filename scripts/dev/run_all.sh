#!/usr/bin/env bash

###############################################################################
# run_all.sh
#
# Executes the complete demonstration pipeline in dependency order and records
# enough provenance that the run can be identified, audited and repeated.
#
# The guide asks its readers to make their analyses reproducible. This script
# is how the guide does the same for its own results: every figure on the site
# and in the manuscript comes from one execution of this file, at a recorded
# commit, against a recorded data release, with recorded tool versions.
#
# USAGE
#   bash scripts/dev/run_all.sh [DATA_TAG]
#
#   DATA_TAG   release tag to fetch the demonstration data from.
#              Default: v0.3-data
#
# WHAT IT DOES
#   1. Refuses to run on a dirty working tree, so the recorded commit is true.
#   2. Downloads the demonstration data from the release and verifies every
#      checksum against the released MANIFEST.tsv.
#   3. Captures the version of every tool it is about to use.
#   4. Runs the 33 analysis scripts in dependency order, one log each,
#      stopping at the first failure.
#   5. Copies the committable artefacts into the section results folders.
#   6. Writes RUN_MANIFEST.txt.
#
# WHAT IT DELIBERATELY DOES NOT DO
#   It does not claim that re-running reproduces identical bytes. It does not,
#   and for eigenvector files it cannot, because the sign of an eigenvector is
#   arbitrary. What is reproducible is the analysis and its reported figures.
###############################################################################

set -euo pipefail

DATA_TAG="${1:-v0.3-data}"
REPO="mgentiluomo/how-to-gwas-pdac"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

die_early() { printf '\nFAILED: %s\n' "$*" >&2; exit 1; }

# Locally installed tools take precedence, matching scripts/dev/test.sh.
[ -d tools/bin ] && PATH="$(cd tools/bin && pwd):$PATH" && export PATH

# ---------------------------------------------------------------- 1. commit --
# A result is only attributable to a commit if the tree matches that commit.
if ! command -v git >/dev/null 2>&1; then
  die_early "git not found. The run must be attributable to a commit."
fi
# env/software_versions.md is written by this script, so a previous run always
# leaves it modified. It is excluded here: the runner owns that file.
DIRTY="$(git status --porcelain | grep -v ' env/software_versions.md$' || true)"
if [ -n "$DIRTY" ]; then
  printf '%s\n' "$DIRTY"
  die_early "the working tree has uncommitted changes. Commit or stash them first, \
otherwise the recorded commit does not describe what actually ran."
fi
COMMIT="$(git rev-parse HEAD)"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
printf "Commit  %s\nBranch  %s\n\n" "$COMMIT" "$BRANCH"


LOG_DIR="logs"
MANIFEST="RUN_MANIFEST.txt"
STARTED="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"

mkdir -p "$LOG_DIR" demo_data
MASTER_LOG="$LOG_DIR/00_run_all.log"
: > "$MASTER_LOG"

say() { printf '%s\n' "$*" | tee -a "$MASTER_LOG"; }
die() { printf '\nFAILED: %s\n' "$*" | tee -a "$MASTER_LOG" >&2; exit 1; }

say "==============================================================="
say " Canonical run"
say " started        $STARTED"
say " host           $(hostname)"
say " user           $(whoami)"
say " project root   $PROJECT_ROOT"
say " data release   $DATA_TAG"
say "==============================================================="


# ------------------------------------------------------------------ 2. data --
say ""
say "--- Data ---"
BASE_URL="https://github.com/${REPO}/releases/download/${DATA_TAG}"
FILES=(pdac_demo.bed pdac_demo.bim pdac_demo.fam sample_ancestry.tsv
       phenotype.txt covariates.txt survival.txt truth.tsv MANIFEST.tsv)

for f in "${FILES[@]}"; do
  if [ -s "demo_data/$f" ]; then
    say "  present  $f"
  else
    say "  fetching $f"
    curl -fSL --retry 3 -o "demo_data/$f" "$BASE_URL/$f" \
      || die "could not download $f from release $DATA_TAG"
  fi
done

# Verify what the release itself claims. A silently truncated download is the
# kind of failure that produces plausible wrong numbers rather than an error.
say ""
say "  verifying checksums against the released MANIFEST.tsv"
BAD=0
while IFS=$'\t' read -r file md5 bytes; do
  [ "$file" = "file" ] && continue
  [ -z "${file:-}" ] && continue
  got="$(md5sum "demo_data/$file" | cut -d' ' -f1)"
  if [ "$got" = "$md5" ]; then
    say "    ok        $file"
  else
    say "    MISMATCH  $file  expected $md5  got $got"
    BAD=1
  fi
done < demo_data/MANIFEST.tsv
[ "$BAD" -eq 0 ] || die "checksum mismatch in the downloaded data"

# The genotypes are not covered by the manifest; record their checksums so this
# run can be compared with any other.
say ""
say "  genotype checksums, recorded for comparison"
for f in pdac_demo.bed pdac_demo.bim pdac_demo.fam sample_ancestry.tsv; do
  say "    $(md5sum "demo_data/$f" | cut -d' ' -f1)  $f"
done

# -------------------------------------------------------------- 3. versions --
say ""
say "--- Tool versions ---"
command -v plink2 >/dev/null 2>&1 || die "plink2 not found on PATH"
command -v Rscript >/dev/null 2>&1 || die "Rscript not found on PATH"
# The relatedness step uses PLINK 1.9 for the PI_HAT comparison. Some systems
# install it under another name; symlink it into tools/bin, which is on PATH.
command -v plink >/dev/null 2>&1 || die "PLINK 1.9 not found on PATH as 'plink'. \
The relatedness step needs it. If it is installed under another name, run: \
mkdir -p tools/bin && ln -sf \"\$(command -v plink19)\" tools/bin/plink"

PLINK_V="$(plink2 --version 2>&1 | head -1)"
PLINK1_V="$(plink --version 2>&1 | head -1)"
R_V="$(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
OS_V="$(. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME" || uname -sr)"
say "  $PLINK_V"
say "  $R_V"
say "  $OS_V"

PKG_TABLE="$(Rscript -e '
pkgs <- readLines("env/r_packages.txt")
pkgs <- trimws(pkgs); pkgs <- pkgs[nzchar(pkgs)]
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  cat(sprintf("| %s | %s |\n", p, v))
}' 2>/dev/null)"
say ""
say "  R packages:"
printf '%s\n' "$PKG_TABLE" | sed 's/^/    /' | tee -a "$MASTER_LOG"

if printf '%s' "$PKG_TABLE" | grep -q "NOT INSTALLED"; then
  die "an R package listed in env/r_packages.txt is not installed"
fi

# env/software_versions.md is a template in the repository and the manuscript
# promises exact versions. Write it from what actually ran.
cat > env/software_versions.md <<EOF
# Software versions

Written automatically by \`scripts/dev/run_all.sh\`. These are the versions that
produced the published results; do not edit by hand.

| Item | Value |
|---|---|
| Run started | $STARTED |
| Commit | \`$COMMIT\` |
| Data release | \`$DATA_TAG\` |
| Operating system | $OS_V |
| PLINK 2 | $PLINK_V |
| PLINK 1.9 | $PLINK1_V |
| R | $R_V |

## R packages

| Package | Version |
|---|---|
$PKG_TABLE
EOF
say ""
say "  wrote env/software_versions.md"

# ------------------------------------------------------------------ 4. run --
# Dependency order. Entries may carry arguments; most scripts need none because
# their defaults already chain together. Two R scripts are absent on purpose:
# 01_initial_qc_plots.R and 09_qc_report.R are invoked by their own shell steps.
# Section 1A is excluded on purpose: it demonstrates raw
# array output handling on its own synthetic example files and does not depend
# on the phenotype, and one of its steps writes into demo_data/.
STEPS=(
  "01B_genotyping_qc/scripts/01_initial_qc_stats.sh"
  "01B_genotyping_qc/scripts/02_sample_callrate.sh"
  "01B_genotyping_qc/scripts/03_sex_check.sh"
  "01B_genotyping_qc/scripts/04_heterozygosity.sh"
  "01B_genotyping_qc/scripts/05_variant_callrate.sh"
  "01B_genotyping_qc/scripts/06_hardy_weinberg.sh"
  "01B_genotyping_qc/scripts/06_hwe_plot.R results/qc/pdac_demo_06_hwe.hardy results/qc/pdac_demo_06_hwe"
  "01B_genotyping_qc/scripts/07_relatedness.sh"
  "01B_genotyping_qc/scripts/08_maf_filter.sh"
  "01B_genotyping_qc/scripts/09_qc_summary.sh"

  "02_population_stratification/scripts/01_ld_pruning.sh"
  "02_population_stratification/scripts/02_pca_all.sh"
  "02_population_stratification/scripts/03_pca_plots.R"
  "02_population_stratification/scripts/04_define_analysis_set.sh"
  "02_population_stratification/scripts/05_pca_within_eur.sh"
  "02_population_stratification/scripts/06_hwe_within_ancestry.sh"

  "03B_statistical_power/scripts/01_power_curves.R"
  "03B_statistical_power/scripts/02_power_events.R"

  "04A_association_binary/scripts/01_make_covariates.R"
  "04A_association_binary/scripts/02_association.sh"
  "04A_association_binary/scripts/03_qq_lambda.R"
  "04A_association_binary/scripts/04_manhattan.R"
  "04A_association_binary/scripts/05_summary.R"

  "05_meta_analysis/scripts/01_stratified_gwas.sh"
  "05_meta_analysis/scripts/02_harmonise.R"
  "05_meta_analysis/scripts/03_meta_analysis.R"
  "05_meta_analysis/scripts/04_meta_plots.R"

  "06_finemapping_annotation/scripts/01_region_and_ld.sh"
  "06_finemapping_annotation/scripts/02_finemap_abf.R"
  "06_finemapping_annotation/scripts/03_finemap_plots.R"

  "07_linkage_disequilibrium/scripts/01_ld_decay.sh"
)

say ""
say "--- Pipeline, ${#STEPS[@]} steps ---"
n=0
for entry in "${STEPS[@]}"; do
  n=$((n + 1))
  read -r rel args <<< "$entry"
  path="sections/$rel"
  [ -f "$path" ] || die "step $n missing from the repository: $path"

  tag="$(printf '%02d' "$n")_$(basename "$rel" | tr '/' '_')"
  log="$LOG_DIR/${tag}.log"
  printf '  [%2d/%d] %-58s ' "$n" "${#STEPS[@]}" "$rel" | tee -a "$MASTER_LOG"

  t0=$(date +%s)
  if [[ "$path" == *.R ]]; then
    Rscript "$path" $args > "$log" 2>&1 || { echo "FAILED"; tail -20 "$log"; die "step $n failed, see $log"; }
  else
    bash "$path" $args > "$log" 2>&1 || { echo "FAILED"; tail -20 "$log"; die "step $n failed, see $log"; }
  fi
  t1=$(date +%s)
  printf 'ok  %3ds\n' "$((t1 - t0))" | tee -a "$MASTER_LOG"
done

# ----------------------------------------------------------------- 5. sync --
# The scripts write to a flat results tree. The repository commits the small
# artefacts inside each section. .gitignore decides what is small enough, so
# everything is copied and git filters it.
say ""
say "--- Collecting artefacts into the section folders ---"
declare -A MAP=(
  [qc]="01B_genotyping_qc"
  [pca]="02_population_stratification"
  [power]="03B_statistical_power"
  [assoc]="04A_association_binary"
  [meta]="05_meta_analysis"
  [finemap]="06_finemapping_annotation"
  [ld]="07_linkage_disequilibrium"
)
for area in "${!MAP[@]}"; do
  src="results/$area"
  dst="sections/${MAP[$area]}/results"
  if [ -d "$src" ]; then
    mkdir -p "$dst"
    cp -f "$src"/* "$dst"/ 2>/dev/null || true
    say "  $(printf '%-9s' "$area") -> $dst  ($(ls -1 "$src" | wc -l) files)"
  else
    say "  $(printf '%-9s' "$area") -> MISSING, no $src produced"
  fi
done

# -------------------------------------------------------------- 6. manifest --
FINISHED="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
{
  echo "Canonical run manifest"
  echo "======================"
  echo
  echo "started        $STARTED"
  echo "finished       $FINISHED"
  echo "host           $(hostname)"
  echo "user           $(whoami)"
  echo "commit         $COMMIT"
  echo "branch         $BRANCH"
  echo "data release   $DATA_TAG"
  echo "steps          ${#STEPS[@]}"
  echo
  echo "$PLINK_V"
  echo "$PLINK1_V"
  echo "$R_V"
  echo "$OS_V"
  echo
  echo "data checksums"
  for f in "${FILES[@]}"; do
    echo "  $(md5sum "demo_data/$f" | cut -d' ' -f1)  $f"
  done
  echo
  echo "Per-step logs are in $LOG_DIR/."
} > "$MANIFEST"

say ""
say "==============================================================="
say " Completed  $FINISHED"
say ""
say " Wrote  $MANIFEST"
say "        env/software_versions.md"
say "        $LOG_DIR/  (one log per step)"
say ""
say " Next:  git status, review the diff, then commit."
say "==============================================================="
