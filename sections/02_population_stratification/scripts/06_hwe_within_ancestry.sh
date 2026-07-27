#!/usr/bin/env bash

################################################################################
# Section 2: Population structure — Step 06: HWE within ancestry (the exclusion)
#
# This is the Hardy-Weinberg filter that Step 1B-06 deliberately did not apply.
# It runs here because the test needs a genetically homogeneous group, which
# does not exist until ancestry has been assigned.
#
# RULE
#   Test in controls, within each ancestry group, at 1e-6. A variant is excluded
#   if it fails in ANY group: a genotyping error is a property of the assay, so
#   failing anywhere is evidence against the variant, while failing nowhere in a
#   well-powered group is evidence for it. Both the per-group counts and the
#   union are reported, because the choice between "any" and "all" is a decision
#   the reader should see rather than inherit.
#
#   Aggregating per-variant results across homogeneous ancestry subsets is the
#   approach recommended for multi-ethnic studies; see also RUTH (Kwong et al.
#   2021) for cohorts where discrete groups cannot be defined.
#
# ON pdac_demo   EUR 4, AFR 0, EAS 1, union 5
#                against 14,378 from the pooled test at the same threshold
#
# INPUT   results/qc/pdac_demo_08_filt.*, demo_data/phenotype.txt,
#         demo_data/sample_ancestry.tsv
# OUTPUT  pdac_demo_02_hwe_<group>.hardy, pdac_demo_02_hwe_summary.tsv,
#         pdac_demo_02_hwe_exclude.txt, pdac_demo_02_hwe_filt.bed/bim/fam
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

SEED="${1:-2026}"
QC_DIR="${2:-results/qc}"
DATASET_NAME="pdac_demo"
PHENOTYPE_FILE="${3:-demo_data/phenotype.txt}"
ANCESTRY_FILE="${4:-demo_data/sample_ancestry.tsv}"
OUT_DIR="${5:-results/pca}"
HWE_P="${6:-1e-6}"

mkdir -p "$OUT_DIR"

awk 'NR > 1 && ($3 == 1 || $3 == "1") { print $2 }' "$PHENOTYPE_FILE" \
  > "${OUT_DIR}/${DATASET_NAME}_02_controls.ids"

GROUPS=$(awk '{print $2}' "$ANCESTRY_FILE" | sort -u | grep -v '^$')
printf "group\tcontrols_tested\tvariants_failing_%s\n" "$HWE_P" \
  > "${OUT_DIR}/${DATASET_NAME}_02_hwe_summary.tsv"
: > "${OUT_DIR}/${DATASET_NAME}_02_hwe_exclude.raw"

for G in $GROUPS; do
  awk -v g="$G" 'NR == FNR { c[$1]; next } ($2 == g && $1 in c) { print $1"\t"$1 }' \
    "${OUT_DIR}/${DATASET_NAME}_02_controls.ids" "$ANCESTRY_FILE" \
    > "${OUT_DIR}/${DATASET_NAME}_02_ctrl_${G}.txt"

  plink2 \
    --bfile "${QC_DIR}/${DATASET_NAME}_08_filt" \
    --keep "${OUT_DIR}/${DATASET_NAME}_02_ctrl_${G}.txt" \
    --hardy \
    --out "${OUT_DIR}/${DATASET_NAME}_02_hwe_${G}"

  NKEEP=$(grep -oE '^--keep: [0-9]+ samples remaining' \
          "${OUT_DIR}/${DATASET_NAME}_02_hwe_${G}.log" | grep -oE '[0-9]+' | head -1)
  awk -v t="$HWE_P" 'NR > 1 && $NF + 0 < t { print $2 }' \
    "${OUT_DIR}/${DATASET_NAME}_02_hwe_${G}.hardy" \
    > "${OUT_DIR}/${DATASET_NAME}_02_hwe_fail_${G}.txt"

  printf "%s\t%s\t%d\n" "$G" "${NKEEP:-NA}" \
    "$(wc -l < "${OUT_DIR}/${DATASET_NAME}_02_hwe_fail_${G}.txt")" \
    >> "${OUT_DIR}/${DATASET_NAME}_02_hwe_summary.tsv"

  cat "${OUT_DIR}/${DATASET_NAME}_02_hwe_fail_${G}.txt" \
    >> "${OUT_DIR}/${DATASET_NAME}_02_hwe_exclude.raw"
done

sort -u "${OUT_DIR}/${DATASET_NAME}_02_hwe_exclude.raw" \
  > "${OUT_DIR}/${DATASET_NAME}_02_hwe_exclude.txt"
rm -f "${OUT_DIR}/${DATASET_NAME}_02_hwe_exclude.raw"

column -t "${OUT_DIR}/${DATASET_NAME}_02_hwe_summary.tsv"
echo ""
echo "Union across groups, excluded: $(wc -l < "${OUT_DIR}/${DATASET_NAME}_02_hwe_exclude.txt")"
if [ -s "${QC_DIR}/${DATASET_NAME}_06_hwe_flagged.txt" ]; then
  echo "Pooled test flagged at the same threshold: $(wc -l < "${QC_DIR}/${DATASET_NAME}_06_hwe_flagged.txt")"
  echo "The difference between those two numbers is the Wahlund effect."
fi

plink2 \
  --bfile "${QC_DIR}/${DATASET_NAME}_08_filt" \
  --exclude "${OUT_DIR}/${DATASET_NAME}_02_hwe_exclude.txt" \
  --make-bed \
  --out "${OUT_DIR}/${DATASET_NAME}_02_hwe_filt"

echo ""
echo "This is the dataset carried into association testing."
