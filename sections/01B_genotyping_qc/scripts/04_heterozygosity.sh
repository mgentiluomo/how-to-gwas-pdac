#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 04: Heterozygosity outliers, WITHIN ANCESTRY
#
# WHAT CHANGED AND WHY
#   The previous version computed F on the pooled cohort and flagged samples
#   outside mean(F) +/- 3 SD. In a multi-ancestry sample that does not work.
#   The ancestry groups differ in mean F, so the pooled SD is inflated by the
#   variance BETWEEN groups and the interval becomes too wide to flag anything.
#
#   Measured on pdac_demo:
#     pooled          mean F = 0.0683, SD = 0.0288  ->  2 outliers
#     within ancestry mean F = 0.005 to 0.009       -> 20 outliers
#
#   The pooled F was also biased upward (0.068 instead of ~0), because expected
#   heterozygosity was computed from pooled allele frequencies. Running --het
#   separately per group fixes both the centre and the spread.
#
#   Reference: Peterson RE et al. Cell 2019;179:589-603. PMID:31607513
#
# INPUT   pdac_demo_03_filt.bed/bim/fam, demo_data/sample_ancestry.tsv
# OUTPUT  pdac_demo_04_het_<group>.het, pdac_demo_04_het_outliers.txt,
#         pdac_demo_04_het_summary.tsv, pdac_demo_04_filt.bed/bim/fam
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
DATASET_INPUT="${2:-results/qc}"
DATASET_NAME="pdac_demo"
ANCESTRY_FILE="${3:-demo_data/sample_ancestry.tsv}"
OUT_DIR="${4:-results/qc}"
SD_MULT="${5:-3}"

mkdir -p "$OUT_DIR"

if [ ! -s "$ANCESTRY_FILE" ]; then
  echo "ERROR: ancestry file not found: $ANCESTRY_FILE" >&2
  echo "  Step 04 is ancestry-aware. In a real study the groups come from a" >&2
  echo "  preliminary PCA; here they ship with the demonstration data." >&2
  exit 1
fi

GROUPS=$(awk '{print $2}' "$ANCESTRY_FILE" | sort -u | grep -v '^$')
echo "=== Heterozygosity within ancestry groups: $(echo "$GROUPS" | tr '\n' ' ')"

: > "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt"
printf "group\tn\tmean_F\tsd_F\tlower\tupper\toutliers\n" \
  > "${OUT_DIR}/${DATASET_NAME}_04_het_summary.tsv"

for G in $GROUPS; do
  awk -v g="$G" '$2 == g { print $1"\t"$1 }' "$ANCESTRY_FILE" \
    > "${OUT_DIR}/${DATASET_NAME}_04_keep_${G}.txt"

  # --het inside the group, so BOTH observed and expected heterozygosity use
  # that group's own allele frequencies. Using pooled frequencies would leave
  # the Wahlund inflation in F even after stratifying the threshold.
  plink2 \
    --bfile "${DATASET_INPUT}/${DATASET_NAME}_03_filt" \
    --keep "${OUT_DIR}/${DATASET_NAME}_04_keep_${G}.txt" \
    --het \
    --out "${OUT_DIR}/${DATASET_NAME}_04_het_${G}"

  awk -v g="$G" -v k="$SD_MULT" -v out="${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt" '
    NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "F") fi = i; if ($i == "IID") ii = i } ; next }
    { n++; id[n] = $ii; f[n] = $fi + 0; s += $fi; ss += ($fi)^2 }
    END {
      m = s / n
      sd = sqrt((ss - n * m * m) / (n - 1))
      lo = m - k * sd; hi = m + k * sd
      c = 0
      for (i = 1; i <= n; i++)
        if (f[i] < lo || f[i] > hi) { print id[i]"\t"id[i] >> out; c++ }
      printf "%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n", g, n, m, sd, lo, hi, c
    }
  ' "${OUT_DIR}/${DATASET_NAME}_04_het_${G}.het" \
    >> "${OUT_DIR}/${DATASET_NAME}_04_het_summary.tsv"
done

column -t "${OUT_DIR}/${DATASET_NAME}_04_het_summary.tsv"
N_OUT=$(sort -u "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt" | grep -c . || true)
echo ""
echo "Heterozygosity outliers, union across groups: ${N_OUT}"
echo "A mean F close to zero in every group is the check that this worked."

sort -u "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt" \
  > "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.uniq.txt"

plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_03_filt" \
  --remove "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.uniq.txt" \
  --make-bed \
  --out "${OUT_DIR}/${DATASET_NAME}_04_filt"

echo ""
echo "[NEXT] bash scripts/01B_genotyping_qc/05_variant_callrate.sh"
