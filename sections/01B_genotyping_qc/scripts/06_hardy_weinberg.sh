#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 06: Hardy-Weinberg, DIAGNOSTIC ONLY
#
# WHAT CHANGED AND WHY
#   This step no longer removes variants. It measures, and passes everything on.
#
#   The HWE exact test assumes one randomly mating population. Applied to a
#   cohort that spans three continental ancestry groups it detects the Wahlund
#   effect, that is, the allele-frequency differences between the constituent
#   populations, and reports them as genotyping failure.
#
#   Measured on pdac_demo, in controls:
#     pooled, 1e-6    14,378 variants
#     pooled, 1e-10    2,496
#     pooled, 1e-15      199
#     pooled, 1e-20       17   <- and these 17 have ZERO overlap with the
#                                 variants that fail within ancestry groups
#     within ancestry  EUR 4, AFR 0, EAS 1  ->  5 in total
#
#   Lowering the threshold does not repair the test, because the pooled and the
#   stratified test are not lenient and stringent versions of one procedure:
#   they measure different things. The exclusion therefore belongs after the
#   ancestry split, and is applied by
#     scripts/02_population_stratification/06_hwe_within_ancestry.sh
#
#   Reference: Peterson RE et al. Cell 2019;179:589-603. PMID:31607513
#
# INPUT   pdac_demo_05_filt.bed/bim/fam, demo_data/phenotype.txt
# OUTPUT  pdac_demo_06_hwe.hardy, pdac_demo_06_hwe_threshold_sweep.tsv,
#         pdac_demo_06_hwe_flagged.txt (reported, NOT removed),
#         pdac_demo_06_filt.* (a copy of the input, so step numbering holds)
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
PHENOTYPE_FILE="${3:-demo_data/phenotype.txt}"
OUT_DIR="${4:-results/qc}"

mkdir -p "$OUT_DIR"

# Controls only: a true disease association legitimately deviates in cases.
awk 'NR > 1 && ($3 == 1 || $3 == "1") { print $1"\t"$2 }' "$PHENOTYPE_FILE" \
  > "${OUT_DIR}/${DATASET_NAME}_06_controls.txt"
echo "Control samples: $(wc -l < "${OUT_DIR}/${DATASET_NAME}_06_controls.txt")"

plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_05_filt" \
  --keep "${OUT_DIR}/${DATASET_NAME}_06_controls.txt" \
  --hardy \
  --out "${OUT_DIR}/${DATASET_NAME}_06_hwe"

# Threshold sweep: shows that no pooled threshold separates structure from error
printf "threshold\tvariants_below\n" > "${OUT_DIR}/${DATASET_NAME}_06_hwe_threshold_sweep.tsv"
for T in 1e-6 1e-10 1e-15 1e-20 1e-30; do
  N=$(awk -v t="$T" 'NR > 1 && $NF + 0 < t' "${OUT_DIR}/${DATASET_NAME}_06_hwe.hardy" | wc -l)
  printf "%s\t%d\n" "$T" "$N" >> "${OUT_DIR}/${DATASET_NAME}_06_hwe_threshold_sweep.tsv"
done
column -t "${OUT_DIR}/${DATASET_NAME}_06_hwe_threshold_sweep.tsv"

awk 'NR > 1 && $NF + 0 < 1e-6 { print $2 }' "${OUT_DIR}/${DATASET_NAME}_06_hwe.hardy" \
  > "${OUT_DIR}/${DATASET_NAME}_06_hwe_flagged.txt"

echo ""
echo "Flagged at 1e-6 in the pooled control series: $(wc -l < "${OUT_DIR}/${DATASET_NAME}_06_hwe_flagged.txt")"
echo "NONE of these are removed here. Compare this figure with the"
echo "within-ancestry count produced after the ancestry split."

# Carry the dataset forward unchanged, so downstream step numbering is unaffected
for EXT in bed bim fam; do
  cp "${DATASET_INPUT}/${DATASET_NAME}_05_filt.${EXT}" \
     "${OUT_DIR}/${DATASET_NAME}_06_filt.${EXT}"
done

echo ""
echo "[NEXT] bash scripts/01B_genotyping_qc/07_relatedness.sh"
