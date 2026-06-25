#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC - Step 09: QC Summary Report
#
# PURPOSE:
#   Create the final genotyping QC report:
#   - Track sample/variant counts through each filter
#   - Summarize key relatedness diagnostics
#   - Write a text summary and one polished multi-page PDF report
#
# INPUT:
#   - All intermediate datasets and summary tables from Steps 01-08
#
# OUTPUT:
#   - pdac_demo_09_qc_counts.tsv   - structured count table
#   - pdac_demo_09_qc_summary.txt  - plain-text filtering report
#   - pdac_demo_09_qc_report.pdf   - combined report with tables and figures
#
# NOTES:
#   The PDF is rendered by scripts/01B_genotyping_qc/09_qc_report.R.
#
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

# Configuration
SEED="${1:-2026}"
DATASET_NAME="pdac_demo"
DATASET_INPUT="${2:-demo_data}"
OUT_DIR="${3:-results/qc}"

SUMMARY_FILE="${OUT_DIR}/${DATASET_NAME}_09_qc_summary.txt"
COUNTS_FILE="${OUT_DIR}/${DATASET_NAME}_09_qc_counts.tsv"
REPORT_FILE="${OUT_DIR}/${DATASET_NAME}_09_qc_report.pdf"
REPORT_SCRIPT="${SCRIPT_DIR}/09_qc_report.R"
RELATEDNESS_DIAGNOSTICS="${OUT_DIR}/${DATASET_NAME}_07_relatedness_pruning_diagnostics.tsv"

mkdir -p "$OUT_DIR"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript was not found in PATH."
  echo "Run Step 5 again: bash scripts/dev/tools_setup.sh"
  exit 1
fi

if [ ! -f "$REPORT_SCRIPT" ]; then
  echo "ERROR: Report renderer not found: $REPORT_SCRIPT"
  exit 1
fi

# ============================================================================
# Helper functions
# ============================================================================

get_counts() {
  local bfile=$1
  local nsamp="NA"
  local nvar="NA"

  if [ -f "${bfile}.fam" ]; then
    nsamp=$(wc -l < "${bfile}.fam")
  fi
  if [ -f "${bfile}.bim" ]; then
    nvar=$(wc -l < "${bfile}.bim")
  fi

  echo "${nsamp} ${nvar}"
}

count_data_rows() {
  local file=$1
  if [ -f "$file" ]; then
    awk 'END { print (NR > 0 ? NR - 1 : 0) }' "$file"
  else
    echo "NA"
  fi
}

is_integer() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

format_pct() {
  local numerator=$1
  local denominator=$2
  awk -v n="$numerator" -v d="$denominator" 'BEGIN {
    if (d > 0) printf "%.1f", 100 * n / d;
    else printf "NA";
  }'
}

read_relatedness_metric() {
  local metric=$1
  if [ -f "$RELATEDNESS_DIAGNOSTICS" ]; then
    awk -v metric="$metric" '
      $1 == metric {
        print $2
        found = 1
      }
      END {
        if (!found) print "NA"
      }
    ' "$RELATEDNESS_DIAGNOSTICS"
  else
    echo "NA"
  fi
}

read_relatedness_pruning_pairs() {
  if [ -f "$RELATEDNESS_DIAGNOSTICS" ]; then
    awk '
      $1 ~ /^king_pruning_pairs_gt_/ {
        print $2
        found = 1
      }
      END {
        if (!found) print "NA"
      }
    ' "$RELATEDNESS_DIAGNOSTICS"
  else
    echo "NA"
  fi
}

# ============================================================================
# STEP 1: Collect sample and variant counts
# ============================================================================

echo ""
echo "=== Collecting QC metrics ==="
echo ""

read NSAMP_RAW NVAR_RAW <<< "$(get_counts "${DATASET_INPUT}/${DATASET_NAME}")"
read NSAMP_02 NVAR_02 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_02_filt")"
read NSAMP_03 NVAR_03 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_03_filt")"
read NSAMP_04 NVAR_04 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_04_filt")"
read NSAMP_05 NVAR_05 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_05_filt")"
read NSAMP_06 NVAR_06 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_06_filt")"
read NSAMP_07 NVAR_07 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_07_filt")"
read NSAMP_08 NVAR_08 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_08_filt")"

if [ "$NSAMP_RAW" = "NA" ]; then
  NSAMP_RAW=$(count_data_rows "${OUT_DIR}/${DATASET_NAME}_01_qc.smiss")
fi

if [ "$NVAR_RAW" = "NA" ]; then
  NVAR_RAW=$(count_data_rows "${OUT_DIR}/${DATASET_NAME}_01_qc.vmiss")
fi

if [ "$NSAMP_RAW" = "NA" ] && is_integer "$NSAMP_02"; then
  SAMPLE_CALLRATE_REMOVED=$(count_data_rows "${OUT_DIR}/${DATASET_NAME}_02_samples_callrate.mindrem.id")
  if is_integer "$SAMPLE_CALLRATE_REMOVED"; then
    NSAMP_RAW=$((NSAMP_02 + SAMPLE_CALLRATE_REMOVED))
  fi
fi

if [ "$NVAR_RAW" = "NA" ] && is_integer "$NVAR_02"; then
  # Step 02 removes samples only, so the variant count is unchanged.
  NVAR_RAW="$NVAR_02"
fi

for metric_name in NSAMP_RAW NVAR_RAW NSAMP_02 NVAR_02 NSAMP_03 NVAR_03 NSAMP_04 NVAR_04 NSAMP_05 NVAR_05 NSAMP_06 NVAR_06 NSAMP_07 NVAR_07 NSAMP_08 NVAR_08; do
  metric_value="${!metric_name}"
  if ! is_integer "$metric_value"; then
    echo "ERROR: Could not resolve numeric QC count: ${metric_name}=${metric_value}"
    echo "Check DATASET_INPUT (${DATASET_INPUT}) and OUT_DIR (${OUT_DIR})."
    exit 1
  fi
done

echo "Raw: ${NSAMP_RAW} samples, ${NVAR_RAW} variants"

SAMPLE_RETENTION=$(format_pct "$NSAMP_08" "$NSAMP_RAW")
VARIANT_RETENTION=$(format_pct "$NVAR_08" "$NVAR_RAW")
SAMPLE_LOSS_TOTAL=$((NSAMP_RAW - NSAMP_08))
VARIANT_LOSS_TOTAL=$((NVAR_RAW - NVAR_08))
SAMPLE_LOSS_PCT=$(format_pct "$SAMPLE_LOSS_TOTAL" "$NSAMP_RAW")
VARIANT_LOSS_PCT=$(format_pct "$VARIANT_LOSS_TOTAL" "$NVAR_RAW")

REL_PRUNING_PAIRS=$(read_relatedness_pruning_pairs)
REL_UNIQUE_SAMPLES=$(read_relatedness_metric "unique_samples_in_king_pruning_pairs")
REL_COMPONENTS=$(read_relatedness_metric "relatedness_components")
REL_LARGEST_COMPONENT=$(read_relatedness_metric "largest_component_samples")
REL_PLINK2_REMOVED=$(read_relatedness_metric "plink2_king_cutoff_removed_samples")
REL_PHENO_REMOVED=$(read_relatedness_metric "phenotype_aware_removed_samples")
REL_REMOVED_PCT=$(read_relatedness_metric "phenotype_aware_removed_percent")
REL_CONTROLS_REMOVED=$(read_relatedness_metric "controls_removed")
REL_CASES_REMOVED=$(read_relatedness_metric "cases_removed")

# ============================================================================
# STEP 2: Write structured counts and plain-text summary
# ============================================================================

echo ""
echo "=== Writing Step 09 summary inputs ==="
echo ""

cat > "$COUNTS_FILE" << EOF
step_index	step	samples	variants	samples_lost	variants_lost
0	Raw data	${NSAMP_RAW}	${NVAR_RAW}	0	0
1	Initial stats	${NSAMP_RAW}	${NVAR_RAW}	0	0
2	Sample call rate	${NSAMP_02}	${NVAR_02}	$((NSAMP_RAW - NSAMP_02))	0
3	Sex check	${NSAMP_03}	${NVAR_03}	$((NSAMP_02 - NSAMP_03))	0
4	Heterozygosity	${NSAMP_04}	${NVAR_04}	$((NSAMP_03 - NSAMP_04))	0
5	Variant call rate	${NSAMP_05}	${NVAR_05}	0	$((NVAR_04 - NVAR_05))
6	Hardy-Weinberg	${NSAMP_06}	${NVAR_06}	0	$((NVAR_05 - NVAR_06))
7	Relatedness	${NSAMP_07}	${NVAR_07}	$((NSAMP_06 - NSAMP_07))	0
8	MAF filter	${NSAMP_08}	${NVAR_08}	0	$((NVAR_07 - NVAR_08))
EOF

cat > "$SUMMARY_FILE" << EOF
PDAC GWAS Genotyping Quality Control Summary
============================================

Dataset: ${DATASET_NAME}
Analysis date: $(date)
Genome build: GRCh38

QC Pipeline Summary
-------------------

Step                          Samples    Variants   Samples lost   Variants lost
00. Raw data                  ${NSAMP_RAW}      ${NVAR_RAW}     0              0
01. Initial stats             ${NSAMP_RAW}      ${NVAR_RAW}     0              0
02. Sample call rate          ${NSAMP_02}      ${NVAR_02}     $((NSAMP_RAW - NSAMP_02))             0
03. Sex check                 ${NSAMP_03}      ${NVAR_03}     $((NSAMP_02 - NSAMP_03))              0
04. Heterozygosity            ${NSAMP_04}      ${NVAR_04}     $((NSAMP_03 - NSAMP_04))              0
05. Variant call rate         ${NSAMP_05}      ${NVAR_05}     0              $((NVAR_04 - NVAR_05))
06. Hardy-Weinberg            ${NSAMP_06}      ${NVAR_06}     0              $((NVAR_05 - NVAR_06))
07. Relatedness (KING)        ${NSAMP_07}      ${NVAR_07}     $((NSAMP_06 - NSAMP_07))             0
08. MAF filter (>=1%)         ${NSAMP_08}      ${NVAR_08}     0              $((NVAR_07 - NVAR_08))

Final dataset: ${NSAMP_08} samples and ${NVAR_08} variants
Total sample loss: ${SAMPLE_LOSS_TOTAL} (${SAMPLE_LOSS_PCT}%)
Total variant loss: ${VARIANT_LOSS_TOTAL} (${VARIANT_LOSS_PCT}%)
Sample retention: ${SAMPLE_RETENTION}%
Variant retention: ${VARIANT_RETENTION}%

QC Thresholds Applied
---------------------

Sample-level:
- Sample call rate: --mind 0.02
- Sex check: X chromosome F-statistic discordance
- Heterozygosity: autosomal F +/- 3 SD
- Relatedness: KING kinship > 0.1875
- Relatedness pruning strategy: phenotype-aware; prefer removing controls over cases

Variant-level:
- Variant call rate: --geno 0.05
- Hardy-Weinberg equilibrium: --hwe 1e-6 in controls
- Minor allele frequency: --maf 0.01

Relatedness Diagnostics
-----------------------

- KING pairs above pruning threshold: ${REL_PRUNING_PAIRS}
- Unique samples in those pairs: ${REL_UNIQUE_SAMPLES}
- Relatedness graph components: ${REL_COMPONENTS}
- Largest relatedness component: ${REL_LARGEST_COMPONENT} samples
- PLINK2 --king-cutoff removed: ${REL_PLINK2_REMOVED} samples
- Phenotype-aware pruning removed: ${REL_PHENO_REMOVED} samples (${REL_REMOVED_PCT}%)
- Controls removed: ${REL_CONTROLS_REMOVED}
- Cases removed: ${REL_CASES_REMOVED}

Interpretation Notes
--------------------

This tutorial uses a rare-cancer case-control setting. Case samples are valuable,
so relatedness pruning is phenotype-aware: in a related case-control pair, the
control is removed where possible and the case is preserved.

If relatedness removes a large fraction of the dataset, inspect:
- ${DATASET_NAME}_07_pruning_components.tsv
- ${DATASET_NAME}_07_relatedness_pruning_diagnostics.tsv
- ${DATASET_NAME}_07_removal_decisions.tsv

The final dataset is ready for population stratification analysis and PCA.
EOF

echo "Summary table: $COUNTS_FILE"
echo "Text summary:  $SUMMARY_FILE"

# ============================================================================
# STEP 3: Render combined PDF report
# ============================================================================

echo ""
echo "=== Rendering combined QC PDF report ==="
echo ""

Rscript "$REPORT_SCRIPT" \
  "$DATASET_NAME" \
  "$SEED" \
  "$COUNTS_FILE" \
  "$OUT_DIR" \
  "$SUMMARY_FILE" \
  "$REPORT_FILE"

# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo ""
echo "=============================================================="
echo "GENOTYPING QC PIPELINE COMPLETE"
echo "=============================================================="
echo ""
echo "Final clean dataset:"
echo "  - ${OUT_DIR}/${DATASET_NAME}_08_filt.bed/bim/fam"
echo "  - ${NSAMP_08} samples"
echo "  - ${NVAR_08} variants"
echo ""
echo "Step 09 report outputs:"
echo "  - ${COUNTS_FILE}"
echo "  - ${SUMMARY_FILE}"
echo "  - ${REPORT_FILE}"
echo ""
echo "Ready for downstream analysis:"
echo "  - Section 2: Population stratification (PCA)"
echo "  - Section 3: Imputation (optional)"
echo "  - Section 4: Association testing"
echo ""
