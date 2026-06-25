#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 07: Relatedness & Sample Pruning (KING)
# 
# PURPOSE:
#   Identify and remove related samples. Relatedness violates the assumption
#   of independent samples in GWAS, leading to inflated test statistics and
#   false positive associations.
#
# INPUT:
#   - pdac_demo_06_filt.bed/bim/fam (from Step 06, HWE-filtered)
#   - phenotype.txt with columns FID IID PHENO (1=control, 2=case)
#
# OUTPUT:
#   - pdac_demo_07_filt.bed/bim/fam — pruned to unrelated samples
#   - pdac_demo_07_king_cutoff.king.cutoff.in.id — PLINK2 KING keep list
#   - pdac_demo_07_king_cutoff.king.cutoff.out.id — PLINK2 KING removal list
#   - pdac_demo_07_ibd_pi_hat.genome — PLINK1.9 IBD/PI_HAT report
#   - pdac_demo_07_ibd_pi_hat_annotated.tsv — PI_HAT pairs with category labels
#   - pdac_demo_07_ibd_pi_hat_summary.tsv — PI_HAT relationship category counts
#   - pdac_demo_07_related_pairs.kin0 — KING table for related sample pairs
#   - pdac_demo_07_related_pairs_annotated.tsv — related pairs with category labels
#   - pdac_demo_07_relatedness_summary.tsv — relatedness category counts
#   - pdac_demo_07_king_vs_pihat_comparison.tsv — paired KING/PI_HAT comparison
#   - pdac_demo_07_king_vs_pihat_summary.tsv — KING/PI_HAT agreement summary
#   - pdac_demo_07_pruning_components.tsv — relatedness graph component summary
#   - pdac_demo_07_relatedness_pruning_diagnostics.tsv — key pruning diagnostics
#   - pdac_demo_07_removal_decisions.tsv — phenotype-aware pruning decisions
#   - pdac_demo_07_removal_summary.tsv — case/control removal summary
#   - pdac_demo_07_removal_phenotype_counts.tsv — final counts from removal list
#   - pdac_demo_07_removal_list.txt — samples removed for relatedness
#
# THRESHOLD:
#   Removal cutoff: KING kinship coefficient > 0.1875
#     This primarily removes duplicate/twin and first-degree relationships.
#   Reporting cutoff: KING kinship coefficient >= 0.0442
#     This also reports weaker pair categories for teaching/review.
#   PI_HAT reporting cutoff: PI_HAT >= 0.0884
#     This is roughly comparable to KING 0.0442 on ideal pedigree expectations.
#
# METHOD:
#   1. LD-prune variants (keep ~30K for kinship inference; too many cause noise)
#   2. Run PLINK2 --king-cutoff 0.1875 as a standard reference check
#   3. Run PLINK1.9 --genome as a PI_HAT teaching/comparison report
#   4. Create an annotated KING pair table and category summary
#   5. Compare KING kinship and PI_HAT for QC review
#   6. Create a phenotype-aware removal list from KING pairs above 0.1875
#   7. Prefer removing controls over cases in case-control related pairs
#
# NOTES:
#   - PDAC may have related samples (family cohorts)
#   - In rare-cancer GWAS, cases are precious; remove controls first where possible
#   - KING algorithm is fast and accurate for unrelated to 3rd-degree relatives
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
DATASET_INPUT="${2:-results/qc}"
DATASET_NAME="pdac_demo"
OUT_DIR="${3:-results/qc}"
PHENOTYPE_FILE="${4:-demo_data/phenotype.txt}"
REMOVE_THRESHOLD="${5:-0.1875}"
REPORT_THRESHOLD="${6:-0.0442}"
PIHAT_REPORT_THRESHOLD="${7:-0.0884}"

mkdir -p "$OUT_DIR"

if ! command -v plink2 >/dev/null 2>&1; then
  echo "✗ plink2 was not found in PATH."
  echo "  Run Step 5 again: bash scripts/dev/tools_setup.sh"
  exit 1
fi

if ! command -v plink >/dev/null 2>&1; then
  echo "✗ plink 1.9 was not found in PATH."
  echo "  Step 07 uses plink for the PI_HAT report."
  echo "  Run Step 5 again: bash scripts/dev/tools_setup.sh"
  exit 1
fi

if [ ! -f "$PHENOTYPE_FILE" ]; then
  echo "✗ Phenotype file not found: $PHENOTYPE_FILE"
  echo "  Step 07 uses phenotype-aware pruning to preserve rare cancer cases."
  echo "  Expected format: FID IID PHENO, where 1=control and 2=case."
  exit 1
fi

if ! awk '
  function clean(value) {
    value = tolower(value)
    gsub(/\r/, "", value)
    return value
  }
  function known_pheno(value) {
    value = clean(value)
    return value == "1" || value == "2" || value == "control" || value == "case" || value == "unaffected" || value == "affected"
  }
  FNR == 1 {
    first = $1
    second = $2
    third = $3
    gsub(/^#/, "", first)
    gsub(/^#/, "", second)
    gsub(/^#/, "", third)
    if (first == "FID" && (second == "IID" || second == "ID") && third == "PHENO") {
      next
    }
  }
  NF < 3 {
    bad_format++
    next
  }
  {
    records++
    if (!known_pheno($3)) {
      bad_pheno++
    }
  }
  END {
    if (records == 0) {
      print "✗ Phenotype file has no sample rows: FID IID PHENO expected." > "/dev/stderr"
      exit 1
    }
    if (bad_format > 0 || bad_pheno > 0) {
      print "✗ Phenotype file format problem." > "/dev/stderr"
      print "  Expected columns: FID IID PHENO" > "/dev/stderr"
      print "  Expected PHENO coding: 1=control, 2=case" > "/dev/stderr"
      print "  Rows with fewer than 3 columns: " bad_format > "/dev/stderr"
      print "  Rows with unknown PHENO values: " bad_pheno > "/dev/stderr"
      exit 1
    }
  }
' "$PHENOTYPE_FILE"; then
  exit 1
fi

# ============================================================================
# STEP 1: LD-based SNP pruning
# ============================================================================
echo ""
echo "=== LD-pruning variants for kinship inference ==="
echo ""

# Why prune SNPs for kinship?
#   Too many variants make kinship estimates noisy (redundant information).
#   Standard practice: keep ~30K-50K independent variants.
#   Parameters:
#   --indep-pairwise <window_size> <window_step> <r2_threshold>
#   - 50 bp window (50 SNPs at a time)
#   - 5 SNP step (move by 5 SNPs)
#   - r2 > 0.2 (remove if more correlated than this)

plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_06_filt" \
  --indep-pairwise 50 5 0.2 \
  --out "${OUT_DIR}/${DATASET_NAME}_07_prune"

echo "✓ LD-pruning complete. Kept $(wc -l < "${OUT_DIR}/${DATASET_NAME}_07_prune.prune.in") independent variants"

# ============================================================================
# STEP 2: Define relatedness thresholds and run standard KING cutoff
# ============================================================================
echo ""
echo "=== Defining KING relatedness thresholds ==="
echo ""

# Why KING algorithm?
#   KING (Kinship-based INference for Gwas) is robust to population stratification.
#   Computes KING kinship coefficients:
#   - Unrelated: near 0
#   - First-degree relatives: near 0.25
#   - Duplicates/twins: near 0.5
#   - Threshold 0.1875 catches duplicate/twin and first-degree relationships
#
echo "✓ Removal threshold: KING kinship > ${REMOVE_THRESHOLD}"
echo "✓ Reporting threshold: KING kinship >= ${REPORT_THRESHOLD}"
echo "✓ PI_HAT reporting threshold: PI_HAT >= ${PIHAT_REPORT_THRESHOLD}"
echo ""
echo "=== Running standard PLINK2 KING cutoff reference ==="
echo ""

# PLINK2 --king-cutoff is the standard quick pruning command. It is useful as a
# reference check, but its greedy sample removal is phenotype-blind. For this
# rare-cancer tutorial, final removal decisions are made from the KING pair table
# below so that case-control pairs preferentially keep the case.
plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_06_filt" \
  --extract "${OUT_DIR}/${DATASET_NAME}_07_prune.prune.in" \
  --king-cutoff "${REMOVE_THRESHOLD}" \
  --out "${OUT_DIR}/${DATASET_NAME}_07_king_cutoff"

echo "✓ PLINK2 KING cutoff reference complete"

# ============================================================================
# STEP 3: Create PLINK1.9 IBD/PI_HAT report
# ============================================================================
echo ""
echo "=== Computing PLINK1.9 IBD/PI_HAT report ==="
echo ""

# PLINK1.9 --genome estimates identity-by-descent (IBD) sharing and writes
# PI_HAT = P(IBD=2) + 0.5 * P(IBD=1). This is useful for teaching and for
# comparison with older GWAS QC workflows. PI_HAT and KING kinship are on
# different scales, so final pruning below remains based on KING pairs.
plink \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_06_filt" \
  --extract "${OUT_DIR}/${DATASET_NAME}_07_prune.prune.in" \
  --genome \
  --min "${PIHAT_REPORT_THRESHOLD}" \
  --out "${OUT_DIR}/${DATASET_NAME}_07_ibd_pi_hat"

PIHAT_GENOME="${OUT_DIR}/${DATASET_NAME}_07_ibd_pi_hat.genome"
PIHAT_ANNOTATED="${OUT_DIR}/${DATASET_NAME}_07_ibd_pi_hat_annotated.tsv"
PIHAT_SUMMARY="${OUT_DIR}/${DATASET_NAME}_07_ibd_pi_hat_summary.tsv"

if [ ! -f "$PIHAT_GENOME" ]; then
  echo "✗ Expected PLINK1.9 PI_HAT report not found: $PIHAT_GENOME"
  exit 1
fi

awk '
  BEGIN {
    OFS = "\t"
    print "FID1", "IID1", "FID2", "IID2", "Z0", "Z1", "Z2", "PI_HAT", "RELATIONSHIP_CATEGORY"
  }
  function relationship_category(pi_hat) {
    pi_hat += 0
    if (pi_hat >= 0.70) return "Duplicate_or_twin"
    if (pi_hat >= 0.354) return "First_degree"
    if (pi_hat >= 0.177) return "Second_degree"
    if (pi_hat >= 0.0884) return "Third_degree"
    return "Unrelated_or_distant"
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "FID1") fid1 = i
      if ($i == "IID1") iid1 = i
      if ($i == "FID2") fid2 = i
      if ($i == "IID2") iid2 = i
      if ($i == "Z0") z0 = i
      if ($i == "Z1") z1 = i
      if ($i == "Z2") z2 = i
      if ($i == "PI_HAT") pi_hat = i
    }
    if (!fid1 || !iid1 || !fid2 || !iid2 || !z0 || !z1 || !z2 || !pi_hat) {
      print "Could not parse required .genome columns." > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    print $fid1, $iid1, $fid2, $iid2, $z0, $z1, $z2, $pi_hat, relationship_category($pi_hat)
  }
' "$PIHAT_GENOME" > "$PIHAT_ANNOTATED"

awk -F "\t" '
  BEGIN {
    OFS = "\t"
    cats[1] = "Duplicate_or_twin"
    cats[2] = "First_degree"
    cats[3] = "Second_degree"
    cats[4] = "Third_degree"
  }
  NR == 1 { next }
  {
    all_pair_count++
    pair_count[$9]++
    sample1 = $1 "/" $2
    sample2 = $3 "/" $4
    related_sample[sample1] = 1
    related_sample[sample2] = 1
    category_sample[$9, sample1] = 1
    category_sample[$9, sample2] = 1
  }
  END {
    for (sample in related_sample) {
      related_sample_count++
    }
    for (key in category_sample) {
      split(key, parts, SUBSEP)
      category_sample_count[parts[1]]++
    }

    print "relationship_category", "n_pairs", "n_unique_samples"
    for (i = 1; i <= 4; i++) {
      cat = cats[i]
      print cat, pair_count[cat] + 0, category_sample_count[cat] + 0
    }
    print "All_PI_HAT_pairs_reported", all_pair_count + 0, related_sample_count + 0
  }
' "$PIHAT_ANNOTATED" > "$PIHAT_SUMMARY"

echo "✓ PI_HAT report: ${PIHAT_GENOME}"
echo "✓ Annotated PI_HAT report: ${PIHAT_ANNOTATED}"
echo "✓ PI_HAT summary: ${PIHAT_SUMMARY}"
echo ""
cat "$PIHAT_SUMMARY"

# ============================================================================
# STEP 4: Create related-pair table and category summary
# ============================================================================
echo ""
echo "=== Summarizing relatedness categories ==="
echo ""

# For interpretation, we create a separate table of sample pairs with KING
# kinship >= REPORT_THRESHOLD. This includes approximate third-degree-or-closer
# pairs when REPORT_THRESHOLD is 0.0442.
#
# Common KING kinship interpretation:
#   >= 0.354  duplicate/twin
#   >= 0.177  first-degree
#   >= 0.0884 second-degree
#   >= 0.0442 third-degree

plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_06_filt" \
  --extract "${OUT_DIR}/${DATASET_NAME}_07_prune.prune.in" \
  --make-king-table \
  --king-table-filter "${REPORT_THRESHOLD}" \
  --out "${OUT_DIR}/${DATASET_NAME}_07_related_pairs"

RELATED_PAIRS_KIN0="${OUT_DIR}/${DATASET_NAME}_07_related_pairs.kin0"
RELATED_PAIRS_ANNOTATED="${OUT_DIR}/${DATASET_NAME}_07_related_pairs_annotated.tsv"
RELATEDNESS_SUMMARY="${OUT_DIR}/${DATASET_NAME}_07_relatedness_summary.tsv"
PRUNING_PAIRS="${OUT_DIR}/${DATASET_NAME}_07_related_pairs_for_pruning.tsv"
REMOVAL_DECISIONS="${OUT_DIR}/${DATASET_NAME}_07_removal_decisions.tsv"
REMOVAL_SUMMARY="${OUT_DIR}/${DATASET_NAME}_07_removal_summary.tsv"
KING_PIHAT_COMPARISON="${OUT_DIR}/${DATASET_NAME}_07_king_vs_pihat_comparison.tsv"
KING_PIHAT_SUMMARY="${OUT_DIR}/${DATASET_NAME}_07_king_vs_pihat_summary.tsv"
PRUNING_COMPONENTS="${OUT_DIR}/${DATASET_NAME}_07_pruning_components.tsv"
PRUNING_DIAGNOSTICS="${OUT_DIR}/${DATASET_NAME}_07_relatedness_pruning_diagnostics.tsv"

if [ ! -f "$RELATED_PAIRS_KIN0" ]; then
  echo "✗ Expected KING related-pair table not found: $RELATED_PAIRS_KIN0"
  exit 1
fi

awk -v remove_threshold="$REMOVE_THRESHOLD" '
  BEGIN {
    OFS = "\t"
    print "FID1", "IID1", "FID2", "IID2", "PHENO1", "PHENO2", "KINSHIP", "RELATIONSHIP_CATEGORY", "USED_FOR_PRUNING"
  }
  function sample_key(fid, iid) {
    return fid "/" iid
  }
  function phenotype_label(value) {
    value = tolower(value)
    if (value == "2" || value == "case" || value == "affected") return "case"
    if (value == "1" || value == "control" || value == "unaffected") return "control"
    return "unknown"
  }
  function relationship_category(k) {
    k += 0
    if (k >= 0.354) return "Duplicate_or_twin"
    if (k >= 0.177) return "First_degree"
    if (k >= 0.0884) return "Second_degree"
    if (k >= 0.0442) return "Third_degree"
    return "Unrelated_or_distant"
  }
  FNR == NR {
    if (FNR == 1) {
      first = $1
      second = $2
      third = $3
      gsub(/^#/, "", first)
      gsub(/^#/, "", second)
      if (first == "FID" && (second == "IID" || second == "ID")) {
        next
      }
    }
    if (NF >= 3) {
      pheno[sample_key($1, $2)] = phenotype_label($3)
    }
    next
  }
  FNR == 1 {
    for (i = 1; i <= NF; i++) {
      col = $i
      gsub(/^#/, "", col)
      if (col == "FID1") fid1 = i
      if (col == "ID1" || col == "IID1") iid1 = i
      if (col == "FID2") fid2 = i
      if (col == "ID2" || col == "IID2") iid2 = i
      if (col == "KINSHIP") kinship = i
    }
    if (!iid1 || !iid2 || !kinship) {
      print "Could not parse required .kin0 columns." > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    sample_fid1 = fid1 ? $fid1 : 0
    sample_fid2 = fid2 ? $fid2 : 0
    sample1 = sample_key(sample_fid1, $iid1)
    sample2 = sample_key(sample_fid2, $iid2)
    pheno1 = sample1 in pheno ? pheno[sample1] : "unknown"
    pheno2 = sample2 in pheno ? pheno[sample2] : "unknown"
    k = $kinship + 0
    used_for_pruning = k > remove_threshold ? "yes" : "no"
    print sample_fid1, $iid1, sample_fid2, $iid2, pheno1, pheno2, k, relationship_category(k), used_for_pruning
  }
' "$PHENOTYPE_FILE" "$RELATED_PAIRS_KIN0" > "$RELATED_PAIRS_ANNOTATED"

awk '
  BEGIN {
    OFS = "\t"
    cats[1] = "Duplicate_or_twin"
    cats[2] = "First_degree"
    cats[3] = "Second_degree"
    cats[4] = "Third_degree"
  }
  function relationship_category(k) {
    k += 0
    if (k >= 0.354) return "Duplicate_or_twin"
    if (k >= 0.177) return "First_degree"
    if (k >= 0.0884) return "Second_degree"
    if (k >= 0.0442) return "Third_degree"
    return "Unrelated_or_distant"
  }
  FNR == NR {
    total_samples++
    next
  }
  FNR == 1 {
    for (i = 1; i <= NF; i++) {
      col = $i
      gsub(/^#/, "", col)
      if (col == "FID1") fid1 = i
      if (col == "ID1" || col == "IID1") iid1 = i
      if (col == "FID2") fid2 = i
      if (col == "ID2" || col == "IID2") iid2 = i
      if (col == "KINSHIP") kinship = i
    }
    if (!iid1 || !iid2 || !kinship) {
      print "Could not parse required .kin0 columns." > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    sample_fid1 = fid1 ? $fid1 : 0
    sample_fid2 = fid2 ? $fid2 : 0
    sample1 = sample_fid1 "/" $iid1
    sample2 = sample_fid2 "/" $iid2
    cat = relationship_category($kinship)
    pair_count[cat]++
    related_sample[sample1] = 1
    related_sample[sample2] = 1
    category_sample[cat, sample1] = 1
    category_sample[cat, sample2] = 1
  }
  END {
    for (sample in related_sample) {
      related_sample_count++
    }
    for (key in category_sample) {
      split(key, parts, SUBSEP)
      category_sample_count[parts[1]]++
    }

    print "relationship_category", "n_pairs", "n_unique_samples"
    for (i = 1; i <= 4; i++) {
      cat = cats[i]
      print cat, pair_count[cat] + 0, category_sample_count[cat] + 0
    }
    print "No_third_degree_or_closer_pair", "NA", total_samples - related_sample_count
    print "Total_samples_before_relatedness_filter", "NA", total_samples
  }
' "${DATASET_INPUT}/${DATASET_NAME}_06_filt.fam" "$RELATED_PAIRS_KIN0" > "$RELATEDNESS_SUMMARY"

echo "✓ Related-pair table: ${RELATED_PAIRS_KIN0}"
echo "✓ Annotated related-pair table: ${RELATED_PAIRS_ANNOTATED}"
echo "✓ Relatedness summary: ${RELATEDNESS_SUMMARY}"
echo ""
cat "$RELATEDNESS_SUMMARY"

# ============================================================================
# STEP 5: Compare KING kinship and PLINK1.9 PI_HAT
# ============================================================================
echo ""
echo "=== Comparing KING kinship and PLINK1.9 PI_HAT ==="
echo ""

# KING kinship and PI_HAT are related but not interchangeable. In ideal outbred
# pedigrees, PI_HAT is approximately 2 * KING kinship. This comparison table is
# for QC review and teaching; the final removal list below still uses KING.
awk -F "\t" '
  BEGIN {
    OFS = "\t"
    print "FID1", "IID1", "FID2", "IID2", "KING_KINSHIP", "KING_CATEGORY", "PI_HAT", "PI_HAT_CATEGORY", "EXPECTED_PI_HAT_FROM_KING", "ABS_PIHAT_MINUS_2KING", "USED_FOR_PRUNING", "COMPARISON_NOTE"
  }
  function pair_key(fid1, iid1, fid2, iid2,    sample1, sample2) {
    sample1 = fid1 "/" iid1
    sample2 = fid2 "/" iid2
    return sample1 <= sample2 ? sample1 SUBSEP sample2 : sample2 SUBSEP sample1
  }
  function store_display(key, fid1, iid1, fid2, iid2,    sample1, sample2) {
    if (key in display_fid1) return
    sample1 = fid1 "/" iid1
    sample2 = fid2 "/" iid2
    if (sample1 <= sample2) {
      display_fid1[key] = fid1
      display_iid1[key] = iid1
      display_fid2[key] = fid2
      display_iid2[key] = iid2
    } else {
      display_fid1[key] = fid2
      display_iid1[key] = iid2
      display_fid2[key] = fid1
      display_iid2[key] = iid1
    }
  }
  FNR == NR {
    if (FNR == 1) next
    key = pair_key($1, $2, $3, $4)
    store_display(key, $1, $2, $3, $4)
    king_kinship[key] = $7
    king_category[key] = $8
    used_for_pruning[key] = $9
    seen[key] = 1
    next
  }
  {
    if (FNR == 1) next
    key = pair_key($1, $2, $3, $4)
    store_display(key, $1, $2, $3, $4)
    pihat[key] = $8
    pihat_category[key] = $9
    seen[key] = 1
  }
  END {
    for (key in seen) {
      k = key in king_kinship ? king_kinship[key] : "NA"
      kcat = key in king_category ? king_category[key] : "NA"
      p = key in pihat ? pihat[key] : "NA"
      pcat = key in pihat_category ? pihat_category[key] : "NA"
      used = key in used_for_pruning ? used_for_pruning[key] : "no"

      if (k != "NA") {
        expected = sprintf("%.6f", 2 * k)
      } else {
        expected = "NA"
      }

      if (k != "NA" && p != "NA") {
        abs_diff = p - (2 * k)
        if (abs_diff < 0) abs_diff *= -1
        abs_diff = sprintf("%.6f", abs_diff)
        note = kcat == pcat ? "same_category" : "different_category_check_values"
      } else if (k != "NA") {
        abs_diff = "NA"
        note = "reported_by_king_only"
      } else {
        abs_diff = "NA"
        note = "reported_by_pihat_only"
      }

      print display_fid1[key], display_iid1[key], display_fid2[key], display_iid2[key], k, kcat, p, pcat, expected, abs_diff, used, note
    }
  }
' "$RELATED_PAIRS_ANNOTATED" "$PIHAT_ANNOTATED" > "$KING_PIHAT_COMPARISON"

awk -F "\t" '
  BEGIN {
    OFS = "\t"
  }
  NR == 1 { next }
  {
    total_pairs++
    note_count[$12]++
    if ($5 != "NA" && $7 != "NA") {
      compared_pairs++
      sum_abs_diff += $10
      if ($10 > max_abs_diff) max_abs_diff = $10
    }
    if ($11 == "yes") {
      pruning_pairs++
      if ($7 != "NA") pruning_pairs_with_pihat++
      else pruning_pairs_without_pihat++
    }
  }
  END {
    print "metric", "value"
    print "total_pairs_in_either_report", total_pairs + 0
    print "pairs_reported_by_both", compared_pairs + 0
    print "same_category", note_count["same_category"] + 0
    print "different_category_check_values", note_count["different_category_check_values"] + 0
    print "reported_by_king_only", note_count["reported_by_king_only"] + 0
    print "reported_by_pihat_only", note_count["reported_by_pihat_only"] + 0
    print "pruning_pairs", pruning_pairs + 0
    print "pruning_pairs_with_pihat", pruning_pairs_with_pihat + 0
    print "pruning_pairs_without_pihat", pruning_pairs_without_pihat + 0
    if (compared_pairs > 0) {
      print "mean_abs_pihat_minus_2king", sprintf("%.6f", sum_abs_diff / compared_pairs)
      print "max_abs_pihat_minus_2king", sprintf("%.6f", max_abs_diff)
    } else {
      print "mean_abs_pihat_minus_2king", "NA"
      print "max_abs_pihat_minus_2king", "NA"
    }
  }
' "$KING_PIHAT_COMPARISON" > "$KING_PIHAT_SUMMARY"

echo "✓ KING vs PI_HAT comparison: ${KING_PIHAT_COMPARISON}"
echo "✓ KING vs PI_HAT summary: ${KING_PIHAT_SUMMARY}"
echo ""
cat "$KING_PIHAT_SUMMARY"

# ============================================================================
# STEP 6: Create phenotype-aware removal list
# ============================================================================
echo ""
echo "=== Creating phenotype-aware relatedness removal list ==="
echo ""

# Rare cancer rule:
#   - If a related pair is case-control, keep the case and remove the control.
#   - If both samples have the same phenotype, remove the sample with more
#     related pairs; if tied, use a deterministic ID tie-breaker.
#   - A case is removed only when paired with another case or when unavoidable.

awk -F "\t" 'NR > 1 && $9 == "yes" { print }' "$RELATED_PAIRS_ANNOTATED" |
  sort -t "$(printf '\t')" -k7,7gr > "$PRUNING_PAIRS"

awk -F "\t" '
  BEGIN {
    OFS = "\t"
  }
  function sample_key(fid, iid) {
    return fid "/" iid
  }
  function normalize_pheno(pheno) {
    pheno = tolower(pheno)
    gsub(/\r/, "", pheno)
    if (pheno == "2" || pheno == "case" || pheno == "affected") return "case"
    if (pheno == "1" || pheno == "control" || pheno == "unaffected") return "control"
    return "unknown"
  }
  function add_sample(sample, pheno) {
    if (!(sample in parent)) {
      parent[sample] = sample
      sample_pheno[sample] = normalize_pheno(pheno)
    }
  }
  function find(sample,    root, current, next_sample) {
    root = sample
    while (parent[root] != root) {
      root = parent[root]
    }
    current = sample
    while (parent[current] != current) {
      next_sample = parent[current]
      parent[current] = root
      current = next_sample
    }
    return root
  }
  function union_samples(sample1, sample2,    root1, root2) {
    root1 = find(sample1)
    root2 = find(sample2)
    if (root1 != root2) {
      parent[root2] = root1
    }
  }
  {
    sample1 = sample_key($1, $2)
    sample2 = sample_key($3, $4)
    add_sample(sample1, $5)
    add_sample(sample2, $6)
    union_samples(sample1, sample2)
    pair_sample1[++n_pairs] = sample1
    pair_kinship[n_pairs] = $7 + 0
  }
  END {
    print "component_id", "n_samples", "n_pairs", "n_controls", "n_cases", "n_unknown", "max_kinship", "mean_kinship"

    for (sample in parent) {
      root = find(sample)
      if (!(root in component_id)) {
        component_id[root] = ++n_components
      }
      comp = component_id[root]
      component_samples[comp]++
      pheno = sample_pheno[sample]
      if (pheno == "control") component_controls[comp]++
      else if (pheno == "case") component_cases[comp]++
      else component_unknown[comp]++
    }

    for (i = 1; i <= n_pairs; i++) {
      comp = component_id[find(pair_sample1[i])]
      component_pairs[comp]++
      component_sum_kinship[comp] += pair_kinship[i]
      if (pair_kinship[i] > component_max_kinship[comp]) {
        component_max_kinship[comp] = pair_kinship[i]
      }
    }

    for (comp = 1; comp <= n_components; comp++) {
      if (component_pairs[comp] > 0) {
        mean_kinship = component_sum_kinship[comp] / component_pairs[comp]
      } else {
        mean_kinship = 0
      }
      print comp, component_samples[comp] + 0, component_pairs[comp] + 0, component_controls[comp] + 0, component_cases[comp] + 0, component_unknown[comp] + 0, sprintf("%.6f", component_max_kinship[comp] + 0), sprintf("%.6f", mean_kinship)
    }
  }
' "$PRUNING_PAIRS" > "${PRUNING_COMPONENTS}.unsorted"

{
  head -n 1 "${PRUNING_COMPONENTS}.unsorted"
  tail -n +2 "${PRUNING_COMPONENTS}.unsorted" | sort -t "$(printf '\t')" -k2,2nr -k3,3nr
} > "$PRUNING_COMPONENTS"
rm -f "${PRUNING_COMPONENTS}.unsorted"

UNKNOWN_PRUNING_PHENOTYPES=$(awk -F "\t" '
  ($5 == "unknown" || $5 == "" || $6 == "unknown" || $6 == "") { count++ }
  END { print count + 0 }
' "$PRUNING_PAIRS")

if [ "$UNKNOWN_PRUNING_PHENOTYPES" -gt 0 ]; then
  echo "✗ Some related pairs used for pruning do not have case/control labels."
  echo "  This step is phenotype-aware, so every sample in a pruned pair must be labelled."
  echo "  Check that ${PHENOTYPE_FILE} contains FID IID PHENO for all samples in ${DATASET_INPUT}/${DATASET_NAME}_06_filt.fam."
  echo ""
  echo "Examples:"
  awk -F "\t" '($5 == "unknown" || $5 == "" || $6 == "unknown" || $6 == "") { print; shown++; if (shown == 5) exit }' "$PRUNING_PAIRS"
  exit 1
fi

: > "${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt"

awk -F "\t" \
  -v decisions="$REMOVAL_DECISIONS" \
  -v removal_list="${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt" \
  -v summary="$REMOVAL_SUMMARY" '
  BEGIN {
    OFS = "\t"
    print "REMOVED_FID", "REMOVED_IID", "REMOVED_PHENO", "KEPT_FID", "KEPT_IID", "KEPT_PHENO", "KINSHIP", "RELATIONSHIP_CATEGORY", "REASON" > decisions
  }
  function sample_key(fid, iid) {
    return fid "/" iid
  }
  function normalize_pheno(pheno) {
    pheno = tolower(pheno)
    gsub(/\r/, "", pheno)
    if (pheno == "2" || pheno == "case" || pheno == "affected") return "case"
    if (pheno == "1" || pheno == "control" || pheno == "unaffected") return "control"
    return "unknown"
  }
  function pheno_rank(pheno) {
    pheno = normalize_pheno(pheno)
    if (pheno == "case") return 3
    if (pheno == "control") return 2
    return 1
  }
  function choose_removal(sample1, pheno1, sample2, pheno2) {
    rank1 = pheno_rank(pheno1)
    rank2 = pheno_rank(pheno2)
    if (rank1 != rank2) {
      return rank1 < rank2 ? sample1 : sample2
    }
    if (degree[sample1] != degree[sample2]) {
      return degree[sample1] >= degree[sample2] ? sample1 : sample2
    }
    return sample1 > sample2 ? sample1 : sample2
  }
  function decision_reason(removed_pheno, kept_pheno) {
    if (kept_pheno == "case" && removed_pheno == "control") {
      return "case_control_pair_keep_case_remove_control"
    }
    if (kept_pheno == "case" && removed_pheno == "unknown") {
      return "case_unknown_pair_keep_case_remove_unknown"
    }
    if (kept_pheno == "control" && removed_pheno == "unknown") {
      return "control_unknown_pair_keep_control_remove_unknown"
    }
    if (removed_pheno == kept_pheno) {
      return "same_phenotype_remove_higher_related_degree_or_id"
    }
    return "phenotype_priority_remove_lower_priority_sample"
  }
  FNR == NR {
    sample1 = sample_key($1, $2)
    sample2 = sample_key($3, $4)
    pheno[sample1] = normalize_pheno($5)
    pheno[sample2] = normalize_pheno($6)
    degree[sample1]++
    degree[sample2]++
    next
  }
  {
    pruning_pairs_evaluated++
    sample1 = sample_key($1, $2)
    sample2 = sample_key($3, $4)
    pheno1 = (sample1 in pheno && pheno[sample1] != "") ? pheno[sample1] : normalize_pheno($5)
    pheno2 = (sample2 in pheno && pheno[sample2] != "") ? pheno[sample2] : normalize_pheno($6)
    if ((sample1 in removed) || (sample2 in removed)) {
      pruning_pairs_already_resolved++
      next
    }

    sample_to_remove = choose_removal(sample1, pheno1, sample2, pheno2)
    if (sample_to_remove == sample1) {
      sample_to_keep = sample2
      removed_fid = $1
      removed_iid = $2
      kept_fid = $3
      kept_iid = $4
    } else {
      sample_to_keep = sample1
      removed_fid = $3
      removed_iid = $4
      kept_fid = $1
      kept_iid = $2
    }

    removed[sample_to_remove] = 1
    removed_pheno = (sample_to_remove in pheno && pheno[sample_to_remove] != "") ? pheno[sample_to_remove] : "unknown"
    kept_pheno = (sample_to_keep in pheno && pheno[sample_to_keep] != "") ? pheno[sample_to_keep] : "unknown"
    removed_pheno = normalize_pheno(removed_pheno)
    kept_pheno = normalize_pheno(kept_pheno)
    reason = decision_reason(removed_pheno, kept_pheno)

    removed_count[removed_pheno]++
    reason_count[reason]++
    print removed_fid, removed_iid, removed_pheno, kept_fid, kept_iid, kept_pheno, $7, $8, reason > decisions
  }
  END {
    for (sample in removed) {
      split(sample, parts, "/")
      print parts[1], parts[2] > removal_list
      total_removed++
    }

    print "metric", "value" > summary
    print "pruning_pairs_evaluated", pruning_pairs_evaluated + 0 > summary
    print "pruning_pairs_already_resolved_by_prior_removal", pruning_pairs_already_resolved + 0 > summary
    print "controls_removed", removed_count["control"] + 0 > summary
    print "cases_removed", removed_count["case"] + 0 > summary
    print "unknown_phenotype_removed", removed_count["unknown"] + 0 > summary
    print "total_removed", total_removed + 0 > summary
    for (reason in reason_count) {
      print reason, reason_count[reason] > summary
    }
  }
' "$PRUNING_PAIRS" "$PRUNING_PAIRS"

sort -u -k1,1 -k2,2 "${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt" \
  > "${OUT_DIR}/${DATASET_NAME}_07_removal_list.sorted.txt"
mv "${OUT_DIR}/${DATASET_NAME}_07_removal_list.sorted.txt" \
  "${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt"

NREMOVE=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt")
SUMMARY_TOTAL_REMOVED=$(awk '$1 == "total_removed" {print $2}' "$REMOVAL_SUMMARY")
NPRUNE_PAIRS_FOR_CHECK=$(wc -l < "$PRUNING_PAIRS")
NREMOVAL_DECISIONS=$(awk 'NR > 1 { count++ } END { print count + 0 }' "$REMOVAL_DECISIONS")

if [ "$NREMOVE" -gt "$NPRUNE_PAIRS_FOR_CHECK" ]; then
  echo "✗ Removal list has more samples than pruning pairs."
  echo "  Pruning pairs: ${NPRUNE_PAIRS_FOR_CHECK}"
  echo "  Removal list lines: ${NREMOVE}"
  echo "  This should not happen when removing at most one sample per pair."
  echo "  Please inspect ${REMOVAL_DECISIONS}."
  exit 1
fi

if [ "$NREMOVE" -ne "$NREMOVAL_DECISIONS" ]; then
  echo "✗ Removal decision table does not match removal list count."
  echo "  Removal decisions: ${NREMOVAL_DECISIONS}"
  echo "  Removal list lines: ${NREMOVE}"
  echo "  Please inspect ${REMOVAL_DECISIONS}."
  exit 1
fi

if [ "$NREMOVE" -ne "$SUMMARY_TOTAL_REMOVED" ]; then
  echo "✗ Removal list count does not match removal summary."
  echo "  Removal list lines: ${NREMOVE}"
  echo "  Summary total_removed: ${SUMMARY_TOTAL_REMOVED}"
  echo "  Please inspect ${REMOVAL_DECISIONS}."
  exit 1
fi

REMOVAL_PHENO_COUNTS="${OUT_DIR}/${DATASET_NAME}_07_removal_phenotype_counts.tsv"

awk '
  BEGIN {
    OFS = "\t"
  }
  function sample_key(fid, iid) {
    return fid "/" iid
  }
  function phenotype_label(value) {
    value = tolower(value)
    gsub(/\r/, "", value)
    if (value == "2" || value == "case" || value == "affected") return "case"
    if (value == "1" || value == "control" || value == "unaffected") return "control"
    return "unknown"
  }
  FNR == NR {
    first = $1
    second = $2
    gsub(/^#/, "", first)
    gsub(/^#/, "", second)
    if (FNR == 1 && first == "FID" && (second == "IID" || second == "ID")) {
      next
    }
    if (NF >= 3) {
      pheno[sample_key($1, $2)] = phenotype_label($3)
    }
    next
  }
  {
    key = sample_key($1, $2)
    label = key in pheno ? pheno[key] : "unknown"
    removed_count[label]++
    total_removed++
  }
  END {
    print "phenotype", "n_removed"
    print "control", removed_count["control"] + 0
    print "case", removed_count["case"] + 0
    print "unknown", removed_count["unknown"] + 0
    print "total", total_removed + 0
  }
' "$PHENOTYPE_FILE" "${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt" > "$REMOVAL_PHENO_COUNTS"

PHENO_TOTAL_REMOVED=$(awk '$1 == "total" {print $2}' "$REMOVAL_PHENO_COUNTS")
if [ "$NREMOVE" -ne "$PHENO_TOTAL_REMOVED" ]; then
  echo "✗ Removal phenotype count does not match removal list count."
  echo "  Removal list lines: ${NREMOVE}"
  echo "  Phenotype count total: ${PHENO_TOTAL_REMOVED}"
  echo "  Please inspect ${REMOVAL_PHENO_COUNTS}."
  exit 1
fi

KING_CUTOFF_REMOVE="${OUT_DIR}/${DATASET_NAME}_07_king_cutoff.king.cutoff.out.id"
if [ -f "$KING_CUTOFF_REMOVE" ]; then
  NPLINK2_KING_REMOVED=$(awk '
    NR == 1 {
      first = $1
      second = $2
      gsub(/^#/, "", first)
      gsub(/^#/, "", second)
      if (first == "FID" || first == "IID" || second == "IID") {
        next
      }
    }
    NF >= 1 { count++ }
    END { print count + 0 }
  ' "$KING_CUTOFF_REMOVE")
else
  NPLINK2_KING_REMOVED=0
fi

NSAMP_BEFORE_DIAG=$(wc -l < "${DATASET_INPUT}/${DATASET_NAME}_06_filt.fam")
NPRUNE_PAIRS=$(wc -l < "$PRUNING_PAIRS")
NPRUNE_SAMPLES=$(awk -F "\t" '
  {
    samples[$1 "/" $2] = 1
    samples[$3 "/" $4] = 1
  }
  END {
    for (sample in samples) count++
    print count + 0
  }
' "$PRUNING_PAIRS")
NPRUNE_COMPONENTS=$(awk 'NR > 1 { count++ } END { print count + 0 }' "$PRUNING_COMPONENTS")
LARGEST_COMPONENT_SAMPLES=$(awk -F "\t" 'NR > 1 && $2 > max { max = $2 } END { print max + 0 }' "$PRUNING_COMPONENTS")
LARGEST_COMPONENT_PAIRS=$(awk -F "\t" 'NR > 1 && $2 > max_samples { max_samples = $2; max_pairs = $3 } END { print max_pairs + 0 }' "$PRUNING_COMPONENTS")
RELATEDNESS_LOSS_PCT=$(awk -v removed="$NREMOVE" -v total="$NSAMP_BEFORE_DIAG" 'BEGIN { if (total > 0) printf "%.1f", 100 * removed / total; else print "NA" }')
NCONTROL_REMOVED_DIAG=$(awk '$1 == "control" {print $2}' "$REMOVAL_PHENO_COUNTS")
NCASE_REMOVED_DIAG=$(awk '$1 == "case" {print $2}' "$REMOVAL_PHENO_COUNTS")
NUNKNOWN_REMOVED_DIAG=$(awk '$1 == "unknown" {print $2}' "$REMOVAL_PHENO_COUNTS")

{
  echo -e "metric\tvalue"
  echo -e "samples_before_relatedness_filter\t${NSAMP_BEFORE_DIAG}"
  echo -e "king_pruning_pairs_gt_${REMOVE_THRESHOLD}\t${NPRUNE_PAIRS}"
  echo -e "unique_samples_in_king_pruning_pairs\t${NPRUNE_SAMPLES}"
  echo -e "relatedness_components\t${NPRUNE_COMPONENTS}"
  echo -e "largest_component_samples\t${LARGEST_COMPONENT_SAMPLES}"
  echo -e "largest_component_pairs\t${LARGEST_COMPONENT_PAIRS}"
  echo -e "plink2_king_cutoff_removed_samples\t${NPLINK2_KING_REMOVED}"
  echo -e "phenotype_aware_removed_samples\t${NREMOVE}"
  echo -e "phenotype_aware_removed_percent\t${RELATEDNESS_LOSS_PCT}"
  echo -e "controls_removed\t${NCONTROL_REMOVED_DIAG}"
  echo -e "cases_removed\t${NCASE_REMOVED_DIAG}"
  echo -e "unknown_phenotype_removed\t${NUNKNOWN_REMOVED_DIAG}"
} > "$PRUNING_DIAGNOSTICS"

echo "✓ Phenotype-aware removal decisions: ${REMOVAL_DECISIONS}"
echo "✓ Removal summary: ${REMOVAL_SUMMARY}"
echo "✓ Removal phenotype counts: ${REMOVAL_PHENO_COUNTS}"
echo "✓ Relatedness pruning components: ${PRUNING_COMPONENTS}"
echo "✓ Relatedness pruning diagnostics: ${PRUNING_DIAGNOSTICS}"
echo "Related samples to remove: ${NREMOVE}"
echo ""
cat "$REMOVAL_SUMMARY"
echo ""
cat "$REMOVAL_PHENO_COUNTS"
echo ""
cat "$PRUNING_DIAGNOSTICS"

if awk -v pct="$RELATEDNESS_LOSS_PCT" 'BEGIN { exit !(pct + 0 >= 15) }'; then
  echo ""
  echo "⚠ Relatedness pruning removed ${RELATEDNESS_LOSS_PCT}% of samples from the Step 06 dataset."
  echo "  This is a large sample loss. Inspect:"
  echo "  - ${PRUNING_COMPONENTS}"
  echo "  - ${PRUNING_DIAGNOSTICS}"
  echo "  - ${REMOVAL_DECISIONS}"
  echo "  If this is real family structure, consider documenting it or using a relatedness-aware association method downstream."
fi

# ============================================================================
# STEP 7: Create kinship-pruned dataset
# ============================================================================
echo ""
echo "=== Creating pruned (unrelated) dataset ==="
echo ""

if [ -s "${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt" ]; then
  plink2 \
    --bfile "${DATASET_INPUT}/${DATASET_NAME}_06_filt" \
    --remove "${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt" \
    --make-bed \
    --out "${OUT_DIR}/${DATASET_NAME}_07_filt"
else
  plink2 \
    --bfile "${DATASET_INPUT}/${DATASET_NAME}_06_filt" \
    --make-bed \
    --out "${OUT_DIR}/${DATASET_NAME}_07_filt"
fi

echo "✓ Relatedness-pruned dataset: ${OUT_DIR}/${DATASET_NAME}_07_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NSAMP_BEFORE=$(wc -l < "${DATASET_INPUT}/${DATASET_NAME}_06_filt.fam")
NSAMP_AFTER=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_07_filt.fam")
NSAMP_REMOVED=$((NSAMP_BEFORE - NSAMP_AFTER))
NCONTROL_REMOVED=$(awk '$1 == "control" {print $2}' "$REMOVAL_PHENO_COUNTS")
NCASE_REMOVED=$(awk '$1 == "case" {print $2}' "$REMOVAL_PHENO_COUNTS")
NUNKNOWN_REMOVED=$(awk '$1 == "unknown" {print $2}' "$REMOVAL_PHENO_COUNTS")
NPHENO_TOTAL_REMOVED=$(awk '$1 == "total" {print $2}' "$REMOVAL_PHENO_COUNTS")
NPHENO_SUM=$((NCONTROL_REMOVED + NCASE_REMOVED + NUNKNOWN_REMOVED))

if [ "$NSAMP_REMOVED" -ne "$NREMOVE" ]; then
  echo "✗ PLINK removed a different number of samples than the removal list contains."
  echo "  Removal list count: ${NREMOVE}"
  echo "  Dataset count difference: ${NSAMP_REMOVED}"
  echo "  Check ${OUT_DIR}/${DATASET_NAME}_07_removal_list.txt and ${OUT_DIR}/${DATASET_NAME}_07_filt.log."
  exit 1
fi

if [ "$NSAMP_REMOVED" -ne "$NPHENO_TOTAL_REMOVED" ]; then
  echo "✗ Case/control removal counts do not add up to the removed sample count."
  echo "  Dataset count difference: ${NSAMP_REMOVED}"
  echo "  Phenotype count total: ${NPHENO_TOTAL_REMOVED}"
  echo "  Check ${REMOVAL_PHENO_COUNTS}."
  exit 1
fi

if [ "$NSAMP_REMOVED" -ne "$NPHENO_SUM" ]; then
  echo "✗ Printed case/control/unknown counts do not add up to the removed sample count."
  echo "  Dataset count difference: ${NSAMP_REMOVED}"
  echo "  control + case + unknown: ${NPHENO_SUM}"
  echo "  Check ${REMOVAL_PHENO_COUNTS}."
  exit 1
fi

echo "Samples before relatedness filter: ${NSAMP_BEFORE}"
echo "Samples after relatedness filter:  ${NSAMP_AFTER}"
echo "Related samples removed:          ${NSAMP_REMOVED}"
echo "  - Controls removed:             ${NCONTROL_REMOVED}"
echo "  - Cases removed:                ${NCASE_REMOVED}"
echo "  - Unknown phenotype removed:    ${NUNKNOWN_REMOVED}"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the MAF (minor allele frequency) filter:"
echo ""
echo "  bash scripts/01B_genotyping_qc/08_maf_filter.sh"
echo ""
echo "This will:"
echo "  - Apply MAF threshold (--maf, context-dependent)"
echo "  - For PDAC: use --maf 0.01 to retain rare variants"
echo "  - Create final clean dataset for association testing"
echo ""
