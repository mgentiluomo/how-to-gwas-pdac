#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 09: QC Summary & Decision Tree
# 
# PURPOSE:
#   Create a comprehensive summary of the QC pipeline:
#   - Track sample/variant counts through each filter
#   - Compute filtering efficiency and power implications
#   - Generate a decision tree visualization for documentation
#
# INPUT:
#   - All intermediate datasets from Steps 01-08
#
# OUTPUT:
#   - pdac_demo_qc_summary.txt — detailed filtering report
#   - pdac_demo_qc_decision_tree.pdf — flowchart showing QC path
#   - pdac_demo_qc_retention_rates.pdf — visualization of filtering impact
#
# NOTES:
#   This script uses R to create reproducible visualizations
#   All filtering thresholds are documented with rationale
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

mkdir -p "$OUT_DIR"

# ============================================================================
# STEP 1: Collect sample and variant counts
# ============================================================================
echo ""
echo "=== Collecting QC metrics ==="
echo ""

# Function to extract counts
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

# Raw data
read NSAMP_RAW NVAR_RAW <<< "$(get_counts "${DATASET_INPUT}/${DATASET_NAME}")"
echo "Raw: ${NSAMP_RAW} samples, ${NVAR_RAW} variants"

# After each step
read NSAMP_02 NVAR_02 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_02_filt")"
read NSAMP_03 NVAR_03 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_03_filt")"
read NSAMP_04 NVAR_04 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_04_filt")"
read NSAMP_05 NVAR_05 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_05_filt")"
read NSAMP_06 NVAR_06 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_06_filt")"
read NSAMP_07 NVAR_07 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_07_filt")"
read NSAMP_08 NVAR_08 <<< "$(get_counts "${OUT_DIR}/${DATASET_NAME}_08_filt")"

# ============================================================================
# STEP 2: Create detailed summary report
# ============================================================================
echo ""
echo "=== Creating summary report ==="
echo ""

cat > "${OUT_DIR}/${DATASET_NAME}_qc_summary.txt" << EOF
================================================================================
PDAC GWAS Genotyping Quality Control Summary
================================================================================

Dataset: ${DATASET_NAME}
Analysis Date: $(date)
Genome Build: GRCh38

================================================================================
QC PIPELINE SUMMARY
================================================================================

Step                          Samples    Variants   Samples Lost  Variants Lost
─────────────────────────────────────────────────────────────────────────────
00. Raw data                  ${NSAMP_RAW}      ${NVAR_RAW}
01. Initial stats             ${NSAMP_RAW}      ${NVAR_RAW}        —             —
02. Sample call rate          ${NSAMP_02}      ${NVAR_02}        $((NSAMP_RAW - NSAMP_02))        0
03. Sex check                 ${NSAMP_03}      ${NVAR_03}        $((NSAMP_02 - NSAMP_03))        0
04. Heterozygosity            ${NSAMP_04}      ${NVAR_04}        $((NSAMP_03 - NSAMP_04))        0
05. Variant call rate         ${NSAMP_05}      ${NVAR_05}        0        $((NVAR_04 - NVAR_05))
06. Hardy-Weinberg            ${NSAMP_06}      ${NVAR_06}        0        $((NVAR_05 - NVAR_06))
07. Relatedness (IBD)         ${NSAMP_07}      ${NVAR_07}        $((NSAMP_06 - NSAMP_07))        0
08. MAF filter (≥1%)          ${NSAMP_08}      ${NVAR_08}        0        $((NVAR_07 - NVAR_08))
─────────────────────────────────────────────────────────────────────────────
FINAL DATASET                 ${NSAMP_08}      ${NVAR_08}        $((NSAMP_RAW - NSAMP_08))        $((NVAR_RAW - NVAR_08))

================================================================================
RETENTION RATES
================================================================================

Sample retention:  $((NSAMP_08 * 100 / NSAMP_RAW))% of raw samples
Variant retention: $((NVAR_08 * 100 / NVAR_RAW))% of raw variants

PDAC-specific notes:
  - Sample-level filters: call rate, sex check, heterozygosity, relatedness
  - Variant-level filters: call rate, HWE (controls), MAF ≥ 1%
  - Low MAF threshold (1% vs standard 5%) preserves rare-variant signal
  - Final dataset: ready for population stratification analysis (PCA)

================================================================================
QC THRESHOLDS APPLIED
================================================================================

Sample-level:
  --mind 0.02         (sample call rate: keep samples with <2% missing)
  Sex check           (F-statistic on chrX: remove discordant samples)
  Heterozygosity      (F ± 3 SD on autosomes: remove outliers)
  --king-cutoff 0.1875 (kinship: remove related pairs, PI_HAT > 0.1875)

Variant-level:
  --geno 0.05         (variant call rate: keep variants with <5% missing)
  --hwe 1e-6          (Hardy-Weinberg p-value in controls: p > 1e-6)
  --maf 0.01          (minor allele frequency: keep variants with MAF ≥ 1%)

================================================================================
QUALITY CONTROL DECISION RATIONALE
================================================================================

1. CALL RATE (sample level, --mind 0.02):
   Removes samples with >2% missing genotypes, indicating platform or lab issues.
   Threshold: 2% is standard; PDAC keeps this for sample retention.

2. SEX CHECK:
   Identifies and removes sex-discordant samples (genetic sex ≠ reported sex).
   Catches sample swaps and data entry errors using X chromosome heterozygosity.

3. HETEROZYGOSITY:
   Removes samples with extreme inbreeding coefficients (F ± 3 SD).
   Identifies contaminated samples (high F) or population outliers (low F).

4. VARIANT CALL RATE (--geno 0.05):
   Removes variants with >5% missing genotypes across the sample.
   Indicates assay design problems or systematic genotyping failures.

5. HARDY-WEINBERG EQUILIBRIUM (--hwe 1e-6):
   Tests HWE in controls only; disease variants may deviate in cases.
   Removes variants with p < 1e-6, indicating potential genotyping errors.
   Threshold 1e-6 is conservative; provides Bonferroni-like protection.

6. RELATEDNESS (--king-cutoff 0.1875):
   Identifies related pairs (2nd-degree relatives: siblings, parent-child).
   Removes one sample per related pair to maintain statistical independence.
   Important for disease cohorts (family recruitment common in cancer GWAS).

7. MINOR ALLELE FREQUENCY (--maf 0.01):
   **PDAC-specific decision: use 1% threshold (not standard 5%)**
   Rationale:
     - Rare variants carry larger effect sizes in complex disease
     - Small sample size: each variant is valuable
     - Pancreatic cancer is rare: causal variants may be rare
     - Common-variant + rare-variant mixed model
   Standard GWAS uses --maf 0.05 (assumes common-disease common-variant model)
   PDAC context requires retention of rare variants for power

================================================================================
POWER AND DOWNSTREAM IMPLICATIONS
================================================================================

Sample Loss Impact (~$((NSAMP_RAW - NSAMP_08)) samples removed, $((100 * (NSAMP_RAW - NSAMP_08) / NSAMP_RAW))% loss):
  - Reduces effective sample size for association tests
  - For EUR subset: from ~762 to ${NSAMP_08} samples available for analysis
  - Among cases (EUR): ~$((254 - (NSAMP_RAW - NSAMP_08))) cases remain
  - **Power impact:** Lower statistical power to detect effects of given size
  - Mitigation: Rare variants (retained via low MAF) may compensate

Variant Loss Impact (~$((NVAR_RAW - NVAR_08)) variants removed, $((100 * (NVAR_RAW - NVAR_08) / NVAR_RAW))% loss):
  - Most loss from variant call rate + HWE: removes noisy/unreliable variants
  - MAF filter removes ~$((NVAR_07 - NVAR_08)) common variants
  - **Result:** Higher-quality variant set for association testing
  - Remaining variants are reliable for downstream imputation and analysis

Recommendations:
  - Confirm final sample/variant counts match expected ranges (see above)
  - Document all threshold choices for manuscript Methods section
  - Keep this summary in supplementary materials for reproducibility
  - If power is insufficient post-association, consider: meta-analysis, imputation

================================================================================
NEXT STEPS
================================================================================

1. Verify final dataset files exist:
   - ${DATASET_NAME}_08_filt.bed/bim/fam

2. Population stratification analysis (Section 2):
   - PCA on final clean genotypes to identify ancestry subgroups
   - Exclude non-EUR individuals for EUR-specific GWAS

3. Imputation (Section 3, optional):
   - Prepare pre-imputation QC-pass dataset
   - Use tools like EAGLE (phasing) + minimac (imputation)

4. Association testing (Section 4):
   - Run logistic/linear regression on final QC-pass genotypes
   - Adjust for principal components (PC1-10) and covariates

================================================================================
EOF

cat "${OUT_DIR}/${DATASET_NAME}_qc_summary.txt"
echo ""
echo "✓ Summary saved to: ${OUT_DIR}/${DATASET_NAME}_qc_summary.txt"

# ============================================================================
# STEP 3: Create R-based visualizations
# ============================================================================
echo ""
echo "=== Creating visualizations (R) ==="
echo ""

COUNTS_FILE="${OUT_DIR}/${DATASET_NAME}_qc_counts.tsv"
cat > "$COUNTS_FILE" << EOF
step	samples	variants
Raw	${NSAMP_RAW}	${NVAR_RAW}
Sample Call Rate	${NSAMP_02}	${NVAR_02}
Sex Check	${NSAMP_03}	${NVAR_03}
Heterozygosity	${NSAMP_04}	${NVAR_04}
Variant Call Rate	${NSAMP_05}	${NVAR_05}
Hardy-Weinberg	${NSAMP_06}	${NVAR_06}
Relatedness	${NSAMP_07}	${NVAR_07}
MAF Filter	${NSAMP_08}	${NVAR_08}
EOF

cat > "${OUT_DIR}/${DATASET_NAME}_qc_visualizations.R" << 'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
dataset_name <- args[1]
seed <- as.numeric(args[2])
counts_file <- args[3]
out_dir <- args[4]

set.seed(seed)

counts <- read.table(counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
steps <- gsub(" ", "\n", counts$step)
samples <- counts$samples
variants <- counts$variants

# Create combined visualization
pdf(file = file.path(out_dir, paste0(dataset_name, "_qc_retention_rates.pdf")), width = 12, height = 8)

# Panel 1: Sample retention
par(mfrow = c(1, 2), mar = c(10, 5, 3, 2))

barplot(samples, names.arg = steps, ylim = c(0, max(samples) * 1.1),
        main = "Sample Retention Through QC Pipeline",
        ylab = "Number of Samples", col = "steelblue", border = "white", cex.axis = 0.9, las = 2)
grid(ny = NA, nx = 0, col = "gray80")

# Panel 2: Variant retention
barplot(variants / 1000, names.arg = steps, ylim = c(0, max(variants) / 1000 * 1.1),
        main = "Variant Retention Through QC Pipeline",
        ylab = "Number of Variants (thousands)", col = "forestgreen", border = "white", cex.axis = 0.9, las = 2)
grid(ny = NA, nx = 0, col = "gray80")

dev.off()

cat(sprintf("Retention visualization saved to: %s_qc_retention_rates.pdf\n", dataset_name))

# Create decision tree (text-based for simplicity)
pdf(file = file.path(out_dir, paste0(dataset_name, "_qc_decision_tree.pdf")), width = 11, height = 8.5)

plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, 
     main = "PDAC GWAS QC Decision Tree", xlab = "", ylab = "", cex.main = 1.5)

# Add text describing the pipeline
text_y <- 9.5
text(5, text_y, sprintf("Raw Genotypes (%s samples, %s variants)", samples[1], variants[1]),
     cex = 1, font = 2, adj = 0.5)

steps_text <- c(
  "↓ Sample Call Rate (--mind 0.02): remove samples with >2% missing",
  "↓ Sex Check (chrX F-stat): remove sex-discordant samples",
  "↓ Heterozygosity (F ± 3 SD): remove contamination/outlier samples",
  "↓ Variant Call Rate (--geno 0.05): remove variants with >5% missing",
  "↓ Hardy-Weinberg (p > 1e-6 in controls): remove error-prone variants",
  "↓ Relatedness (PI_HAT > 0.1875): prune related samples",
  "↓ MAF Filter (--maf 0.01): retain rare variants (PDAC strategy)"
)

for (i in seq_along(steps_text)) {
  text_y <- text_y - 0.9
  text(1, text_y, steps_text[i], cex = 0.85, adj = 0)
}

text(5, 0.8, sprintf("Clean Dataset Ready for Association Testing (%s samples, %s variants)", samples[length(samples)], variants[length(variants)]),
     cex = 1, font = 2, adj = 0.5, col = "darkgreen")

dev.off()

cat(sprintf("Decision tree saved to: %s_qc_decision_tree.pdf\n", dataset_name))
EOF

Rscript "${OUT_DIR}/${DATASET_NAME}_qc_visualizations.R" "${DATASET_NAME}" "${SEED}" "$COUNTS_FILE" "$OUT_DIR"

# ============================================================================
# FINAL SUMMARY
# ============================================================================
echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║          GENOTYPING QC PIPELINE COMPLETE                      ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "✓ Final clean dataset:"
echo "  - ${OUT_DIR}/${DATASET_NAME}_08_filt.bed/bim/fam"
echo "  - ${NSAMP_08} samples"
echo "  - ${NVAR_08} variants"
echo ""
echo "✓ Summary and decision tree:"
echo "  - ${OUT_DIR}/${DATASET_NAME}_qc_summary.txt"
echo "  - ${OUT_DIR}/${DATASET_NAME}_qc_retention_rates.pdf"
echo "  - ${OUT_DIR}/${DATASET_NAME}_qc_decision_tree.pdf"
echo ""
echo "Ready for downstream analysis:"
echo "  → Section 2: Population stratification (PCA)"
echo "  → Section 3: Imputation (optional)"
echo "  → Section 4: Association testing"
echo ""
