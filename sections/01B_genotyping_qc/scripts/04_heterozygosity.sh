#!/bin/bash

################################################################################
# Section 1B: Genotyping QC — Step 04: Heterozygosity & Inbreeding Outliers
# 
# PURPOSE:
#   Remove samples with extreme heterozygosity that indicate:
#   - Sample contamination (unusually high heterozygosity)
#   - Inbreeding or population substructure (unusually low heterozygosity)
#   - Quality outliers from hidden stratification
#
# INPUT:
#   - pdac_demo_03_filt.bed/bim/fam (from Step 03, sex-checked)
#
# OUTPUT:
#   - pdac_demo.04_het.het — heterozygosity (F statistic) per sample
#   - pdac_demo_04_het_outliers.txt — list of outliers to remove
#   - pdac_demo_04_het_outliers.pdf — visualization (R)
#
# METHOD:
#   1. Compute autosomal heterozygosity for each sample
#   2. Identify outliers as samples with F > mean(F) + 3*SD or < mean(F) - 3*SD
#   3. This 3-SD threshold catches ~0.27% of samples in normal distribution
#   4. Visualize distribution and flag samples
#
# NOTES:
#   - Applied after sex check (which uses X chromosome)
#   - Uses autosomes only (more stable, larger SNP number)
#   - Combines with later PCA (Section 2) for population-level stratification
#
################################################################################

set -e

# Configuration
SEED="${1:-2026}"
DATASET_INPUT="${2:-.}"
DATASET_NAME="pdac_demo"

# ============================================================================
# STEP 1: Compute heterozygosity
# ============================================================================
echo ""
echo "=== Computing heterozygosity (autosomal F statistic) ==="
echo ""

# Why --het?
#   Heterozygosity (F) reflects individual-level genetic patterns:
#   - F = inbreeding coefficient
#   - Positive F: fewer heterozygotes (inbreeding, population substructure)
#   - Negative F: more heterozygotes (outbreeding, admixture, contamination)
#
# Why autosomes only?
#   More variants (430K on autosomes vs 6K on our simulated chrX)
#   More stable and less sensitive to single variants
#   Sex-specific heterozygosity patterns confound X analysis

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_03_filt \
  --het \
  --out ${DATASET_NAME}_04_het

echo "✓ Heterozygosity computed: ${DATASET_NAME}_04_het.het"

# ============================================================================
# STEP 2: Identify outliers using R
# ============================================================================
echo ""
echo "=== Identifying heterozygosity outliers (R) ==="
echo ""

# Create R script to detect outliers and visualize
cat > ${DATASET_NAME}_04_het_outliers.R << 'EOF'
#!/usr/bin/env Rscript

# Configuration
args <- commandArgs(trailingOnly = TRUE)
het_file <- args[1]
output_prefix <- args[2]
seed <- as.numeric(args[3])

set.seed(seed)

# Read heterozygosity data
het <- read.table(het_file, header = TRUE, stringsAsFactors = FALSE)

# Calculate F statistics
# NOTE: Het column represents observed heterozygosity
# We compute inbreeding coefficient F manually if needed, or use PLINK2's output
het$F <- 1 - (het$HET / (het$N_NM))  # het$HET = observed HET, het$N_NM = # non-missing

# Identify outliers: mean(F) ± 3 * SD(F)
mean_f <- mean(het$F, na.rm = TRUE)
sd_f <- sd(het$F, na.rm = TRUE)
lower_bound <- mean_f - 3 * sd_f
upper_bound <- mean_f + 3 * sd_f

# Flag outliers
het$outlier <- het$F < lower_bound | het$F > upper_bound

# Summary statistics
n_total <- nrow(het)
n_outliers <- sum(het$outlier, na.rm = TRUE)
percent_outliers <- 100 * n_outliers / n_total

cat(sprintf("Heterozygosity Summary:\n"))
cat(sprintf("  Mean F: %.4f\n", mean_f))
cat(sprintf("  SD F:   %.4f\n", sd_f))
cat(sprintf("  Lower outlier bound (F < %.4f): %.4f\n", lower_bound, lower_bound))
cat(sprintf("  Upper outlier bound (F > %.4f): %.4f\n", upper_bound, upper_bound))
cat(sprintf("  Total samples: %d\n", n_total))
cat(sprintf("  Outliers (3-SD): %d (%.2f%%)\n\n", n_outliers, percent_outliers))

# Write outlier list (for removal)
outliers <- het[het$outlier, c("FID", "IID")]
write.table(
  outliers,
  file = paste0(output_prefix, "_outliers.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = " "
)
cat(sprintf("Outlier list written to: %s_outliers.txt\n\n", output_prefix))

# Create visualization: histogram + density + outliers highlighted
pdf(file = paste0(output_prefix, "_outliers.pdf"), width = 10, height = 6)

hist(het$F, breaks = 50, main = "Distribution of Heterozygosity (F statistic)",
     xlab = "Inbreeding Coefficient (F)", ylab = "Frequency", col = "lightblue", border = "white")
abline(v = mean_f, col = "blue", lwd = 2, lty = 1, label = "Mean")
abline(v = lower_bound, col = "red", lwd = 2, lty = 2, label = "Outlier threshold")
abline(v = upper_bound, col = "red", lwd = 2, lty = 2)
legend("topright", c("Mean", "Outlier bounds (±3 SD)"), 
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# Q-Q plot
qqnorm(het$F, main = "Q-Q Plot: Heterozygosity")
qqline(het$F, col = "red")

dev.off()
cat(sprintf("Plot saved to: %s_outliers.pdf\n", output_prefix))
EOF

# Run R script
Rscript ${DATASET_NAME}_04_het_outliers.R \
  ${DATASET_NAME}_04_het.het \
  ${DATASET_NAME}_04_het \
  ${SEED}

# ============================================================================
# STEP 3: Create dataset without outliers
# ============================================================================
echo ""
echo "=== Removing heterozygosity outliers ==="
echo ""

plink2 \
  --bfile ${DATASET_INPUT}/${DATASET_NAME}_03_filt \
  --remove ${DATASET_NAME}_04_het_outliers.txt \
  --make-bed \
  --out ${DATASET_NAME}_04_filt

echo "✓ Het-filtered dataset: ${DATASET_NAME}_04_filt.bed/bim/fam"

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NSAMP_BEFORE=$(tail -n +2 ${DATASET_INPUT}/${DATASET_NAME}_03_filt.fam | wc -l)
NSAMP_AFTER=$(tail -n +2 ${DATASET_NAME}_04_filt.fam | wc -l)
NSAMP_REMOVED=$((NSAMP_BEFORE - NSAMP_AFTER))

echo "Samples before het filter:   ${NSAMP_BEFORE}"
echo "Samples after het filter:    ${NSAMP_AFTER}"
echo "Het outliers removed:        ${NSAMP_REMOVED}"
echo ""
echo "Visualization: ${DATASET_NAME}_04_het_outliers.pdf"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the variant call rate filter:"
echo ""
echo "  bash 05_variant_callrate.sh"
echo ""
echo "This will:"
echo "  - Filter variants by genotyping call rate (--geno 0.05)"
echo "  - Remove variants with >5% missing genotypes"
echo ""
