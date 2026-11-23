#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 03: Sex Check
# 
# PURPOSE:
#   Identify and flag samples with sex discordance (genetic sex does not match
#   reported sex). This catches:
#   - Sample swaps
#   - Data entry errors
#   - Potential contamination (XY in supposed females, XX in males)
#
# INPUT:
#   - pdac_demo_02_filt.bed/bim/fam (from Step 02, sample call rate filtered)
#
# OUTPUT:
#   - pdac_demo_03_sexcheck.sexcheck — sex concordance report
#   - pdac_demo_03_sexcheck_discordant.txt — list of discordant samples to remove
#   - pdac_demo_03_sexcheck_summary.tsv — thresholds, counts by reported sex,
#     and the samples between the two thresholds
#   - pdac_demo_03_sexcheck_F.png — distribution of X F by reported sex,
#     with the two thresholds drawn
#
# METHOD:
#   Compute X chromosome heterozygosity (F statistic on chrX).
#   - Females (XX): low heterozygosity (F near 0)
#   - Males (XY): high heterozygosity (F near 1)
#   Samples with F inconsistent with reported sex are flagged.
#
# NOTES:
#   - Requires sex information in .fam file (column 4: 1=M, 2=F)
#   - Our demo dataset has simulated chrX; real data uses genotyped chrX
#   - Thresholds are passed explicitly: max-female-xf 0.2, min-male-xf 0.8.
#     These are the PLINK 1.9 defaults, a CONVENTION, not an evaluated
#     threshold. PLINK 2 has no defaults: without them it warns, sets
#     min-male-xf to 1 and no female threshold, and asks for thresholds read
#     from the data. The script therefore checks the convention against the
#     observed distribution: it counts the samples falling between 0.2 and
#     0.8, and draws the distribution by reported sex with both lines. Where
#     samples sit between the thresholds, read the plot and decide before
#     excluding anyone.
#   - The pseudo-autosomal regions are split off first (--split-par hg38):
#     left on X, they are diploid in males and pull the male F towards the
#     female range. In the demo data all X variants are outside the PARs, so
#     the split moves none; the step is kept because real data need it.
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

mkdir -p "$OUT_DIR"

# ============================================================================
# STEP 1: Compute X-chromosome genotyping frequency (F statistic)
# ============================================================================
echo ""
echo "=== Computing X chromosome heterozygosity (sex check) ==="
echo ""

# Why --check-sex?
#   PLINK2 computes the inbreeding coefficient F on the X chromosome only.
#   - Females (XX) genotypes: F near 0 (mostly heterozygous relative to population)
#   - Males (XY) genotypes: F near 1 (all X alleles homozygous, haploid)
#
# The rationale:
#   Males have one X chromosome, so all X markers are "homozygous" (haploid).
#   Females have two X chromosomes, so they can be heterozygous.
#   If a sample labeled "female" shows high F (close to 1), it's likely mislabeled
#   or swapped. Conversely, a sample labeled "male" with low F is suspicious.

MAX_FEMALE_XF="0.2"   # convention: PLINK 1.9 default
MIN_MALE_XF="0.8"     # convention: PLINK 1.9 default

# Split the pseudo-autosomal regions off X (GRCh38 boundaries 2,781,479 and
# 155,701,383) before computing F on X.
plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_02_filt" \
  --split-par hg38 \
  --make-bed \
  --out "${OUT_DIR}/${DATASET_NAME}_03_splitpar"
grep -i "split-par" "${OUT_DIR}/${DATASET_NAME}_03_splitpar.log" | grep -iv "^  --" || true

plink2 \
  --bfile "${OUT_DIR}/${DATASET_NAME}_03_splitpar" \
  --check-sex max-female-xf="${MAX_FEMALE_XF}" min-male-xf="${MIN_MALE_XF}" \
  --out "${OUT_DIR}/${DATASET_NAME}_03_sexcheck"

echo "✓ Sex check complete: ${OUT_DIR}/${DATASET_NAME}_03_sexcheck.sexcheck"

# ============================================================================
# STEP 2: Identify discordant samples
# ============================================================================
echo ""
echo "=== Identifying sex-discordant samples ==="
echo ""

# Extract discordant samples (STATUS = "PROBLEM"): reported sex differs from
# the sex called with the thresholds above, or F lies between them (no call).

awk '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      gsub(/^#/, "", $i)
      if ($i == "FID") fid = i
      if ($i == "IID") iid = i
      if ($i == "STATUS") status = i
    }
    next
  }
  status && $status == "PROBLEM" { print $fid, $iid }
' "${OUT_DIR}/${DATASET_NAME}_03_sexcheck.sexcheck" > "${OUT_DIR}/${DATASET_NAME}_03_sexcheck_discordant.txt"

NDISCORDANT=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_03_sexcheck_discordant.txt")

echo "Discordant samples identified: ${NDISCORDANT}"

if [ $NDISCORDANT -gt 0 ]; then
  echo ""
  echo "Samples to review:"
  head -20 "${OUT_DIR}/${DATASET_NAME}_03_sexcheck_discordant.txt"
fi

# ============================================================================
# STEP 3: Create dataset without discordant samples
# ============================================================================
echo ""
echo "=== Removing discordant samples ==="
echo ""

if [ -s "${OUT_DIR}/${DATASET_NAME}_03_sexcheck_discordant.txt" ]; then
  plink2 \
    --bfile "${OUT_DIR}/${DATASET_NAME}_03_splitpar" \
    --remove "${OUT_DIR}/${DATASET_NAME}_03_sexcheck_discordant.txt" \
    --make-bed \
    --out "${OUT_DIR}/${DATASET_NAME}_03_filt"
else
  plink2 \
    --bfile "${OUT_DIR}/${DATASET_NAME}_03_splitpar" \
    --make-bed \
    --out "${OUT_DIR}/${DATASET_NAME}_03_filt"
fi

echo "✓ Sex-filtered dataset: ${OUT_DIR}/${DATASET_NAME}_03_filt.bed/bim/fam"
rm -f "${OUT_DIR}/${DATASET_NAME}_03_splitpar".{bed,bim,fam}

# ============================================================================
# STEP 4: Check the conventional thresholds against the observed distribution
# ============================================================================
echo ""
echo "=== Thresholds against the observed X F distribution ==="
echo ""
Rscript --vanilla - "${OUT_DIR}/${DATASET_NAME}_03_sexcheck.sexcheck" \
  "${OUT_DIR}/${DATASET_NAME}_03_sexcheck" "$MAX_FEMALE_XF" "$MIN_MALE_XF" <<'RSCRIPT'
a <- commandArgs(trailingOnly = TRUE)
d <- read.table(a[1], header = TRUE, comment.char = "", check.names = FALSE)
names(d)[1] <- sub("^#", "", names(d)[1])
lo <- as.numeric(a[3]); hi <- as.numeric(a[4])
rep_sex <- ifelse(d$PEDSEX == 1, "male", ifelse(d$PEDSEX == 2, "female", "unknown"))
q <- function(x) if (length(x)) sprintf("%.3f", x) else "NA"
s <- data.frame(
  quantity = c("max_female_xf (convention, PLINK 1.9 default)", "min_male_xf (convention, PLINK 1.9 default)",
               "reported males", "reported females",
               "reported males with F >= min_male_xf", "reported females with F <= max_female_xf",
               "samples with F between the thresholds (no call)",
               "lowest F among reported males called male", "highest F among reported females called female",
               "discordant (STATUS PROBLEM)"),
  value = c(lo, hi, sum(rep_sex == "male"), sum(rep_sex == "female"),
            sum(rep_sex == "male" & d$F >= hi), sum(rep_sex == "female" & d$F <= lo),
            sum(d$F > lo & d$F < hi),
            q(min(d$F[rep_sex == "male" & d$SNPSEX == 1])), q(max(d$F[rep_sex == "female" & d$SNPSEX == 2])),
            sum(d$STATUS == "PROBLEM")))
write.table(s, paste0(a[2], "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(s, row.names = FALSE)
png(paste0(a[2], "_F.png"), width = 1600, height = 1000, res = 200)
br <- seq(min(-0.2, min(d$F)), max(1.05, max(d$F)) + 0.05, by = 0.025)
hf <- hist(d$F[rep_sex == "female"], breaks = br, plot = FALSE)
hm <- hist(d$F[rep_sex == "male"], breaks = br, plot = FALSE)
ymax <- max(c(hf$counts, hm$counts, 1))
plot(hf, col = adjustcolor("#D55E00", 0.6), border = NA, ylim = c(0, ymax),
     xlab = "X-chromosome inbreeding coefficient F", ylab = "Samples",
     main = "Sex check: X F by reported sex")
plot(hm, col = adjustcolor("#0072B2", 0.6), border = NA, add = TRUE)
abline(v = c(lo, hi), lty = 2)
text(c(lo, hi), ymax * 0.95, c(paste("max female", lo), paste("min male", hi)), pos = c(4, 2), cex = 0.8)
pr <- d[d$STATUS == "PROBLEM", ]
if (nrow(pr)) rug(pr$F, col = "black", lwd = 2)
legend(x = 0.32, y = ymax * 0.6, c("reported female", "reported male", "discordant (rug)"),
       fill = c(adjustcolor("#D55E00", 0.6), adjustcolor("#0072B2", 0.6), NA),
       border = NA, lty = c(NA, NA, 1), bty = "n", cex = 0.8)
invisible(dev.off())
n_mid <- sum(d$F > lo & d$F < hi)
cat(if (n_mid == 0) "\nNo sample lies between the thresholds: on these data the convention and the\ndistribution agree, and the calls do not depend on where inside the gap the lines sit.\n"
    else sprintf("\n%d sample(s) lie between the thresholds: inspect the plot before excluding.\n", n_mid))
RSCRIPT

# ============================================================================
# SUMMARY & NEXT STEP
# ============================================================================
echo ""
echo "=== Summary ==="
echo ""

NSAMP_BEFORE=$(wc -l < "${DATASET_INPUT}/${DATASET_NAME}_02_filt.fam")
NSAMP_AFTER=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_03_filt.fam")
NSAMP_REMOVED=$((NSAMP_BEFORE - NSAMP_AFTER))

echo "Samples before sex check:      ${NSAMP_BEFORE}"
echo "Samples after sex check:       ${NSAMP_AFTER}"
echo "Sex-discordant samples removed: ${NSAMP_REMOVED}"
echo ""

# ============================================================================
# NEXT STEP
# ============================================================================
echo "=== NEXT STEP ==="
echo ""
echo "Run the heterozygosity check:"
echo ""
echo "  bash scripts/01B_genotyping_qc/04_heterozygosity.sh"
echo ""
echo "This will:"
echo "  - Compute per-sample heterozygosity (F coefficient on autosomes)"
echo "  - Identify outliers (contamination, inbreeding, population outliers)"
echo "  - Generate visualization plots"
echo ""
