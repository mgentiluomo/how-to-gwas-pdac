#!/usr/bin/env bash

################################################################################
# Section 7: Linkage disequilibrium — 01: LD decay by ancestry
#
# Measures, within each ancestry group separately, how quickly linkage
# disequilibrium decays with physical distance across the ABO region, and how
# many variants are statistically indistinguishable from the causal one.
#
# WHY THIS IS DONE PER GROUP AND NEVER POOLED
#   LD is a property of a population. Computing r-squared on a pooled
#   multi-ancestry sample measures neither group: allele-frequency differences
#   between groups generate correlations between variants that are not in LD
#   within any of them, so the pooled estimate is inflated by exactly the
#   structure this guide spends Section 4 removing. There is no threshold that
#   repairs this, for the same reason there is none for the pooled
#   Hardy-Weinberg test.
#
# WHY A COMMON MAF FILTER IS APPLIED FIRST
#   r-squared is bounded by allele frequency: two variants with very different
#   frequencies cannot reach r-squared = 1 even when perfectly associated. A
#   group-specific MAF filter would therefore compare different variant sets
#   between groups and confound frequency with LD. The filter here is applied
#   at 5% within each group, and the number of variants retained per group is
#   reported so that the comparison can be judged.
#
# INPUT   results/pca/pdac_demo_02_hwe_filt.bed/bim/fam  (the analysis-ready set)
#         demo_data/sample_ancestry.tsv
# OUTPUT  pdac_demo_07_ld_decay.tsv        mean r2 by distance bin and group
#         pdac_demo_07_ld_causal.tsv       variants in LD with the causal one
#         pdac_demo_07_ld_decay.png        the decay curves
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "$SCRIPT_DIR/../../scripts/dev" ]; then
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
else
  PROJECT_ROOT="$(pwd)"
fi
cd "$PROJECT_ROOT"

DATA="${1:-results/pca/pdac_demo_02_hwe_filt}"
ANC="${2:-demo_data/sample_ancestry.tsv}"
OUT="${3:-results/ld}"

# The region: the causal locus plus 400 kb on each side.
CHR=9
FROM=132900000
TO=133700000
CAUSAL="9:133273682:A:T"

WINDOW_KB=300     # widest pair distance considered
MAF=0.05          # applied identically in every group, see header

mkdir -p "$OUT"
GROUPS=$(awk '{print $2}' "$ANC" | sort -u | grep -v '^$')

echo ">>> Region: chr${CHR}:${FROM}-${TO}, causal variant ${CAUSAL}"
echo ">>> Groups: $(echo "$GROUPS" | tr '\n' ' ')"
echo ""

for G in $GROUPS; do
    awk -v g="$G" '$2 == g { print $1"\t"$1 }' "$ANC" > "${OUT}/_keep_${G}.txt"

    plink2 --bfile "$DATA" \
           --keep "${OUT}/_keep_${G}.txt" \
           --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
           --maf "$MAF" \
           --make-bed --out "${OUT}/_region_${G}" > /dev/null

    # --r2-unphased: r-squared from unphased genotypes, which is what array data
    # gives you. --ld-window-r2 0 keeps every pair, including the uninformative
    # ones, because dropping them would bias the mean upward in exactly the
    # groups where LD is lowest.
    plink2 --bfile "${OUT}/_region_${G}" \
           --r2-unphased \
           --ld-window-kb "$WINDOW_KB" \
           --ld-window 999999 \
           --ld-window-r2 0 \
           --out "${OUT}/_ld_${G}" > /dev/null

    N=$(wc -l < "${OUT}/_region_${G}.bim")
    echo "    ${G}: ${N} variants with MAF >= ${MAF} in the region"
done

echo ""
echo ">>> Summarising"

Rscript - "$OUT" "$CAUSAL" "$(echo "$GROUPS" | tr '\n' ' ')" <<'RSCRIPT'
args   <- commandArgs(trailingOnly = TRUE)
OUT    <- args[1]; CAUSAL <- args[2]
GROUPS <- strsplit(trimws(args[3]), "\\s+")[[1]]

breaks <- c(0, 10, 25, 50, 100, 200, 300)
labels <- paste0(head(breaks, -1), "-", tail(breaks, -1), " kb")

decay <- NULL; causal <- NULL
for (g in GROUPS) {
  d <- read.table(file.path(OUT, paste0("_ld_", g, ".vcor")),
                  header = TRUE, comment.char = "", check.names = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1])
  d$dist_kb <- abs(d$POS_B - d$POS_A) / 1000
  d$bin <- cut(d$dist_kb, breaks = breaks, labels = labels, right = FALSE)

  m <- tapply(d$UNPHASED_R2, d$bin, mean)
  n <- tapply(d$UNPHASED_R2, d$bin, length)
  decay <- rbind(decay, data.frame(group = g, bin = labels,
                                   mean_r2 = round(as.numeric(m[labels]), 4),
                                   n_pairs = as.numeric(n[labels])))

  # pairs involving the causal variant
  k <- d[d$ID_A == CAUSAL | d$ID_B == CAUSAL, ]
  causal <- rbind(causal, data.frame(
    group        = g,
    pairs_tested = nrow(k),
    r2_above_0.8 = sum(k$UNPHASED_R2 > 0.8),
    r2_above_0.5 = sum(k$UNPHASED_R2 > 0.5),
    max_r2       = round(max(k$UNPHASED_R2), 3)))
}

write.table(decay,  file.path(OUT, "pdac_demo_07_ld_decay.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(causal, file.path(OUT, "pdac_demo_07_ld_causal.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nMean r-squared by distance\n")
w <- reshape(decay[, c("group", "bin", "mean_r2")], idvar = "bin",
             timevar = "group", direction = "wide")
print(w, row.names = FALSE)

cat("\nVariants in LD with the causal variant\n")
print(causal, row.names = FALSE)

png(file.path(OUT, "pdac_demo_07_ld_decay.png"),
    width = 1400, height = 1000, res = 160)
mid <- (head(breaks, -1) + tail(breaks, -1)) / 2
cols <- c("#1f4e79", "#c55a11", "#2e7d32")
plot(NA, xlim = c(0, max(mid)), ylim = c(0, max(decay$mean_r2) * 1.05),
     xlab = "Distance between variants (kb)",
     ylab = expression(paste("mean ", r^2)),
     main = "LD decay in the ABO region, by ancestry group")
for (i in seq_along(GROUPS)) {
  s <- decay[decay$group == GROUPS[i], ]
  lines(mid, s$mean_r2, col = cols[i], lwd = 2, type = "b", pch = 19)
}
legend("topright", toupper(GROUPS), col = cols[seq_along(GROUPS)],
       lwd = 2, pch = 19, bty = "n")
invisible(dev.off())
RSCRIPT

rm -f "${OUT}"/_keep_*.txt "${OUT}"/_region_* "${OUT}"/_ld_*

echo ""
echo ">>> Written to ${OUT}/"
echo "[NEXT] The fine-mapping section uses these numbers to explain its result."
