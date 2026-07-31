#!/usr/bin/env Rscript

################################################################################
# Section 6: Fine mapping — Step 03: seeing the credible set
#
# PURPOSE:
#   Two figures. The regional plot shows where the signal is and how the
#   posterior probability is distributed across it; the comparison shows what
#   combining ancestries does to the resolution.
#
#   A resolved locus looks like one or two variants carrying most of the
#   probability. An unresolved one looks like probability spread thinly over
#   tens of variants, and no method will improve on that: the information is not
#   in the data.
#
# INPUT:
#   - results/finemap/pdac_demo_06_pp_eur.tsv
#   - results/finemap/pdac_demo_06_pp_meta.tsv   (optional)
#
# OUTPUT:
#   - results/finemap/pdac_demo_06_regional.png
#   - results/finemap/pdac_demo_06_credible_comparison.png
################################################################################

out_dir <- "results/finemap"
# The generative truth is read from the released truth table rather than
# hardcoded. A constant written here goes stale the moment the phenotype is
# regenerated, and a stale constant of exactly this kind is why this dataset
# was rebuilt.
.truth <- read.delim("demo_data/truth.tsv", stringsAsFactors = FALSE)
.truth <- .truth[grepl(":", .truth$variant), ]
TRUTH   <- .truth$variant[1]
blue <- "#0072B2"; orange <- "#D55E00"; green <- "#009E73"; grey <- "#4B6584"

e <- read.delim(file.path(out_dir, "pdac_demo_06_pp_eur.tsv"), stringsAsFactors = FALSE)
has_meta <- file.exists(file.path(out_dir, "pdac_demo_06_pp_meta.tsv"))
if (has_meta) m <- read.delim(file.path(out_dir, "pdac_demo_06_pp_meta.tsv"),
                              stringsAsFactors = FALSE)

# --- regional plot: association above, posterior below ------------------------
png(file.path(out_dir, "pdac_demo_06_regional.png"),
    width = 1800, height = 1600, res = 200)
par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 1), mgp = c(2.5, 0.7, 0))

col_e <- ifelse(e$in_cs == "TRUE" | e$in_cs == TRUE, orange, grey)
plot(e$POS / 1e6, -log10(e$P), pch = 19, cex = 0.8, col = col_e,
     xlab = "", ylab = expression(-log[10](italic(P))),
     main = "Association in the European analysis set")
abline(h = -log10(5e-8), col = orange, lty = 2)
tr <- e[e$ID == TRUTH, ]
if (nrow(tr)) points(tr$POS / 1e6, -log10(tr$P), pch = 5, cex = 2, lwd = 2)
legend("topleft", bty = "n", cex = 0.75, pch = c(19, 19, 5),
       col = c(orange, grey, "black"),
       legend = c("in the 95% credible set", "outside it", "simulated causal variant"))

plot(e$POS / 1e6, e$PP, type = "h", lwd = 2, col = col_e,
     xlab = "Position on chromosome 9 (Mb)",
     ylab = "Posterior probability of being causal",
     main = "Where the probability actually sits")
if (nrow(tr)) points(tr$POS / 1e6, tr$PP, pch = 5, cex = 2, lwd = 2)
dev.off()

# --- what combining ancestries does ------------------------------------------
if (has_meta) {
  ord_e <- sort(e$PP, decreasing = TRUE)
  ord_m <- sort(m$PP, decreasing = TRUE)
  png(file.path(out_dir, "pdac_demo_06_credible_comparison.png"),
      width = 1700, height = 1000, res = 200)
  par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 1), mgp = c(2.5, 0.7, 0))

  plot(cumsum(ord_e), type = "s", lwd = 2, col = orange, xlim = c(0, 80),
       ylim = c(0, 1), xlab = "Variants, ordered by posterior probability",
       ylab = "Cumulative posterior probability",
       main = "How many variants to reach 95%")
  lines(cumsum(ord_m), type = "s", lwd = 2, col = green)
  abline(h = 0.95, col = "grey50", lty = 3)
  n_e <- sum(e$in_cs == "TRUE" | e$in_cs == TRUE)
  n_m <- sum(m$in_cs == "TRUE" | m$in_cs == TRUE)
  legend("bottomright", bty = "n", cex = 0.75, lwd = 2, col = c(orange, green),
         legend = c(sprintf("European set: %d variants", n_e),
                    sprintf("three strata: %d variant(s)", n_m)))

  b <- barplot(c(max(e$PP), max(m$PP)), names.arg = c("European", "Three strata"),
               col = c(orange, green), ylim = c(0, 1.05),
               ylab = "Posterior probability of the top variant",
               main = "Confidence in a single variant")
  text(b, c(max(e$PP), max(m$PP)), sprintf("%.3f", c(max(e$PP), max(m$PP))),
       pos = 3, cex = 0.8)
  dev.off()
}

cat("Written:\n  ", file.path(out_dir, "pdac_demo_06_regional.png"), "\n")
if (has_meta) cat("  ", file.path(out_dir, "pdac_demo_06_credible_comparison.png"), "\n")
cat("\nEuropean credible set: ", sum(e$in_cs == "TRUE" | e$in_cs == TRUE),
    " variants, top PP ", sprintf("%.3f", max(e$PP)), "\n", sep = "")
if (has_meta) cat("Three strata:          ", sum(m$in_cs == "TRUE" | m$in_cs == TRUE),
    " variant(s), top PP ", sprintf("%.3f", max(m$PP)), "\n", sep = "")
