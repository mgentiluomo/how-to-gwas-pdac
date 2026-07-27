################################################################################
# Section 4A: Association testing — Step 05: results summary and regional plot
#
# PURPOSE:
#   Reduce a genome-wide scan to the small number of things that are actually
#   reported: the analysis set, the inflation factor, the variants above each
#   threshold, and a regional view of the top locus.
#
#   The regional plot is the one that tells you whether a peak is real. A
#   genuine association appears as a cluster: the lead variant surrounded by
#   correlated neighbours with progressively weaker signal. A single isolated
#   point with nothing around it is usually a genotyping artefact, and should be
#   inspected in the cluster plots before it is believed.
#
# INPUT:
#   - results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid
#   - results/assoc/pdac_demo_04A_lambda.tsv
#
# OUTPUT:
#   - results/assoc/pdac_demo_04A_top_hits.tsv
#   - results/assoc/pdac_demo_04A_regional.png
#   - results/assoc/pdac_demo_04A_summary.txt
################################################################################

out_dir <- "results/assoc"
window  <- 250000   # base pairs either side of the lead variant

f <- Sys.glob(file.path(out_dir, "pdac_demo_04A_gwas.*.glm.logistic.hybrid"))
if (!length(f)) f <- Sys.glob(file.path(out_dir, "pdac_demo_04A_gwas.*.glm.logistic"))
if (!length(f)) stop("No association output found in ", out_dir)

d <- read.table(f[1], header = TRUE, comment.char = "", check.names = FALSE)
names(d)[1] <- sub("^#", "", names(d)[1])
if ("TEST" %in% names(d)) d <- d[d$TEST == "ADD", ]
d <- d[is.finite(d$P) & d$P > 0, ]
d <- d[order(d$P), ]

# Odds ratio and confidence interval from the log-odds scale.
if (!"OR" %in% names(d) && "BETA" %in% names(d)) d$OR <- exp(d$BETA)
if ("BETA" %in% names(d) && "SE" %in% names(d)) {
  d$OR_L95 <- exp(d$BETA - 1.96 * d$SE)
  d$OR_U95 <- exp(d$BETA + 1.96 * d$SE)
}

keep_cols <- intersect(c("CHROM", "POS", "ID", "A1", "A1_FREQ", "OBS_CT",
                         "OR", "OR_L95", "OR_U95", "P", "FIRTH?", "ERRCODE"),
                       names(d))
top <- head(d[, keep_cols], 25)
write.table(top, file.path(out_dir, "pdac_demo_04A_top_hits.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- regional plot around the lead variant -----------------------------------
lead <- d[1, ]
reg  <- d[d$CHROM == lead$CHROM &
          d$POS > lead$POS - window & d$POS < lead$POS + window, ]

png(file.path(out_dir, "pdac_demo_04A_regional.png"),
    width = 1700, height = 1200, res = 200)
plot(reg$POS / 1e6, -log10(reg$P), pch = 19, cex = 0.8, col = "#4B6584",
     xlab = sprintf("Position on chromosome %s (Mb)", lead$CHROM),
     ylab = expression(-log[10](italic(P))),
     main = sprintf("Regional association around %s", lead$ID))
points(lead$POS / 1e6, -log10(lead$P), pch = 18, cex = 1.8, col = "#D55E00")
abline(h = -log10(5e-8), col = "#D55E00", lty = 2)
text(lead$POS / 1e6, -log10(lead$P), labels = lead$ID, pos = 4, cex = 0.7)
dev.off()

# --- text summary ------------------------------------------------------------
lam_file <- file.path(out_dir, "pdac_demo_04A_lambda.tsv")
lam <- if (file.exists(lam_file)) read.table(lam_file, header = TRUE, sep = "\t") else NULL
getv <- function(k) if (!is.null(lam) && k %in% lam$quantity) lam$value[lam$quantity == k] else NA

con <- file(file.path(out_dir, "pdac_demo_04A_summary.txt"), open = "wt")
w <- function(...) cat(..., "\n", sep = "", file = con)

w("Association summary, European analysis set")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("")
w("Cases: ", getv("cases"), "   Controls: ", getv("controls"))
w("Variants tested: ", nrow(d))
w("lambda GC: ", getv("lambda_gc"), "   lambda_1000: ", getv("lambda_1000"))
w("")
w("P < 5e-8: ", sum(d$P < 5e-8), " variants")
w("P < 1e-5: ", sum(d$P < 1e-5), " variants")
w("")
w("Lead variant: ", lead$ID)
w("  effect allele: ", lead$A1, "   frequency: ", sprintf("%.4f", lead$A1_FREQ))
w("  OR: ", sprintf("%.3f", lead$OR),
  " (95% CI ", sprintf("%.3f", lead$OR_L95), " to ", sprintf("%.3f", lead$OR_U95), ")")
w("  P:  ", format(lead$P, digits = 4))
w("")
w("Variants within ", window / 1000, " kb of the lead variant: ", nrow(reg))
close(con)

cat("Lead variant:", lead$ID, "\n")
cat("  OR", sprintf("%.3f", lead$OR),
    sprintf("(95%% CI %.3f-%.3f)", lead$OR_L95, lead$OR_U95),
    " P =", format(lead$P, digits = 3), "\n")
cat("\nWritten:\n")
cat("  ", file.path(out_dir, "pdac_demo_04A_top_hits.tsv"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_04A_regional.png"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_04A_summary.txt"), "\n")
cat("\nA note on the effect estimate. In this demonstration the true simulated\n")
cat("odds ratio at the causal variant is 2.40. The value recovered above is\n")
cat("larger. That is the winner's curse: reading an effect size from the same\n")
cat("underpowered scan that discovered it selects for upward fluctuations.\n")
cat("Effect sizes for follow-up should come from replication or from a\n")
cat("combined analysis, never from discovery alone.\n\n")
