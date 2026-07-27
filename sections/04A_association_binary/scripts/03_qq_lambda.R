################################################################################
# Section 4A: Association testing — Step 03: QQ plot and genomic inflation
#
# PURPOSE:
#   Before looking at any result, check whether the test statistics as a whole
#   behave as they should under the null. This is the step that catches residual
#   population structure, cryptic relatedness and batch effects.
#
#   The quantile-quantile plot compares observed P values with the uniform
#   distribution expected if nothing were associated. Points should follow the
#   diagonal, with a departure only in the extreme upper tail, where true
#   associations live. A curve that lifts off the diagonal along its whole
#   length indicates inflation affecting every variant, which is a property of
#   the analysis, not of biology.
#
# HOW LAMBDA IS COMPUTED:
#   lambda = median(qchisq(1 - P, df = 1)) / qchisq(0.5, df = 1)
#
#   Compute it this way. Do not compute it by sorting test statistics in awk or
#   a spreadsheet: that mistake was made once in the construction of this very
#   dataset and produced a convincing inflation signal that did not exist.
#
# INTERPRETING LAMBDA IN A SMALL STUDY:
#   Lambda scales with sample size under genuine polygenicity, so a value must
#   be read relative to the size of the study. In a few hundred cases it is
#   estimated with considerable uncertainty, and 1.05 can be entirely compatible
#   with a well-controlled analysis. lambda_1000, rescaled to 1,000 cases and
#   1,000 controls, is the comparable figure across studies and is reported here
#   as well.
#
# INPUT:
#   - results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid
#
# OUTPUT:
#   - results/assoc/pdac_demo_04A_qq.png
#   - results/assoc/pdac_demo_04A_lambda.tsv
################################################################################

out_dir <- "results/assoc"

f <- Sys.glob(file.path(out_dir, "pdac_demo_04A_gwas.*.glm.logistic.hybrid"))
if (!length(f)) f <- Sys.glob(file.path(out_dir, "pdac_demo_04A_gwas.*.glm.logistic"))
if (!length(f)) stop("No association output found in ", out_dir)

d <- read.table(f[1], header = TRUE, comment.char = "", check.names = FALSE)
names(d)[1] <- sub("^#", "", names(d)[1])
if ("TEST" %in% names(d)) d <- d[d$TEST == "ADD", ]

# Drop variants PLINK could not test, for example those that became monomorphic
# once the analysis set was restricted to one ancestry group.
n_all <- nrow(d)
d <- d[is.finite(d$P) & d$P > 0, ]
cat("Variants in file:            ", n_all, "\n")
cat("Variants with a valid P:     ", nrow(d), "\n")
cat("Dropped (no valid P value):  ", n_all - nrow(d), "\n")

# --- lambda ------------------------------------------------------------------
chisq  <- qchisq(1 - d$P, df = 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

# Case and control counts, taken from the log written by Step 02.
log_file <- file.path(out_dir, "pdac_demo_04A_gwas.log")
n_case <- NA; n_ctrl <- NA
if (file.exists(log_file)) {
  ln <- grep("cases and .* controls remaining", readLines(log_file), value = TRUE)
  if (length(ln)) {
    num <- as.numeric(regmatches(ln[1], gregexpr("[0-9]+", ln[1]))[[1]])
    if (length(num) >= 2) { n_case <- num[1]; n_ctrl <- num[2] }
  }
}

lambda_1000 <- NA
if (!is.na(n_case) && !is.na(n_ctrl)) {
  lambda_1000 <- 1 + (lambda - 1) * (1 / n_case + 1 / n_ctrl) / (1 / 1000 + 1 / 1000)
}

cat("\nGenomic inflation factor\n")
cat("  lambda:      ", sprintf("%.4f", lambda), "\n")
if (!is.na(lambda_1000)) cat("  lambda_1000: ", sprintf("%.4f", lambda_1000), "\n")

write.table(
  data.frame(quantity = c("variants_tested", "cases", "controls",
                          "lambda_gc", "lambda_1000"),
             value = c(nrow(d), n_case, n_ctrl,
                       round(lambda, 4), round(lambda_1000, 4))),
  file.path(out_dir, "pdac_demo_04A_lambda.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

# --- QQ plot -----------------------------------------------------------------
p   <- sort(d$P)
n   <- length(p)
obs <- -log10(p)
exp <- -log10(ppoints(n))

# Thin the bulk of the distribution: plotting 400,000 overlapping points in the
# uninformative part of the plot produces a large file and no extra information.
keep <- c(seq_len(min(n, 5000)),
          seq(5001, n, by = max(1, floor(n / 20000))))
keep <- unique(pmin(keep, n))

png(file.path(out_dir, "pdac_demo_04A_qq.png"),
    width = 1500, height = 1500, res = 200)
plot(exp[keep], obs[keep], pch = 19, cex = 0.4, col = "#0072B2",
     xlab = expression(Expected~-log[10](italic(P))),
     ylab = expression(Observed~-log[10](italic(P))),
     main = "Quantile-quantile plot")
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", bty = "n",
       legend = c(sprintf("lambda = %.3f", lambda),
                  if (!is.na(lambda_1000)) sprintf("lambda[1000] = %.3f", lambda_1000)))
dev.off()

cat("\nWritten:\n")
cat("  ", file.path(out_dir, "pdac_demo_04A_qq.png"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_04A_lambda.tsv"), "\n")
cat("\n=== NEXT STEP ===\n\n")
cat("  Rscript scripts/04A_association_binary/04_manhattan.R\n\n")
