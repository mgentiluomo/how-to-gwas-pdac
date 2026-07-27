################################################################################
# Section 5: Meta-analysis — Step 03: combining the strata
#
# PURPOSE:
#   Combine the harmonised per-stratum estimates into one result, and quantify
#   how much the strata disagree.
#
# THE FIXED-EFFECT MODEL:
#   Assumes every stratum estimates the same underlying effect, and that they
#   differ only by sampling error. Each estimate is weighted by its precision,
#   w = 1 / SE^2:
#
#     beta_meta = sum(w_i * beta_i) / sum(w_i)
#     SE_meta   = sqrt(1 / sum(w_i))
#
#   This is inverse-variance weighting, the standard for GWAS meta-analysis when
#   all contributing studies report effects on the same scale, here the log odds
#   ratio. Note what the weighting does: a stratum contributes according to its
#   precision, which depends on its effective sample size and on the allele
#   frequency in that population, not on its headline number of participants.
#
# HETEROGENEITY:
#   Cochran's Q tests whether the differences between strata exceed what
#   sampling error alone would produce. I^2 expresses the same information as
#   the percentage of total variation attributable to real differences rather
#   than chance.
#
#   Heterogeneity is not a nuisance to be smoothed away. Across ancestries it is
#   expected and informative: different linkage-disequilibrium structure means a
#   marker tags a causal variant to different degrees in different populations,
#   so the apparent effect differs even when the underlying biology is identical.
#   A locus with strong heterogeneity is a locus where fine-mapping should be
#   done per ancestry, not pooled.
#
# THE RANDOM-EFFECTS MODEL:
#   Allows the true effect to differ between strata, estimating the between-study
#   variance tau^2 by the DerSimonian-Laird method and widening the confidence
#   interval accordingly. It is more conservative wherever heterogeneity exists.
#   It is reported here alongside the fixed-effect result rather than instead of
#   it, because the comparison is what tells the reader whether the choice of
#   model changes the conclusion.
#
# INPUT:
#   - results/meta/pdac_demo_05_harmonised.tsv.gz
#
# OUTPUT:
#   - results/meta/pdac_demo_05_meta.tsv.gz
#   - results/meta/pdac_demo_05_meta_summary.txt
#   - results/meta/pdac_demo_05_top_hits.tsv
#
# Base R only. In practice most consortia use METAL or GWAMA; the arithmetic is
# written out here so that the reader can see what those tools do.
################################################################################

out_dir <- "results/meta"
strata  <- c("eur", "afr", "eas")

h <- read.table(gzfile(file.path(out_dir, "pdac_demo_05_harmonised.tsv.gz")),
                header = TRUE, sep = "\t", check.names = FALSE,
                stringsAsFactors = FALSE)
cat("Harmonised variants read:", nrow(h), "\n")

B <- as.matrix(h[, paste0("BETA_", strata)])
S <- as.matrix(h[, paste0("SE_",   strata)])

ok <- is.finite(B) & is.finite(S) & S > 0
k  <- rowSums(ok)                       # number of strata contributing

W  <- ifelse(ok, 1 / S^2, 0)
B0 <- ifelse(ok, B, 0)

sumW  <- rowSums(W)
beta  <- rowSums(W * B0) / sumW
se    <- sqrt(1 / sumW)
z     <- beta / se
p     <- 2 * pnorm(-abs(z))

# --- heterogeneity -----------------------------------------------------------
Q  <- rowSums(ifelse(ok, W * (B - beta)^2, 0))
df <- k - 1
p_het <- ifelse(df > 0, pchisq(Q, df = df, lower.tail = FALSE), NA)
I2    <- ifelse(df > 0, pmax(0, 100 * (Q - df) / Q), NA)

# --- random effects, DerSimonian-Laird ---------------------------------------
sumW2 <- rowSums(W^2)
C     <- sumW - sumW2 / sumW
tau2  <- ifelse(df > 0, pmax(0, (Q - df) / C), 0)

Wr    <- ifelse(ok, 1 / (S^2 + tau2), 0)
sumWr <- rowSums(Wr)
beta_r <- rowSums(Wr * B0) / sumWr
se_r   <- sqrt(1 / sumWr)
p_r    <- 2 * pnorm(-abs(beta_r / se_r))

res <- data.frame(
  ID = h$ID, CHROM = h$CHROM, POS = h$POS, A1 = h$A1,
  n_strata = k,
  BETA = beta, SE = se, OR = exp(beta), P = p,
  Q = Q, P_het = p_het, I2 = I2,
  BETA_re = beta_r, SE_re = se_r, P_re = p_r,
  strand_ambiguous = h$strand_ambiguous,
  stringsAsFactors = FALSE
)
for (g in strata) {
  res[[paste0("BETA_", g)]] <- h[[paste0("BETA_", g)]]
  res[[paste0("SE_",   g)]] <- h[[paste0("SE_",   g)]]
  res[[paste0("P_",    g)]] <- h[[paste0("P_",    g)]]
  res[[paste0("FREQ_", g)]] <- h[[paste0("FREQ_", g)]]
}

res <- res[order(res$P), ]

gz <- gzfile(file.path(out_dir, "pdac_demo_05_meta.tsv.gz"), "wt")
write.table(res, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

top_cols <- c("ID", "CHROM", "POS", "A1", "n_strata", "OR", "SE", "P",
              "I2", "P_het", "P_re",
              paste0("P_", strata), paste0("FREQ_", strata))
write.table(head(res[, top_cols], 25),
            file.path(out_dir, "pdac_demo_05_top_hits.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- lambda for the meta-analysis --------------------------------------------
lam <- median(qchisq(1 - res$P, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

# --- per-stratum lambda, for comparison --------------------------------------
lam_stratum <- sapply(strata, function(g) {
  pv <- h[[paste0("P_", g)]]
  pv <- pv[is.finite(pv) & pv > 0]
  median(qchisq(1 - pv, df = 1)) / qchisq(0.5, df = 1)
})

lead <- res[1, ]

con <- file(file.path(out_dir, "pdac_demo_05_meta_summary.txt"), open = "wt")
w <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = con) }

w("Meta-analysis summary")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("")
w("Variants meta-analysed: ", nrow(res))
w("Strata: ", paste(strata, collapse = ", "))
w("")
w("Genomic inflation factor")
for (g in strata) w(sprintf("  %-4s  lambda = %.4f", g, lam_stratum[g]))
w(sprintf("  meta  lambda = %.4f", lam))
w("")
w("Significance counts")
w("  P < 5e-8 (fixed effects):  ", sum(res$P < 5e-8))
w("  P < 1e-5 (fixed effects):  ", sum(res$P < 1e-5))
w("  P < 5e-8 (random effects): ", sum(res$P_re < 5e-8))
w("")
w("Lead variant: ", lead$ID)
w(sprintf("  meta OR %.3f, SE %.4f, P %.3e", lead$OR, lead$SE, lead$P))
w(sprintf("  heterogeneity: I2 = %.1f%%, Q = %.2f, P_het = %.3f",
          lead$I2, lead$Q, lead$P_het))
w(sprintf("  random effects: P = %.3e", lead$P_re))
w("  per stratum:")
for (g in strata) {
  w(sprintf("    %-4s  OR %.3f  P %.3e  freq %.3f",
            g, exp(lead[[paste0("BETA_", g)]]), lead[[paste0("P_", g)]],
            lead[[paste0("FREQ_", g)]]))
}
w("")
w("Heterogeneity across all variants")
w("  I2 > 50%: ", sum(res$I2 > 50, na.rm = TRUE),
  sprintf(" (%.1f%%)", 100 * mean(res$I2 > 50, na.rm = TRUE)))
w("  P_het < 0.05: ", sum(res$P_het < 0.05, na.rm = TRUE),
  sprintf(" (%.1f%%)", 100 * mean(res$P_het < 0.05, na.rm = TRUE)))
w("  Under the null, 5% of variants should have P_het < 0.05 by chance.")
close(con)

cat("\nWritten:\n")
cat("  ", file.path(out_dir, "pdac_demo_05_meta.tsv.gz"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_05_meta_summary.txt"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_05_top_hits.tsv"), "\n")
cat("\n=== NEXT STEP ===\n\n")
cat("  Rscript scripts/05_meta_analysis/04_meta_plots.R\n\n")
