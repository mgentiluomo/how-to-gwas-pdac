#!/usr/bin/env Rscript

################################################################################
# Section 8, step 08: observed against expected gain when the effective sample
# size changes (relatives retained; the two relatedness-pruning rules)
#
# PURPOSE
#   Section 4C retains the relatives the pipeline prunes and reports the P of
#   the causal ABO variant with and without them; the sensitivity step of
#   Section 4A compares two pruning rules. A change in P between two analyses
#   of one phenotype mixes what the change in sample size buys with the draw
#   of the individuals added or removed. This step puts the expected change
#   beside the observed one.
#
# WHAT IT DOES
#   For two analyses with effective sizes N0 and N1 (4 / (1/cases + 1/controls))
#   and P values P0 and P1 at the causal variant:
#     observed chi-square: the 1-df quantile of each P;
#     expected chi-square at N1 with the effect and the frequency fixed:
#       chi0 * N1 / N0 (the non-centrality is proportional to N_eff), and its P;
#     for nested analyses (relatives retained): the z contributed by the added
#       individuals, (z1 sqrt(N1) - z0 sqrt(N0)) / sqrt(N1 - N0), against the
#       z expected from them at the same effect, z0 sqrt((N1 - N0) / N0).
#   Reads only committed result tables of the canonical run; no genotypes.
#
# INPUT
#   sections/04C_association_other/results/pdac_demo_regenie_comparison.tsv
#   sections/04A_association_binary/results/pdac_demo_04A_sens_relatedness.tsv
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_expected_gain.tsv, pdac_demo_08_expected_gain.txt
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(normalizePath(file.path(script_dir, "..", "..", "..")))
OUT <- "sections/08_replicates/results"
f4c <- "sections/04C_association_other/results/pdac_demo_regenie_comparison.tsv"
f4a <- "sections/04A_association_binary/results/pdac_demo_04A_sens_relatedness.tsv"
a <- read.delim(f4c, check.names = FALSE, stringsAsFactors = FALSE); b <- read.delim(f4a, check.names = FALSE, stringsAsFactors = FALSE)
neff <- function(ca, co) 4 / (1 / ca + 1 / co)
row <- function(comparison, source, ref, alt, ca0, co0, p0, ca1, co1, p1, nested) {
  n0 <- neff(ca0, co0); n1 <- neff(ca1, co1); c0 <- qchisq(p0, 1, lower.tail = FALSE); c1 <- qchisq(p1, 1, lower.tail = FALSE)
  ce <- c0 * n1 / n0; z0 <- sqrt(c0); z1 <- sqrt(c1)
  data.frame(comparison = comparison, source = source, reference = ref, alternative = alt,
             cases_ref = ca0, controls_ref = co0, n_eff_ref = n0, cases_alt = ca1, controls_alt = co1, n_eff_alt = n1,
             n_eff_ratio = n1 / n0, p_ref = p0, p_alt = p1, chisq_ref = c0, chisq_alt_observed = c1, chisq_alt_expected = ce,
             p_alt_expected = pchisq(ce, 1, lower.tail = FALSE),
             log10_p_gain_observed = log10(p0) - log10(p1), log10_p_gain_expected = log10(p0) - log10(pchisq(ce, 1, lower.tail = FALSE)),
             z_added_observed = if (nested) (z1 * sqrt(n1) - z0 * sqrt(n0)) / sqrt(n1 - n0) else NA,
             z_added_expected = if (nested) z0 * sqrt((n1 - n0) / n0) else NA, stringsAsFactors = FALSE)
}
A <- a[grepl("^A ", a$analysis), ]; B <- a[grepl("^B ", a$analysis), ]
P <- b[grepl("^phenotype-aware", b$rule), ]; D <- b[grepl("king-cutoff", b$rule), ]
if (nrow(A) != 1 || nrow(B) != 1 || nrow(P) != 1 || nrow(D) != 1) stop("rows not found in the input tables")
tab <- rbind(
  row("relatives retained", f4c, A$analysis, B$analysis, A$cases, A$controls, A$P_causal, B$cases, B$controls, B$P_causal, TRUE),
  row("relatedness rule", f4a, D$rule, P$rule, D$eur_cases, D$eur_controls, D$causal_p, P$eur_cases, P$eur_controls, P$causal_p, FALSE))
num <- sapply(tab, is.numeric); tabw <- tab; tabw[num] <- lapply(tabw[num], function(v) signif(v, 4))
write.table(tabw, file.path(OUT, "pdac_demo_08_expected_gain.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
R <- c("SECTION 8: OBSERVED AGAINST EXPECTED GAIN AT THE CAUSAL ABO VARIANT", "")
for (i in seq_len(nrow(tab))) with(tab[i, ], R <<- c(R,
  sprintf("%s (%s)", toupper(comparison), source),
  sprintf("  %-42s %d cases, %d controls, N_eff %.1f, P %.3g, chi-square %.2f", reference, cases_ref, controls_ref, n_eff_ref, p_ref, chisq_ref),
  sprintf("  %-42s %d cases, %d controls, N_eff %.1f, P %.3g, chi-square %.2f", alternative, cases_alt, controls_alt, n_eff_alt, p_alt, chisq_alt_observed),
  sprintf("  N_eff ratio %.3f: expected chi-square %.2f, expected P %.3g; gain in -log10 P expected %.2f, observed %.2f",
          n_eff_ratio, chisq_alt_expected, p_alt_expected, log10_p_gain_expected, log10_p_gain_observed),
  if (!is.na(z_added_observed)) sprintf("  z contributed by the added individuals: %.2f observed, %.2f expected at the same effect", z_added_observed, z_added_expected) else
    "  the two sets are not nested: no z of the added individuals", ""))
R <- c(R, "HOW TO READ THIS",
  "  The expected values hold the effect estimate and the allele frequency of",
  "  the reference analysis fixed and scale the chi-square by N_eff. What the",
  "  observed gain shows beyond the expected one is the draw of the individuals",
  "  that differ between the two analyses, in one phenotype.")
writeLines(R, file.path(OUT, "pdac_demo_08_expected_gain.txt")); cat(R, sep = "\n")
