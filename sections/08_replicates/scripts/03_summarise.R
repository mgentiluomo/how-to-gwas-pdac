#!/usr/bin/env Rscript

################################################################################
# Section 8, step 03: summarise the replicates
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_power.tsv          realised against predicted power, per locus and analysis
#   pdac_demo_08_winners_curse.tsv  estimated / generative OR, all replicates and given P < 5e-8
#   pdac_demo_08_finemap.tsv        ABF 95% credible-set coverage, size and rank of the causal variant
#   pdac_demo_08_i2.tsv             I-squared at the three causal variants (meta-analysis)
#   pdac_demo_08_lambda_summary.tsv distribution of lambda (LD-pruned set)
#   pdac_demo_08_summary.txt        plain-text report
#   pdac_demo_08_i2.png, pdac_demo_08_winners_curse.png
#
# Intervals: Wilson 95% intervals for proportions; 2.5th and 97.5th
# percentiles across replicates for distributions; for the winner's curse on
# the log OR scale, percentile bootstrap over replicates (10,000 resamples,
# seed 20261124). Each resample draws replicates with replacement and
# recomputes the mean over all of them and the mean among the significant
# ones, so the conditional interval carries the uncertainty in how many are
# significant; resamples without a significant replicate are counted and
# left out of the conditional interval.
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(normalizePath(file.path(script_dir, "..", "..", "..")))
OUT <- "sections/08_replicates/results"
d <- read.delim(file.path(OUT, "pdac_demo_08_per_replicate.tsv"), stringsAsFactors = FALSE)
lam <- read.delim(file.path(OUT, "pdac_demo_08_lambda.tsv"))
rz <- read.delim(file.path(OUT, "pdac_demo_08_realised.tsv"))
OR_GEN <- c(G1 = 2.60, G2 = 1.50, G3 = 1.97)
LABEL <- c(G1 = "ABO 9:133273682:A:T", G2 = "TERT/CLPTM1L 5:1286401:C:A", G3 = "13q22.1 13:73362027:C:T")
d <- d[d$present, ]
N <- length(unique(d$replicate))
wilson <- function(x, n) { if (n == 0) return(c(NA, NA)); p <- x / n; z <- 1.959964
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
grp <- split(d, list(d$locus, d$analysis), drop = TRUE)

power <- do.call(rbind, lapply(grp, function(x) { k <- sum(x$significant); ci <- wilson(k, nrow(x))
  data.frame(locus = x$locus[1], analysis = x$analysis[1], replicates = nrow(x),
             predicted_power_mean = mean(x$predicted_power), realised_power = k / nrow(x),
             realised_l95 = ci[1], realised_u95 = ci[2]) }))
wc <- do.call(rbind, lapply(grp, function(x) { g <- OR_GEN[[x$locus[1]]]; s <- x[x$significant, ]
  data.frame(locus = x$locus[1], analysis = x$analysis[1], generative_or = g,
             mean_ratio_all = mean(x$or_hat / g), median_ratio_all = median(x$or_hat / g),
             n_significant = nrow(s),
             mean_ratio_significant = if (nrow(s)) mean(s$or_hat / g) else NA,
             median_ratio_significant = if (nrow(s)) median(s$or_hat / g) else NA,
             mean_logratio_significant = if (nrow(s)) mean(s$log_or_hat / log(g)) else NA,
             mean_ratio_not_significant = if (nrow(s) < nrow(x)) mean(x$or_hat[!x$significant] / g) else NA) }))
## winner's curse on the log OR scale: estimated / generative log OR, with bootstrap intervals
set.seed(20261124L); B_BOOT <- 10000L
wcl <- do.call(rbind, lapply(grp, function(x) { g <- log(OR_GEN[[x$locus[1]]]); lr <- x$log_or_hat / g; sg <- x$significant; n <- nrow(x)
  idx <- matrix(sample.int(n, n * B_BOOT, replace = TRUE), n, B_BOOT)
  L <- matrix(lr[idx], n, B_BOOT); S <- matrix(sg[idx], n, B_BOOT); ns <- colSums(S)
  ball <- colMeans(L); bsig <- (colSums(L * S) / ns)[ns > 0]; bor <- (colSums(exp(L * g) / exp(g) * S) / ns)[ns > 0]
  qq <- function(v, p) if (length(v)) unname(quantile(v, p)) else NA
  data.frame(mean_logratio_all = mean(lr), mean_logratio_all_l95 = qq(ball, 0.025), mean_logratio_all_u95 = qq(ball, 0.975),
             mean_logratio_significant_l95 = if (any(sg)) qq(bsig, 0.025) else NA, mean_logratio_significant_u95 = if (any(sg)) qq(bsig, 0.975) else NA,
             mean_ratio_significant_l95 = if (any(sg)) qq(bor, 0.025) else NA, mean_ratio_significant_u95 = if (any(sg)) qq(bor, 0.975) else NA,
             min_logratio_significant = if (any(sg)) min(lr[sg]) else NA, min_ratio_significant = if (any(sg)) min(x$or_hat[sg]) / exp(g) else NA,
             mean_logratio_not_significant = if (any(!sg)) mean(lr[!sg]) else NA,
             bootstrap_resamples = B_BOOT, resamples_without_significant = sum(ns == 0)) }))
stopifnot(identical(rownames(wcl), rownames(wc))); wc <- cbind(wc, wcl)
fm <- do.call(rbind, lapply(grp, function(x) { s <- x[x$region_significant, ]
  ci <- wilson(sum(x$causal_in_cs), nrow(x)); cis <- wilson(sum(s$causal_in_cs), nrow(s))
  data.frame(locus = x$locus[1], analysis = x$analysis[1],
             coverage_all = mean(x$causal_in_cs), coverage_all_l95 = ci[1], coverage_all_u95 = ci[2],
             n_region_significant = nrow(s),
             coverage_region_significant = if (nrow(s)) mean(s$causal_in_cs) else NA,
             coverage_sig_l95 = cis[1], coverage_sig_u95 = cis[2],
             median_cs_size = median(x$cs_size), median_cs_size_significant = if (nrow(s)) median(s$cs_size) else NA,
             causal_rank1 = mean(x$causal_rank == 1, na.rm = TRUE),
             median_causal_rank = median(x$causal_rank, na.rm = TRUE),
             causal_rank1_significant = if (nrow(s)) mean(s$causal_rank == 1, na.rm = TRUE) else NA) }))
mi <- d[d$analysis == "meta", ]
i2 <- do.call(rbind, lapply(split(mi, mi$locus), function(x) data.frame(
  locus = x$locus[1], replicates = nrow(x), i2_mean = mean(x$I2), i2_median = median(x$I2),
  i2_q025 = q(x$I2, 0.025), i2_q975 = q(x$I2, 0.975), frac_i2_zero = mean(x$I2 == 0),
  frac_i2_gt_50 = mean(x$I2 > 50), frac_phet_lt_0.05 = mean(x$P_het < 0.05),
  frac_random_weaker = mean(x$p_random > x$p))))
lsum <- do.call(rbind, lapply(setdiff(names(lam), "replicate"), function(cn) data.frame(
  lambda = cn, mean = mean(lam[[cn]]), sd = sd(lam[[cn]]), q025 = q(lam[[cn]], 0.025),
  median = median(lam[[cn]]), q975 = q(lam[[cn]], 0.975), min = min(lam[[cn]]), max = max(lam[[cn]]))))
for (nm in c("power", "wc", "fm", "i2", "lsum")) {
  x <- get(nm); num <- sapply(x, is.numeric); x[num] <- lapply(x[num], function(v) signif(v, 4)); assign(nm, x) }
wr <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
wr(power, "pdac_demo_08_power.tsv"); wr(wc, "pdac_demo_08_winners_curse.tsv"); wr(fm, "pdac_demo_08_finemap.tsv")
wr(i2, "pdac_demo_08_i2.tsv"); wr(lsum, "pdac_demo_08_lambda_summary.tsv")

png(file.path(OUT, "pdac_demo_08_i2.png"), width = 1800, height = 700, res = 200); par(mfrow = c(1, 3))
for (L in names(OR_GEN)) { x <- mi$I2[mi$locus == L]
  hist(x, breaks = seq(0, 100, 5), col = "#56B4E9", border = "white", main = LABEL[[L]], cex.main = 0.8,
       xlab = "I-squared (%)", ylab = "Replicates"); abline(v = 50, lty = 2) }
invisible(dev.off())
png(file.path(OUT, "pdac_demo_08_winners_curse.png"), width = 1800, height = 700, res = 200); par(mfrow = c(1, 3))
for (L in names(OR_GEN)) { x <- d[d$locus == L & d$analysis == "eur_4a", ]
  br <- pretty(range(x$or_hat), 30)
  hist(x$or_hat[!x$significant], breaks = br, col = adjustcolor("grey60", 0.7), border = "white",
       main = paste(LABEL[[L]], "(European scan)"), cex.main = 0.8, xlab = "Estimated OR", ylab = "Replicates",
       ylim = c(0, max(hist(x$or_hat, breaks = br, plot = FALSE)$counts)))
  if (any(x$significant)) hist(x$or_hat[x$significant], breaks = br, col = adjustcolor("#D55E00", 0.7), border = "white", add = TRUE)
  abline(v = OR_GEN[[L]], lwd = 2) }
invisible(dev.off())

R <- c("SECTION 8: REPLICATES OF THE SIMULATED PHENOTYPE",
       sprintf("%d replicates; generative ORs: ABO 2.60, TERT/CLPTM1L 1.50, 13q22.1 1.97; P < 5e-8", N),
       sprintf("realised ORs on all 1,461 (mean over replicates): ABO %.3f, TERT %.3f, 13q22.1 %.3f",
               mean(rz$realised_or_G1), mean(rz$realised_or_G2), mean(rz$realised_or_G3)), "",
       "POWER, realised (95% CI) against predicted")
for (i in seq_len(nrow(power))) with(power[i, ], R <<- c(R, sprintf("  %-3s %-7s predicted %5.1f%%  realised %5.1f%% (%.1f to %.1f)",
  locus, analysis, 100 * predicted_power_mean, 100 * realised_power, 100 * realised_l95, 100 * realised_u95)))
R <- c(R, "", "WINNER'S CURSE: estimated / generative OR (mean)")
for (i in seq_len(nrow(wc))) with(wc[i, ], R <<- c(R, sprintf("  %-3s %-7s all %.3f   given P < 5e-8 %s (n %d)   not significant %s",
  locus, analysis, mean_ratio_all, if (is.na(mean_ratio_significant)) "NA" else sprintf("%.3f", mean_ratio_significant),
  n_significant, if (is.na(mean_ratio_not_significant)) "NA" else sprintf("%.3f", mean_ratio_not_significant))))
R <- c(R, "", "WINNER'S CURSE ON THE LOG OR SCALE: estimated / generative log OR (mean, bootstrap 95% interval)")
for (i in seq_len(nrow(wc))) with(wc[i, ], R <<- c(R, sprintf("  %-3s %-7s all %.3f (%.3f to %.3f)   given P < 5e-8 %s (n %d)   smallest significant %s",
  locus, analysis, mean_logratio_all, mean_logratio_all_l95, mean_logratio_all_u95,
  if (is.na(mean_logratio_significant)) "NA" else sprintf("%.3f (%.3f to %.3f)", mean_logratio_significant, mean_logratio_significant_l95, mean_logratio_significant_u95),
  n_significant, if (is.na(min_logratio_significant)) "NA" else sprintf("%.3f", min_logratio_significant))))
R <- c(R, "", "ABF 95% CREDIBLE SETS: coverage of the causal variant")
for (i in seq_len(nrow(fm))) with(fm[i, ], R <<- c(R, sprintf("  %-3s %-7s all %.3f   region significant %s (n %d)   median size %s   causal ranked first %.2f",
  locus, analysis, coverage_all, if (is.na(coverage_region_significant)) "NA" else sprintf("%.3f", coverage_region_significant),
  n_region_significant, median_cs_size, causal_rank1)))
R <- c(R, "", "I-SQUARED AT THE CAUSAL VARIANTS (identical effect in the three strata by construction)")
for (i in seq_len(nrow(i2))) with(i2[i, ], R <<- c(R, sprintf("  %-3s median %.1f%%  95%% range %.1f to %.1f%%  I2 = 0 in %.2f  I2 > 50%% in %.2f  P_het < 0.05 in %.2f",
  locus, i2_median, i2_q025, i2_q975, frac_i2_zero, frac_i2_gt_50, frac_phet_lt_0.05)))
R <- c(R, "", "LAMBDA, LD-pruned set")
for (i in seq_len(nrow(lsum))) with(lsum[i, ], R <<- c(R, sprintf("  %-20s mean %.4f  sd %.4f  2.5-97.5%% %.4f to %.4f", lambda, mean, sd, q025, q975)))
writeLines(R, file.path(OUT, "pdac_demo_08_summary.txt")); cat(R, sep = "\n")
