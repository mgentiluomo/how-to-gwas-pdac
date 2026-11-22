#!/usr/bin/env Rscript

################################################################################
# Section 1B, step 10: what pooling does to two QC statistics
#
# PURPOSE
#   Quantify, from files the earlier steps already wrote, the two claims the
#   QC page makes about pooled statistics in a multi-ancestry sample:
#
#   1. Heterozygosity. Step 04 computes the inbreeding coefficient F within
#      each ancestry group. This step reports the same summary for the pooled
#      sample, to show what a pooled filter would have looked like.
#   2. Hardy-Weinberg. Step 06 tests the pooled controls and flags, without
#      removing, the variants below 1e-6; Section 2 then tests within each
#      group. This step asks how much of the pooled signal is predicted by
#      allele-frequency differences between groups alone (the Wahlund effect),
#      with no Hardy-Weinberg test involved, and counts the variants the
#      pooled test would have discarded for no reason.
#
# METHOD
#   Pooled heterozygosity: F from the step 01 --het run, restricted to the
#   samples that entered step 04. Step 01 ran on all samples, so the expected
#   homozygosity behind F uses allele frequencies from that larger set; the
#   summary rule (mean +/- 3 SD, sample SD) is the one step 04 applies.
#
#   Wahlund prediction: for each variant, the frequency of A1 in the controls
#   of each group (Section 2 .hardy files) is averaged with weights equal to
#   the ancestry composition of the pooled controls tested in step 06, and
#       F_ST = sum_g w_g (p_g - pbar)^2 / (pbar (1 - pbar)).
#   Pooling groups that differ in frequency creates a heterozygote deficit
#   equal to F_ST, so the expected pooled statistic is about N * F_ST^2 and
#   the observed one is N * F^2, with F = 1 - O(HET)/E(HET) from the pooled
#   .hardy file. Expected failures at 1e-6 are computed from the noncentral
#   chi-square with noncentrality N * F_ST^2.
#
# INPUT (all written by earlier steps of the canonical run)
#   results/qc/pdac_demo_01_qc.het
#   results/qc/pdac_demo_03_filt.fam
#   results/qc/pdac_demo_06_hwe.hardy, pdac_demo_06_controls.txt,
#   results/qc/pdac_demo_06_hwe_flagged.txt
#   results/qc/pdac_demo_08_filt.bim
#   results/pca/pdac_demo_02_hwe_{afr,eas,eur}.hardy, pdac_demo_02_hwe_exclude.txt
#   demo_data/sample_ancestry.tsv
#
# OUTPUT
#   results/qc/pdac_demo_10_het_pooled.tsv     metric/value
#   results/qc/pdac_demo_10_wahlund.tsv        metric/value
#   results/qc/pdac_demo_10_wahlund_bins.tsv   observed and expected failures by F_ST
#
# USAGE
#   Rscript 10_pooled_diagnostics.R      (from the project root)
#   Runs after Section 2 step 06, whose within-group files it reads.
################################################################################

suppressPackageStartupMessages(library(data.table))

QC  <- "results/qc"; PCA <- "results/pca"; DS <- "pdac_demo"
ANC <- "demo_data/sample_ancestry.tsv"
THR <- 1e-6; SD_MULT <- 3
GROUPS <- c("afr", "eas", "eur")

need <- c(file.path(QC, paste0(DS, c("_01_qc.het", "_03_filt.fam", "_06_hwe.hardy",
                                    "_06_controls.txt", "_06_hwe_flagged.txt", "_08_filt.bim"))),
          file.path(PCA, paste0(DS, "_02_hwe_", GROUPS, ".hardy")),
          file.path(PCA, paste0(DS, "_02_hwe_exclude.txt")), ANC)
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing input(s):\n  ", paste(miss, collapse = "\n  "))

kv_write <- function(x, f) fwrite(data.table(metric = names(x), value = unlist(x)), f, sep = "\t")
ids_of <- function(f) fread(f, header = FALSE)[[2]]

## ---------------------------------------------------------- heterozygosity ---
het <- fread(file.path(QC, paste0(DS, "_01_qc.het")))
setnames(het, sub("^#", "", names(het)))
in04 <- fread(file.path(QC, paste0(DS, "_03_filt.fam")), header = FALSE)[[2]]
f <- het[IID %in% in04, F]
m <- mean(f); s <- sd(f); lo <- m - SD_MULT * s; hi <- m + SD_MULT * s
het_out <- list(
  n = length(f), mean_F = sprintf("%.4f", m), sd_F = sprintf("%.4f", s),
  lower = sprintf("%.4f", lo), upper = sprintf("%.4f", hi),
  outliers = sum(f < lo | f > hi),
  frequency_basis_samples = nrow(het))
kv_write(het_out, file.path(QC, paste0(DS, "_10_het_pooled.tsv")))

## ------------------------------------------------------------------ HWE ---
rd <- function(f) { x <- fread(f); setnames(x, sub("^#", "", names(x))); x }
pooled <- rd(file.path(QC, paste0(DS, "_06_hwe.hardy")))
pooled[, N := HOM_A1_CT + HET_A1_CT + TWO_AX_CT]
pooled[, F_obs := 1 - `O(HET_A1)` / `E(HET_A1)`]

anc <- fread(ANC, header = FALSE, col.names = c("IID", "anc"))
ctrl <- ids_of(file.path(QC, paste0(DS, "_06_controls.txt")))
comp <- table(factor(anc[IID %in% ctrl, anc], levels = GROUPS))
w <- as.numeric(comp) / sum(comp); names(w) <- GROUPS

freq <- pooled[, .(ID, A1)]
for (g in GROUPS) {
  h <- rd(file.path(PCA, paste0(DS, "_02_hwe_", g, ".hardy")))
  h[, p := (2 * HOM_A1_CT + HET_A1_CT) / (2 * (HOM_A1_CT + HET_A1_CT + TWO_AX_CT))]
  h <- h[, .(ID, A1g = A1, AXg = AX, p)]
  freq <- merge(freq, h, by = "ID", all = FALSE)
  # align to the pooled file's A1
  freq[A1g != A1, p := 1 - p]
  setnames(freq, "p", paste0("p_", g)); freq[, c("A1g", "AXg") := NULL]
}
P <- as.matrix(freq[, paste0("p_", GROUPS), with = FALSE])
pbar <- drop(P %*% w)
fst <- drop(((P - pbar)^2) %*% w) / (pbar * (1 - pbar))
freq[, `:=`(pbar = pbar, F_ST = fst)]
d <- merge(pooled[, .(ID, N, F_obs, P)], freq[, .(ID, pbar, F_ST)], by = "ID")
d <- d[is.finite(F_obs) & is.finite(F_ST) & pbar > 0 & pbar < 1]

d[, obs_stat := N * F_obs^2]
d[, pred_stat := N * F_ST^2]
fit <- lm(obs_stat ~ pred_stat, d)
crit <- qchisq(THR, 1, lower.tail = FALSE)
d[, exp_fail := pchisq(crit, 1, ncp = pred_stat, lower.tail = FALSE)]
d[, fail := P < THR]

bins <- cut(d$F_ST, c(-Inf, 0.05, 0.10, Inf), right = FALSE,
            labels = c("F_ST < 0.05", "0.05 <= F_ST < 0.10", "F_ST >= 0.10"))
bt <- d[, .(variants = .N, observed_failures = sum(fail),
            expected_failures_wahlund = round(sum(exp_fail), 1),
            expected_failures_null = round(.N * THR, 2)), by = .(bin = bins)][order(bin)]
fwrite(bt, file.path(QC, paste0(DS, "_10_wahlund_bins.tsv")), sep = "\t")

## variants the pooled test discards for no reason: flagged in the pooled
## controls, surviving every other filter, and passing within every group
flagged <- fread(file.path(QC, paste0(DS, "_06_hwe_flagged.txt")), header = FALSE)[[1]]
final   <- fread(file.path(QC, paste0(DS, "_08_filt.bim")), header = FALSE)[[2]]
within  <- if (file.size(file.path(PCA, paste0(DS, "_02_hwe_exclude.txt"))) > 0)
  fread(file.path(PCA, paste0(DS, "_02_hwe_exclude.txt")), header = FALSE)[[1]] else character()
flag_final <- intersect(flagged, final)
needless <- setdiff(flag_final, within)

wa <- list(
  pooled_controls = sum(comp),
  pooled_controls_afr = comp[["afr"]], pooled_controls_eas = comp[["eas"]], pooled_controls_eur = comp[["eur"]],
  variants_compared = nrow(d),
  pearson_r = sprintf("%.3f", cor(d$obs_stat, d$pred_stat)),
  slope = sprintf("%.3f", coef(fit)[["pred_stat"]]),
  intercept = sprintf("%.3f", coef(fit)[["(Intercept)"]]),
  r_squared_pct = sprintf("%.1f", 100 * summary(fit)$r.squared),
  pearson_r_F_scale = sprintf("%.3f", cor(d$F_obs, d$F_ST)),
  failures_observed = sum(d$fail),
  failures_expected_wahlund = sprintf("%.1f", sum(d$exp_fail)),
  failures_share_fst_ge_0.10_pct = sprintf("%.1f", 100 * sum(d$fail & d$F_ST >= 0.10) / max(1, sum(d$fail))),
  flagged_pooled = length(flagged),
  flagged_pooled_in_final_set = length(flag_final),
  failing_within_ancestry = length(within),
  discarded_for_no_reason = length(needless))
kv_write(wa, file.path(QC, paste0(DS, "_10_wahlund.tsv")))

cat("Pooled heterozygosity (n ", het_out$n, "): mean F ", het_out$mean_F, ", SD ", het_out$sd_F,
    ", outliers at +/- ", SD_MULT, " SD: ", het_out$outliers, "\n", sep = "")
cat("Wahlund prediction across ", wa$variants_compared, " variants: r ", wa$pearson_r,
    ", slope ", wa$slope, ", R2 ", wa$r_squared_pct, "%\n", sep = "")
print(bt)
cat("Pooled test would discard ", wa$discarded_for_no_reason,
    " variants that pass every other filter and the within-ancestry test\n", sep = "")
cat("Written: ", paste(file.path(QC, paste0(DS, c("_10_het_pooled.tsv", "_10_wahlund.tsv", "_10_wahlund_bins.tsv"))), collapse = ", "), "\n")
