#!/usr/bin/env Rscript

################################################################################
# Statistical power — Step 01: what this study can and cannot detect
#
# PURPOSE:
#   Compute, for the analysis set actually in hand, the smallest effect that
#   could be detected. Done before the association test, this answers the only
#   question that matters at this stage: is the planned analysis capable of
#   finding what it is looking for?
#
# THE ARITHMETIC, IN THREE LINES:
#   Power depends on how much information the design carries about the genetic
#   coefficient. Under an additive model that information is a product:
#
#       information  =  N  x  pi(1 - pi)  x  Var(G),      Var(G) = 2p(1 - p)
#
#   the sample size, the outcome variability, and the genotype variability.
#   Multiplying by the squared effect gives the non-centrality parameter of a
#   one degree of freedom chi-squared statistic, and power is the probability
#   that this exceeds the critical value set by the significance threshold.
#
#   Each term earns its place. Information rises with allele frequency, so rare
#   variants are hard; pi(1 - pi) is maximised at pi = 0.5, so a balanced design
#   is the most efficient; and only the sample size can be bought.
#
# EFFECTIVE SAMPLE SIZE:
#   An unbalanced study behaves like a smaller balanced one. With Ncase cases
#   and Nctrl controls the equivalent balanced size is
#
#       N_eff = 4 x Ncase x Nctrl / (Ncase + Nctrl)
#
#   Report this, not the headline total. Beyond about four controls per case,
#   additional controls buy very little.
#
# INPUT:
#   - results/pca/pdac_demo_02_analysis_set.tsv       (Section 2)
#   - results/meta/pdac_demo_05_strata_counts.tsv     (Section 5, optional)
#
# OUTPUT:
#   - results/power/pdac_demo_07_power_table.tsv
#   - results/power/pdac_demo_07_power_curves.png
#   - results/power/pdac_demo_07_power_summary.txt
################################################################################

out_dir <- "results/power"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha  <- 5e-8                       # genome-wide significance
crit   <- qchisq(1 - alpha, df = 1)  # critical value

neff <- function(ca, co) 4 * ca * co / (ca + co)

# information for one variant, additive model, effective (balanced) sample size
info <- function(n_eff, maf) n_eff * 0.25 * 2 * maf * (1 - maf)

power_of <- function(n_eff, maf, or)
  pchisq(crit, df = 1, ncp = info(n_eff, maf) * log(or)^2, lower.tail = FALSE)

# the effect just detectable at a given power
min_or <- function(n_eff, maf, target = 0.80) {
  lam <- uniroot(function(l) pchisq(crit, 1, ncp = l, lower.tail = FALSE) - target,
                 c(1e-6, 500))$root
  exp(sqrt(lam / info(n_eff, maf)))
}

# --- the analysis set actually in hand ---------------------------------------
as_file <- "results/pca/pdac_demo_02_analysis_set.tsv"
if (file.exists(as_file)) {
  a  <- read.delim(as_file, stringsAsFactors = FALSE)
  ca <- as.numeric(a$value[a$quantity == "analysis_set_cases"])
  co <- as.numeric(a$value[a$quantity == "analysis_set_controls"])
} else { ca <- 229; co <- 396 }
N1 <- neff(ca, co)

# --- all strata combined, if the meta-analysis has been run -------------------
st_file <- "results/meta/pdac_demo_05_strata_counts.tsv"
have_meta <- file.exists(st_file)
if (have_meta) {
  s  <- read.delim(st_file, stringsAsFactors = FALSE)
  N2 <- sum(neff(s$cases, s$controls))
} else N2 <- NA

con <- file(file.path(out_dir, "pdac_demo_07_power_summary.txt"), open = "wt")
w <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = con) }

w("Statistical power of the analysis set")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("Threshold: P < ", alpha, "   target power: 80%")
w("")
w("ANALYSIS SET")
w("  cases: ", ca, "   controls: ", co, "   total: ", ca + co)
w("  effective sample size: ", sprintf("%.1f", N1))
w("  Note that ", ca + co, " people behave like ", sprintf("%.0f", N1),
  " in a balanced design. The imbalance is not free.")
if (have_meta) {
  w("")
  w("ALL STRATA COMBINED")
  for (i in seq_len(nrow(s)))
    w("  ", s$stratum[i], ": ", s$cases[i], " cases, ", s$controls[i],
      " controls, effective ", sprintf("%.1f", neff(s$cases[i], s$controls[i])))
  w("  combined effective sample size: ", sprintf("%.1f", N2))
}
w("")

# --- the table the paper should carry ----------------------------------------
mafs <- c(0.50, 0.30, 0.20, 0.10, 0.05, 0.02)
tab <- data.frame(MAF = mafs,
                  min_OR_analysis_set = sapply(mafs, function(p) min_or(N1, p)))
if (have_meta) tab$min_OR_all_strata <- sapply(mafs, function(p) min_or(N2, p))
tab[-1] <- round(tab[-1], 2)
write.table(tab, file.path(out_dir, "pdac_demo_07_power_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

w("SMALLEST ODDS RATIO DETECTABLE AT 80% POWER")
w(paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n"))
w("")

# --- what this means for effects of a realistic size -------------------------
w("POWER AT EFFECT SIZES THAT ACTUALLY OCCUR")
for (or in c(1.2, 1.3, 1.5, 2.0)) {
  w(sprintf("  OR %.1f at MAF 0.30: %5.1f%%   at MAF 0.10: %5.1f%%",
            or, 100 * power_of(N1, 0.30, or), 100 * power_of(N1, 0.10, or)))
}
w("  Established cancer susceptibility loci sit between 1.1 and 1.5. Read the")
w("  line for 1.3 and then read it again: that is the answer to whether this")
w("  study could have found a typical locus.")
w("")

# --- the simulated causal variants -------------------------------------------
# Read from the released truth table and the observed allele frequencies, not
# from constants. A constant here goes stale the moment the phenotype is
# regenerated, which is exactly how a previous version of this dataset came to
# document an effect it did not contain.
truth <- read.delim("demo_data/truth.tsv", stringsAsFactors = FALSE)
truth <- truth[grepl(":", truth$variant), ]

af <- read.delim("results/qc/pdac_demo_01_qc.afreq", check.names = FALSE,
                 comment.char = "", stringsAsFactors = FALSE)
names(af) <- sub("^#", "", names(af))
truth$MAF <- sapply(truth$variant, function(v) {
  f <- af$ALT_FREQS[match(v, af$ID)]
  if (is.na(f)) NA_real_ else min(f, 1 - f)
})

w("THE SIMULATED CAUSAL VARIANTS")
w("  Power computed at the generative odds ratio and the observed frequency.")
for (i in seq_len(nrow(truth))) {
  v <- truth$variant[i]; or <- truth$generative_OR[i]; maf <- truth$MAF[i]
  if (is.na(maf)) { w("  ", v, ": not present in the QC-passed data"); next }
  w(sprintf("  %s  (%s)", v, truth$locus[i]))
  w(sprintf("    generative OR %.2f at frequency %.3f", or, maf))
  w(sprintf("    power in the analysis set: %5.1f%%", 100 * power_of(N1, maf, or)))
  if (have_meta)
    w(sprintf("    power across all strata:   %5.1f%%", 100 * power_of(N2, maf, or)))
}
w("")
w("  Read the two together. Whether a scan finds a locus is decided here,")
w("  before any analysis is run.")
w("")
w("  Where a variant is underpowered in the analysis set, combining the strata")
w("  raises its power but does not necessarily rescue it. Compare the two")
w("  percentages above for each variant, then look at what Section 5 actually")
w("  recovers. A large relative gain in power can still leave a real effect far")
w("  short of genome-wide significance, and that gap is the honest case for")
w("  consortia: meta-analysis is the only lever available, not a rescue.")
w("")
w("  NOTE: selection bias still applies to whatever crosses the threshold. The")
w("  size of that bias is reported in Section 4A by comparing the recovered")
w("  odds ratio with the generative one in truth.tsv.")
close(con)

# --- figure ------------------------------------------------------------------
p <- seq(0.01, 0.50, by = 0.002)
png(file.path(out_dir, "pdac_demo_07_power_curves.png"),
    width = 1700, height = 1300, res = 200)
par(mar = c(4.4, 4.6, 2.5, 1), mgp = c(2.6, 0.7, 0))
plot(p, sapply(p, function(q) min_or(N1, q)), type = "l", lwd = 2, col = "#D55E00",
     log = "y", ylim = c(1.1, 20), xlab = "Minor allele frequency",
     ylab = "Smallest odds ratio detectable at 80% power",
     main = "What this study could have found")
if (have_meta) lines(p, sapply(p, function(q) min_or(N2, q)), lwd = 2, col = "#009E73")
rect(0.01, 1.1, 0.50, 1.5, col = "#00000010", border = NA)
text(0.28, 1.27, "effect sizes of established cancer susceptibility loci", cex = 0.7)
ok <- !is.na(truth$MAF)
points(truth$MAF[ok], truth$generative_OR[ok], pch = 18, cex = 1.6)
text(truth$MAF[ok], truth$generative_OR[ok], truth$locus[ok], pos = 4, cex = 0.7)
legend("topright", bty = "n", lwd = 2, cex = 0.8,
       col = if (have_meta) c("#D55E00", "#009E73") else "#D55E00",
       legend = if (have_meta)
         c(sprintf("analysis set, N_eff = %.0f", N1),
           sprintf("all strata combined, N_eff = %.0f", N2))
       else sprintf("analysis set, N_eff = %.0f", N1))
dev.off()

cat("\nWritten to ", out_dir, "\n", sep = "")
cat("\n=== NEXT STEP ===\n\n")
cat("  Rscript scripts/03B_statistical_power/02_power_events.R\n\n")
