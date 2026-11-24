#!/usr/bin/env Rscript
## Figure 2. Power: prediction against the replicates, and the detectable effect.
##  a  predicted power (formula of Section 3B) against realised power over 200
##     replicates, with 95% CI, per causal locus, European scan and
##     three-stratum meta-analysis.
##     Source: sections/08_replicates/results/pdac_demo_08_power.tsv
##  b  smallest OR detectable at 80% power, P < 5e-8, against MAF, for the
##     European analysis set and all strata combined. Curves use the Section 3B
##     functions (sourced from 01_power_curves.R), effective sizes read from
##     pdac_demo_07_power_summary.txt (labels) and computed from the case and
##     control counts of pdac_demo_05_strata_counts.tsv; points are the tabulated values of
##     pdac_demo_07_power_table.tsv. The three simulated causal variants are
##     drawn at their generative OR and at their minor allele frequency in the
##     European analysis set (pdac_demo_08_causal_structure.tsv, Section 8 step 07).
##     Bar: the OR at 80% power found by simulation in the European set
##     at one frequency, Wald test, with the risk on the minor or on the major
##     allele (pdac_demo_08_simulated_detectable_or.tsv, Section 8 step 10).

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fig_common.R"))

f_pow <- S("08_replicates", "pdac_demo_08_power.tsv")
f_sum <- S("03B_statistical_power", "pdac_demo_07_power_summary.txt")
f_tab <- S("03B_statistical_power", "pdac_demo_07_power_table.tsv")
f_scr <- "sections/03B_statistical_power/scripts/01_power_curves.R"
f_rd  <- "sections/08_replicates/README.md"
f_cnt <- S("05_meta_analysis", "pdac_demo_05_strata_counts.tsv")
f_frq <- S("08_replicates", "pdac_demo_08_causal_structure.tsv")
f_sim <- S("08_replicates", "pdac_demo_08_simulated_detectable_or.tsv")
sim <- rd(f_sim); sim <- sim[sim$test == "wald", ]
if (!setequal(sim$risk_allele, c("minor", "major"))) stop("unexpected rows in ", f_sim)

pw <- rd(f_pow); tb <- rd(f_tab)
sm <- readLines(need(f_sum))
grab <- function(pat) as.numeric(sub(".*: *([0-9.]+).*", "\\1", grep(pat, sm, value = TRUE)[1]))
N1_lab <- grab("^  effective sample size:"); N2_lab <- grab("combined effective sample size:")
alp <- as.numeric(sub(".*P < ([0-9.e-]+).*", "\\1", grep("^Threshold:", sm, value = TRUE)[1]))
## the Section 3B functions, verbatim from its script (lines defining crit, info, min_or)
src <- readLines(need(f_scr))
alpha <- alp
eval(parse(text = grep("^crit +<-", src, value = TRUE)))
eval(parse(text = src[grep("^info <- function", src)]))
eval(parse(text = src[grep("^neff <- function", src)]))
## effective sizes from the case and control counts, as Section 3B computes them
sc <- rd(f_cnt)
N1 <- with(sc[sc$stratum == "eur", ], neff(cases, controls)); N2 <- sum(neff(sc$cases, sc$controls))
if (abs(N1 - N1_lab) > 0.05 || abs(N2 - N2_lab) > 0.05) stop("effective sizes disagree with pdac_demo_07_power_summary.txt")
i0 <- grep("^min_or <- function", src); eval(parse(text = src[i0:(i0 + 4)]))
## check the curve function against the tabulated values
chk <- max(abs(round(sapply(tb$MAF, function(p) min_or(N1, p)), 2) - tb$min_OR_analysis_set),
           abs(round(sapply(tb$MAF, function(p) min_or(N2, p)), 2) - tb$min_OR_all_strata))
if (chk > 0) stop("Section 3B functions do not reproduce pdac_demo_07_power_table.tsv")

## locus labels and generative ORs from the Section 8 README table
rdm <- readLines(need(f_rd))
row_of <- function(v) strsplit(grep(v, rdm, value = TRUE, fixed = TRUE)[1], "\\|")[[1]]
loc <- data.frame(locus = c("G1", "G2", "G3"), id = c("9:133273682:A:T", "5:1286401:C:A", "13:73362027:C:T"))
loc$name <- sapply(loc$id, function(v) trimws(row_of(v)[2]))
loc$or <- sapply(loc$id, function(v) trimws(row_of(v)[4]))
fq <- rd(f_frq); fq <- fq[!duplicated(fq$locus), ]
loc$maf <- sapply(loc$locus, function(l) { f <- fq$risk_allele_freq[fq$locus == l]; min(f, 1 - f) })

draw <- function() {
  layout(matrix(1:2, 1), widths = c(1.15, 1))
  ## ---- a --------------------------------------------------------------------------
  par_base(mar = c(4.6, 3.6, 2.2, 0.6))
  an <- c("eur_4a", "meta"); xo <- c(eur_4a = -0.17, meta = 0.17)
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0, 103), xaxt = "n", yaxt = "n", xlab = "", ylab = expression("Power at P < 5" %*% 10^-8 ~ "(%)"), yaxs = "i")
  axis(2, at = seq(0, 100, 20))
  abline(h = seq(20, 80, 20), col = "grey92", lwd = 0.6)
  for (j in seq_len(nrow(loc))) for (a in an) {
    r <- pw[pw$locus == loc$locus[j] & pw$analysis == a, ]; x <- j + xo[[a]]
    col <- if (a == "eur_4a") ACCENT else ACCENT2
    segments(x, 100 * r$realised_l95, x, 100 * r$realised_u95, col = col, lwd = 1.2)
    points(x, 100 * r$realised_power, pch = 16, col = col, cex = 1.1)
    points(x, 100 * r$predicted_power_mean, pch = 4, col = INK, cex = 1.1, lwd = 1.1)
  }
  axis(1, at = 1:3, labels = FALSE)
  mtext(paste0(loc$name, "\nOR ", loc$or), side = 1, at = 1:3, line = 1.6, cex = 0.95, col = INK)
  legend(1.3, 80, bty = "n", cex = 0.95, xjust = 0, yjust = 1,
         legend = c("Predicted by the formula", "Realised, European scan (95% CI)", "Realised, meta-analysis (95% CI)"),
         pch = c(4, 16, 16), col = c(INK, ACCENT, ACCENT2), pt.lwd = c(1.1, 1, 1))
  panel_label("a")
  ## ---- b --------------------------------------------------------------------------
  par_base(mar = c(4.6, 3.6, 2.2, 0.8))
  p <- seq(0.01, 0.5, by = 0.0025)
  y1 <- sapply(p, function(q) min_or(N1, q)); y2 <- sapply(p, function(q) min_or(N2, q))
  plot(p, y1, type = "n", log = "y", ylim = c(1.3, 20), xlim = c(0, 0.52), xaxs = "i",
       xlab = "Minor allele frequency", ylab = "Smallest OR detectable, 80% power", yaxt = "n")
  yt <- c(1.5, 2, 3, 5, 10, 20); axis(2, at = yt, labels = yt)
  abline(h = yt, col = "grey92", lwd = 0.6)
  lines(p, y1, lwd = 1.4, col = ACCENT); lines(p, y2, lwd = 1.4, col = ACCENT2)
  points(tb$MAF, tb$min_OR_analysis_set, pch = 16, col = ACCENT, cex = 0.8)
  points(tb$MAF, tb$min_OR_all_strata, pch = 16, col = ACCENT2, cex = 0.8)
  ## OR at 80% power by simulation, European set, Wald test
  sm <- sim[match(c("minor", "major"), sim$risk_allele), ]
  segments(sm$maf[1], sm$or_at_80_power[1], sm$maf[2], sm$or_at_80_power[2], col = INK, lwd = 1.1)
  segments(sm$maf - 0.006, sm$or_at_80_power, sm$maf + 0.006, sm$or_at_80_power, col = INK, lwd = 1.1)
  ## the simulated causal variants: generative OR at the European minor allele frequency
  points(loc$maf, as.numeric(loc$or), pch = 23, bg = "white", col = INK, cex = 1.2, lwd = 1.1)
  text(loc$maf, as.numeric(loc$or), loc$name, pos = ifelse(loc$locus == "G1", 4, 2), offset = 0.5, cex = 0.9, col = INK)
  n1 <- formatC(N1_lab, format = "f", digits = 0, big.mark = ","); n2 <- formatC(N2_lab, format = "f", digits = 0, big.mark = ",")
  legend("topright", bty = "n", cex = 0.95, lwd = c(1.4, 1.4, 1.1, NA), pch = c(NA, NA, NA, 23),
         pt.bg = "white", pt.lwd = 1.1, pt.cex = 1.2, col = c(ACCENT, ACCENT2, INK, INK),
         legend = as.expression(c(bquote("Formula, European set, " * N[eff] ~ .(n1)), bquote("Formula, all strata, " * N[eff] ~ .(n2)),
                                  "Simulation, European set, Wald test", "Simulated causal variant")))
  panel_label("b")
}
save_fig("Fig2_power", 78, draw, c(f_pow, f_sum, f_tab, f_scr, f_rd, f_cnt, f_frq, f_sim))
