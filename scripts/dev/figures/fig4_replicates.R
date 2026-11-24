#!/usr/bin/env Rscript
## Figure 4. Heterogeneity and fine mapping over 200 replicates (Section 8).
##  a  distribution of I-squared at the three causal variants in the
##     three-stratum meta-analysis, where the true effect is identical in every
##     stratum by construction; dashed line at 50%.
##     Sources: sections/08_replicates/results/pdac_demo_08_per_replicate.tsv
##              (column I2, meta rows), pdac_demo_08_i2.tsv (fraction above 50%)
##  b  ABF 95% credible sets: rank of the causal variant (stacked bars) and
##     coverage of the causal variant with 95% CI (points), European scan and
##     meta-analysis, all replicates.
##     Sources: pdac_demo_08_per_replicate.tsv (causal_rank),
##              pdac_demo_08_finemap.tsv (coverage_all and its CI)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fig_common.R"))

f_per <- S("08_replicates", "pdac_demo_08_per_replicate.tsv")
f_i2  <- S("08_replicates", "pdac_demo_08_i2.tsv")
f_fm  <- S("08_replicates", "pdac_demo_08_finemap.tsv")
f_rd  <- "sections/08_replicates/README.md"

per <- rd(f_per); i2 <- rd(f_i2); fm <- rd(f_fm)
rdm <- readLines(need(f_rd))
loc <- data.frame(locus = c("G1", "G2", "G3"), id = c("9:133273682:A:T", "5:1286401:C:A", "13:73362027:C:T"))
loc$name <- sapply(loc$id, function(v) trimws(strsplit(grep(v, rdm, value = TRUE, fixed = TRUE)[1], "\\|")[[1]][2]))
## check the plotted I2 against the summary table
for (l in loc$locus) {
  x <- per$I2[per$locus == l & per$analysis == "meta"]
  if (abs(mean(x > 50) - i2$frac_i2_gt_50[i2$locus == l]) > 1e-9) stop("I2 fractions disagree for ", l)
}

draw <- function() {
  layout(matrix(c(1, 4, 2, 4, 3, 4), 3, byrow = TRUE), widths = c(1, 1.15))
  par(cex = 1)
  ## ---- a: three small histograms ----------------------------------------------------
  br <- seq(0, 100, 5)
  for (j in seq_len(nrow(loc))) {
    l <- loc$locus[j]; x <- per$I2[per$locus == l & per$analysis == "meta"]
    h <- hist(x, breaks = br, plot = FALSE)
    par_base(mar = c(if (j == 3) 3.4 else 1.2, 3.6, if (j == 1) 2.2 else 0.6, 0.8)); par(cex = 1)
    plot(NA, xlim = c(0, 100), ylim = c(0, 200), xaxs = "i", yaxs = "i", xaxt = "n",
         xlab = if (j == 3) expression(I^2 ~ "(%), meta-analysis of three strata") else "", ylab = "")
    axis(1, at = seq(0, 100, 25), labels = j == 3)
    rect(br[-length(br)], 0, br[-1], h$counts, col = adjustcolor(ACCENT, 0.55), border = "white", lwd = 0.4)
    abline(v = 50, lty = 2, col = INK, lwd = 0.9)
    r <- i2[i2$locus == l, ]
    text(97, 185, adj = c(1, 1), cex = 0.95, col = INK,
         sprintf("%s\nI\u00b2 = 0 in %.0f%%; > 50%% in %.1f%%", loc$name[j], 100 * r$frac_i2_zero, 100 * r$frac_i2_gt_50))
    if (j == 2) mtext("Replicates", side = 2, line = 2.4, las = 0, cex = 1)
    if (j == 1) panel_label("a")
  }
  ## ---- b --------------------------------------------------------------------------
  par_base(mar = c(3.6, 3.6, 2.2, 0.8)); par(cex = 1)
  cats <- c("1", "2-3", "4-10", "> 10"); cuts <- c(0, 1, 3, 10, Inf)
  shades <- c(ACCENT, adjustcolor(ACCENT, 0.6), adjustcolor(ACCENT, 0.3), "grey88")
  an <- c(eur_4a = "Eur", meta = "Meta")
  grp <- expand.grid(analysis = names(an), locus = loc$locus, stringsAsFactors = FALSE)
  x <- seq_len(nrow(grp)) + rep(c(0, 0.6, 1.2), each = 2)
  plot(NA, xlim = c(0.4, max(x) + 0.6), ylim = c(0, 1.2), yaxs = "i", xaxt = "n", yaxt = "n", xlab = "",
       ylab = "Proportion of replicates")
  axis(2, at = seq(0, 1, 0.2))
  for (g in seq_len(nrow(grp))) {
    rk <- per$causal_rank[per$locus == grp$locus[g] & per$analysis == grp$analysis[g]]
    pr <- as.numeric(table(cut(rk, cuts, labels = cats))) / length(rk)
    rect(x[g] - 0.38, c(0, cumsum(pr)[-4]), x[g] + 0.38, cumsum(pr), col = shades, border = "white", lwd = 0.4)
    f <- fm[fm$locus == grp$locus[g] & fm$analysis == grp$analysis[g], ]
    segments(x[g], f$coverage_all_l95, x[g], f$coverage_all_u95, col = INK, lwd = 1)
    points(x[g], f$coverage_all, pch = 23, bg = "white", col = INK, cex = 1.1)
  }
  abline(h = 0.95, lty = 2, col = INK, lwd = 0.9)
  axis(1, at = x, labels = an[grp$analysis], tick = FALSE, line = -0.6)
  mtext(loc$name, side = 1, at = tapply(x, grp$locus, mean)[loc$locus], line = 1.4, cex = 1)
  legend("top", ncol = 3, bty = "n", cex = 0.9, inset = c(0, 0), x.intersp = 0.6, text.width = NULL,
         title = "Rank of the causal variant (ABF); coverage of the 95% credible set", title.adj = 0,
         legend = c(paste("rank", cats), "coverage, 95% CI"),
         pch = c(rep(22, 4), 23), pt.bg = c(shades, "white"), col = c(rep("white", 4), INK), pt.cex = c(rep(1.8, 4), 1.1))
  panel_label("b")
}
save_fig("Fig4_replicates", 92, draw, c(f_per, f_i2, f_fm, f_rd))
