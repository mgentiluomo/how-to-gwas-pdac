#!/usr/bin/env Rscript
## Figure 3. Calibration of the scans.
##  a  lambda of the European logistic scan with 2, 5, 10 and 20 principal
##     components; the number chosen by 04A step 01b is filled.
##     Sources: sections/04A_association_binary/results/pdac_demo_04A_npc_scans.tsv,
##              pdac_demo_04A_npc_chosen.txt
##  b  permutation null of lambda in the Cox scan (100 permutations of the
##     outcome), with and without principal components, and the observed value.
##     Sources: sections/04B_association_survival/results/pdac_demo_04B_perm_lambda.tsv,
##              pdac_demo_04B_perm_summary.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fig_common.R"))

f_scan <- S("04A_association_binary", "pdac_demo_04A_npc_scans.tsv")
f_k    <- S("04A_association_binary", "pdac_demo_04A_npc_chosen.txt")
f_perm <- S("04B_association_survival", "pdac_demo_04B_perm_lambda.tsv")
f_psum <- S("04B_association_survival", "pdac_demo_04B_perm_summary.tsv")

sc <- rd(f_scan); k <- as.integer(read.table(need(f_k))[1, 2])
pl <- rd(f_perm); ps <- rd(f_psum)
if (!k %in% sc$n_pcs) stop("chosen number of components not among the scans")

draw <- function() {
  layout(matrix(1:2, 1), widths = c(0.8, 1.2))
  ## ---- a --------------------------------------------------------------------------
  par_base(mar = c(3.6, 3.8, 2.2, 0.8))
  yl <- c(1, max(sc$lambda) + 0.006)
  plot(sc$n_pcs, sc$lambda, type = "n", xlim = c(0, 22.5), ylim = yl, xaxs = "i",
       xlab = "Principal components in the model", ylab = expression(lambda), xaxt = "n")
  axis(1, at = sc$n_pcs)
  abline(h = pretty(yl), col = "grey92", lwd = 0.6)
  lines(sc$n_pcs, sc$lambda, col = MID, lwd = 1)
  points(sc$n_pcs, sc$lambda, pch = ifelse(sc$n_pcs == k, 16, 21), bg = "white", col = ACCENT, cex = 1.2)
  text(sc$n_pcs, sc$lambda, sprintf("%.4f", sc$lambda), pos = 3, offset = 0.5, cex = 0.9, col = INK)
  text(k, sc$lambda[sc$n_pcs == k], "chosen", pos = 1, offset = 0.6, cex = 0.9, col = ACCENT)
  panel_label("a")
  ## ---- b --------------------------------------------------------------------------
  par_base(mar = c(3.6, 3.6, 2.2, 0.8))
  mods <- c(no_pc = "SEX, AGE", primary = "SEX, AGE, within-set PC1-PC10")
  cols <- c(no_pc = MID, primary = ACCENT)
  perm <- lapply(names(mods), function(m) pl$lambda[pl$model == m & grepl("^perm", pl$fit)])
  names(perm) <- names(mods)
  if (any(lengths(perm) != ps$n_perm[match(names(mods), ps$model)])) stop("permutation count mismatch")
  rng <- range(unlist(perm), ps$lambda_observed); br <- seq(floor(rng[1] / 0.005) * 0.005, ceiling(rng[2] / 0.005) * 0.005 + 0.005, 0.005)
  h <- lapply(perm, hist, breaks = br, plot = FALSE)
  ymax <- max(sapply(h, function(x) max(x$counts)))
  plot(NA, xlim = range(br), ylim = c(0, ymax * 1.75), yaxs = "i", yaxt = "n",
       xlab = expression(lambda ~ "of the Cox scan"), ylab = "Permutations")
  axis(2, at = pretty(c(0, ymax)))
  for (m in names(mods)) {
    rect(h[[m]]$breaks[-length(br)], 0, h[[m]]$breaks[-1], h[[m]]$counts,
         col = adjustcolor(cols[[m]], 0.45), border = "white", lwd = 0.4)
    o <- ps$lambda_observed[ps$model == m]
    segments(o, 0, o, ymax * 1.08, col = cols[[m]], lwd = 1.6)
    text(o, ymax * 1.08, sprintf("observed %.4f\nempirical P %.2f", o, ps$empirical_P[ps$model == m]),
         pos = 3, offset = 0.2, cex = 0.85, col = cols[[m]], xpd = NA)
  }
  legend("topleft", inset = c(0.01, 0), bty = "n", cex = 0.95, fill = adjustcolor(cols, 0.45), border = NA,
         legend = sprintf("%s: mean %.4f, SD %.4f", mods, ps$perm_mean[match(names(mods), ps$model)],
                          ps$perm_sd[match(names(mods), ps$model)]))
  panel_label("b")
}
save_fig("Fig3_calibration", 72, draw, c(f_scan, f_k, f_perm, f_psum))
