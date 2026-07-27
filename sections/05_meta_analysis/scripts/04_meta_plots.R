################################################################################
# Section 5: Meta-analysis — Step 04: plots
#
# PURPOSE:
#   Four figures.
#
#     Manhattan          the meta-analysis result
#     Power comparison   what combining the strata actually bought
#     Forest plot        the lead variant, stratum by stratum
#     Heterogeneity      where the strata disagree, and why it matters
#
# INPUT:
#   - results/meta/pdac_demo_05_meta.tsv.gz
#
# OUTPUT:
#   - results/meta/pdac_demo_05_manhattan.png
#   - results/meta/pdac_demo_05_gain.png
#   - results/meta/pdac_demo_05_forest.png
#   - results/meta/pdac_demo_05_heterogeneity.png
################################################################################

out_dir <- "results/meta"
strata  <- c("eur", "afr", "eas")
lab     <- c(eur = "European", afr = "African", eas = "East Asian")
cols_s  <- c(eur = "#0072B2", afr = "#D55E00", eas = "#009E73")

d <- read.table(gzfile(file.path(out_dir, "pdac_demo_05_meta.tsv.gz")),
                header = TRUE, sep = "\t", check.names = FALSE,
                stringsAsFactors = FALSE)
d <- d[is.finite(d$P) & d$P > 0, ]
cat("Variants:", nrow(d), "\n")

# --- 1. Manhattan ------------------------------------------------------------
m <- d[order(d$CHROM, d$POS), ]
chr_len   <- tapply(m$POS, m$CHROM, max)
chr_start <- c(0, cumsum(as.numeric(chr_len))[-length(chr_len)])
names(chr_start) <- names(chr_len)
m$cum <- m$POS + chr_start[as.character(m$CHROM)]
axis_at <- tapply(m$cum, m$CHROM, function(x) (min(x) + max(x)) / 2)
cols <- ifelse(m$CHROM %% 2 == 0, "#4B6584", "#A5B1C2")
logp <- -log10(m$P)
keep <- logp > 1 | seq_len(nrow(m)) %% 4 == 0

png(file.path(out_dir, "pdac_demo_05_manhattan.png"),
    width = 2400, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 3, 1))
plot(m$cum[keep], logp[keep], col = cols[keep], pch = 19, cex = 0.35,
     xaxt = "n", xlab = "Chromosome",
     ylab = expression(-log[10](italic(P))),
     ylim = c(0, max(3, max(logp) * 1.12)),
     main = "Meta-analysis across three ancestry strata")
axis(1, at = axis_at, labels = names(axis_at), cex.axis = 0.7, tick = FALSE)
abline(h = -log10(5e-8), col = "#D55E00", lty = 2)
abline(h = -log10(1e-5), col = "#0072B2", lty = 3)
lead <- m[which.min(m$P), ]
text(lead$cum, -log10(lead$P), labels = lead$ID, pos = 2, cex = 0.65)
dev.off()

# --- 2. what the meta-analysis bought ----------------------------------------
# Every variant, European P against meta-analysis P. Points below the diagonal
# are variants where combining the strata strengthened the evidence.
png(file.path(out_dir, "pdac_demo_05_gain.png"),
    width = 1500, height = 1400, res = 200)
sel <- d$P_eur < 0.01 | d$P < 0.01
plot(-log10(d$P_eur[sel]), -log10(d$P[sel]),
     pch = 19, cex = 0.4, col = "#4B6584",
     xlab = expression(European~stratum~alone:~-log[10](italic(P))),
     ylab = expression(Meta-analysis:~-log[10](italic(P))),
     main = "What combining the strata bought")
abline(0, 1, col = "grey40", lty = 2)
abline(h = -log10(5e-8), col = "#D55E00", lty = 2)
abline(v = -log10(5e-8), col = "#D55E00", lty = 2)
l <- d[which.min(d$P), ]
points(-log10(l$P_eur), -log10(l$P), pch = 18, cex = 1.6, col = "#D55E00")
text(-log10(l$P_eur), -log10(l$P), labels = l$ID, pos = 2, cex = 0.65)
legend("topleft", bty = "n", cex = 0.8, lty = c(2, 2), col = c("grey40", "#D55E00"),
       legend = c("no gain", "genome-wide threshold"))
dev.off()

# --- 3. forest plot, lead variant --------------------------------------------
lead <- d[which.min(d$P), ]
est <- data.frame(
  label = c(lab[strata], "Meta-analysis (fixed)"),
  beta  = c(sapply(strata, function(g) lead[[paste0("BETA_", g)]]), lead$BETA),
  se    = c(sapply(strata, function(g) lead[[paste0("SE_",   g)]]), lead$SE),
  col   = c(cols_s[strata], "black"),
  is_meta = c(rep(FALSE, length(strata)), TRUE),
  stringsAsFactors = FALSE
)
est$or <- exp(est$beta)
est$lo <- exp(est$beta - 1.96 * est$se)
est$hi <- exp(est$beta + 1.96 * est$se)

png(file.path(out_dir, "pdac_demo_05_forest.png"),
    width = 1700, height = 1000, res = 200)
par(mar = c(4.5, 10, 3.5, 2))
yy <- rev(seq_len(nrow(est)))
xlim <- range(c(est$lo, est$hi, 1)) * c(0.9, 1.05)
plot(NA, xlim = xlim, ylim = c(0.5, nrow(est) + 0.5), log = "x",
     yaxt = "n", ylab = "", xlab = "Odds ratio (log scale)",
     main = sprintf("%s, by ancestry stratum", lead$ID))
abline(v = 1, col = "grey60", lty = 3)
abline(v = 2.40, col = "#009E73", lty = 2)     # the simulated truth
for (i in seq_len(nrow(est))) {
  y <- yy[i]
  segments(est$lo[i], y, est$hi[i], y, col = est$col[i], lwd = 2)
  points(est$or[i], y, pch = if (est$is_meta[i]) 18 else 19,
         cex = if (est$is_meta[i]) 1.8 else 1.2, col = est$col[i])
}
axis(2, at = yy, labels = est$label, las = 1, tick = FALSE, cex.axis = 0.85)
mtext(sprintf("I2 = %.0f%%,  P heterogeneity = %.2f", lead$I2, lead$P_het),
      side = 3, line = 0.2, cex = 0.8)
legend("bottomright", bty = "n", cex = 0.75, lty = 2, col = "#009E73",
       legend = "simulated true OR = 2.40")
dev.off()

# --- 4. heterogeneity in the top region --------------------------------------
# The causal variant and its correlated neighbours, coloured by I^2. The pattern
# is the point: the causal variant behaves the same in every population, while
# the variants that merely tag it do not, because linkage disequilibrium differs.
reg <- d[d$CHROM == lead$CHROM &
         d$POS > lead$POS - 250000 & d$POS < lead$POS + 250000, ]
reg <- reg[order(reg$POS), ]

png(file.path(out_dir, "pdac_demo_05_heterogeneity.png"),
    width = 1700, height = 1200, res = 200)
par(mar = c(4.5, 4.5, 3.5, 5))
i2 <- pmin(reg$I2, 100)
pal <- colorRampPalette(c("#0072B2", "#F0E442", "#D55E00"))(101)
plot(reg$POS / 1e6, -log10(reg$P), pch = 19, cex = 1.1,
     col = pal[floor(i2) + 1],
     xlab = sprintf("Position on chromosome %s (Mb)", lead$CHROM),
     ylab = expression(-log[10](italic(P))~meta),
     main = "Signal and heterogeneity around the lead variant")
points(lead$POS / 1e6, -log10(lead$P), pch = 5, cex = 2.2, lwd = 2)
abline(h = -log10(5e-8), col = "#D55E00", lty = 2)
legend("topright", bty = "n", cex = 0.75, pch = c(5, 19, 19, 19), pt.cex = c(1.4, 1, 1, 1),
       col = c("black", pal[1], pal[51], pal[101]),
       legend = c("lead variant", "I2 = 0%", "I2 = 50%", "I2 = 100%"))
dev.off()

cat("\nWritten four figures to", out_dir, "\n")
cat("\nLead variant across strata:\n")
for (i in seq_len(nrow(est))) {
  cat(sprintf("  %-22s OR %.3f (%.3f-%.3f)\n",
              est$label[i], est$or[i], est$lo[i], est$hi[i]))
}
cat(sprintf("\nHeterogeneity at the lead variant: I2 = %.1f%%, P = %.3f\n",
            lead$I2, lead$P_het))
cat(sprintf("European alone: P = %.3e     Meta-analysis: P = %.3e\n",
            lead$P_eur, lead$P))
