#!/usr/bin/env Rscript

################################################################################
#  make_manuscript_figures.R
#
#  Builds Figures 2, 3 and 4 of the manuscript at 300 dpi from the committed
#  results of the tutorial sections, so that every panel is reproducible from
#  the released data and scripts.
#
#    Figure 2  quality control trajectory
#    Figure 3  population structure
#    Figure 4  association results and what the study could detect
#
#  Figure 1 (workflow diagram) is drawn separately; Figure 5 (fine mapping)
#  awaits Section 6.
#
#  Usage:  Rscript make_manuscript_figures.R
################################################################################

out_dir <- "docs/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

blue <- "#0072B2"; orange <- "#D55E00"; green <- "#009E73"; grey <- "#4B6584"

# ==============================================================================
# Figure 2: quality control trajectory
# ==============================================================================
qc <- read.delim("sections/01B_genotyping_qc/results/pdac_demo_09_qc_counts.tsv",
                 stringsAsFactors = FALSE)
qc$label <- c("Raw", "Initial\nstats", "Sample\ncall rate", "Sex\ncheck",
              "Hetero-\nzygosity", "Variant\ncall rate", "Hardy-\nWeinberg",
              "Related-\nness", "MAF", "HWE\nwithin\nancestry")

png(file.path(out_dir, "Figure2_qc_trajectory.png"),
    width = 2400, height = 1500, res = 300)
par(mfrow = c(2, 1), mar = c(3.2, 5, 1.6, 5), mgp = c(3, 0.7, 0))

# samples
plot(seq_len(nrow(qc)), qc$samples, type = "o", pch = 19, col = blue, lwd = 2,
     xaxt = "n", xlab = "", ylab = "Samples retained",
     ylim = c(min(qc$samples) * 0.97, max(qc$samples) * 1.01))
axis(1, at = seq_len(nrow(qc)), labels = qc$label, cex.axis = 0.6, tick = FALSE)
lost <- qc$samples_lost; idx <- which(lost > 0)
text(idx, qc$samples[idx], labels = paste0("-", lost[idx]), pos = 1,
     cex = 0.65, col = orange)
mtext("a", side = 3, line = 0.3, adj = 0, font = 2, cex = 1.1)

# variants
plot(seq_len(nrow(qc)), qc$variants / 1000, type = "o", pch = 19, col = grey, lwd = 2,
     xaxt = "n", xlab = "", ylab = "Variants retained (thousands)",
     ylim = c(min(qc$variants) / 1000 * 0.97, max(qc$variants) / 1000 * 1.01))
axis(1, at = seq_len(nrow(qc)), labels = qc$label, cex.axis = 0.6, tick = FALSE)
lostv <- qc$variants_lost; idxv <- which(lostv > 0)
text(idxv, qc$variants[idxv] / 1000, labels = paste0("-", format(lostv[idxv], big.mark = ",")),
     pos = 1, cex = 0.65, col = orange)
mtext("b", side = 3, line = 0.3, adj = 0, font = 2, cex = 1.1)
dev.off()

# ==============================================================================
# Figure 3: population structure
# ==============================================================================
pcs <- read.table("results/pca/pdac_demo_02_pca_all.eigenvec", header = TRUE,
                  comment.char = "", check.names = FALSE)
names(pcs)[1] <- sub("^#", "", names(pcs)[1])
val <- scan("results/pca/pdac_demo_02_pca_all.eigenval", quiet = TRUE)
val_eur <- scan("results/pca/pdac_demo_02_pca_eur.eigenval", quiet = TRUE)
anc <- read.table("demo_data/sample_ancestry.tsv", col.names = c("IID", "g"),
                  stringsAsFactors = FALSE)
pcs <- merge(pcs, anc, by = "IID"); pcs$g <- toupper(pcs$g)
pal <- c(EUR = blue, AFR = orange, EAS = green)
cols <- pal[pcs$g]
pct <- 100 * val / sum(val); pct_e <- 100 * val_eur / sum(val_eur)

png(file.path(out_dir, "Figure3_population_structure.png"),
    width = 2400, height = 1700, res = 300)
layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
par(mar = c(4, 4.2, 2, 1), mgp = c(2.4, 0.7, 0))

plot(seq_along(pct)[1:10], pct[1:10], type = "b", pch = 19, col = blue,
     xlab = "Principal component", ylab = "Variance explained (%)")
mtext("a", 3, 0.3, adj = 0, font = 2, cex = 1.1)

plot(pcs$PC1, pcs$PC2, col = cols, pch = 19, cex = 0.55,
     xlab = sprintf("PC1 (%.1f%%)", pct[1]), ylab = sprintf("PC2 (%.1f%%)", pct[2]))
legend("topright", legend = names(pal), col = pal, pch = 19, bty = "n", cex = 0.8)
mtext("b", 3, 0.3, adj = 0, font = 2, cex = 1.1)

plot(pcs$PC3, pcs$PC4, col = cols, pch = 19, cex = 0.55,
     xlab = sprintf("PC3 (%.1f%%)", pct[3]), ylab = sprintf("PC4 (%.1f%%)", pct[4]))
mtext("c", 3, 0.3, adj = 0, font = 2, cex = 1.1)

plot(seq_along(pct_e)[1:10], pct_e[1:10], type = "b", pch = 19, col = blue,
     xlab = "Principal component", ylab = "Variance explained (%)")
mtext("d", 3, 0.3, adj = 0, font = 2, cex = 1.1)
dev.off()

# ==============================================================================
# Figure 4: association results and detectable effects
# ==============================================================================
f <- Sys.glob("results/assoc/pdac_demo_04A_gwas.*.glm.logistic.hybrid")
d <- read.table(f[1], header = TRUE, comment.char = "", check.names = FALSE)
names(d)[1] <- sub("^#", "", names(d)[1])
d <- d[d$TEST == "ADD" & is.finite(d$P) & d$P > 0, ]
d$CHR <- suppressWarnings(as.integer(d$CHROM))
d <- d[!is.na(d$CHR) & d$CHR %in% 1:22, ]
d <- d[order(d$CHR, d$POS), ]
len <- tapply(d$POS, d$CHR, max)
st  <- c(0, cumsum(as.numeric(len))[-length(len)]); names(st) <- names(len)
d$cum <- d$POS + st[as.character(d$CHR)]
axat <- tapply(d$cum, d$CHR, function(z) (min(z) + max(z)) / 2)
logp <- -log10(d$P)
keep <- logp > 1 | seq_len(nrow(d)) %% 5 == 0

png(file.path(out_dir, "Figure4_association_and_power.png"),
    width = 2400, height = 2100, res = 300)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))
par(mar = c(4, 4.4, 2, 1), mgp = c(2.5, 0.7, 0))

plot(d$cum[keep], logp[keep], col = ifelse(d$CHR[keep] %% 2 == 0, grey, "#A5B1C2"),
     pch = 19, cex = 0.3, xaxt = "n", xlab = "Chromosome",
     ylab = expression(-log[10](italic(P))), ylim = c(0, max(logp) * 1.15))
axis(1, at = axat, labels = names(axat), cex.axis = 0.55, tick = FALSE)
abline(h = -log10(5e-8), col = orange, lty = 2)
abline(h = -log10(1e-5), col = blue, lty = 3)
lead <- d[which.min(d$P), ]
text(lead$cum, -log10(lead$P), labels = "ABO", pos = 2, cex = 0.7, font = 3)
mtext("a", 3, 0.3, adj = 0, font = 2, cex = 1.1)

# QQ
pv <- sort(d$P); n <- length(pv)
obs <- -log10(pv); exp_ <- -log10(ppoints(n))
lam <- median(qchisq(1 - d$P, 1)) / qchisq(0.5, 1)
kq <- c(seq_len(min(n, 3000)), seq(3001, n, by = max(1, floor(n / 12000))))
kq <- unique(pmin(kq, n))
plot(exp_[kq], obs[kq], pch = 19, cex = 0.35, col = blue,
     xlab = expression(Expected~-log[10](italic(P))),
     ylab = expression(Observed~-log[10](italic(P))))
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", bty = "n", cex = 0.75, legend = sprintf("lambda = %.3f", lam))
mtext("b", 3, 0.3, adj = 0, font = 2, cex = 1.1)

# power curve
alpha <- 5e-8; ca <- qchisq(1 - alpha, 1)
lam80 <- uniroot(function(l) pchisq(ca, 1, ncp = l, lower.tail = FALSE) - 0.80,
                 c(1, 300))$root
minOR <- function(neff, p) exp(sqrt(lam80 / (neff * 0.25 * 2 * p * (1 - p))))
pgrid <- seq(0.02, 0.50, by = 0.005)
plot(pgrid, minOR(580.4, pgrid), type = "l", lwd = 2, col = orange, log = "y",
     ylim = c(1.1, 12), xlab = "Minor allele frequency",
     ylab = "Minimum odds ratio detectable\nat 80% power")
lines(pgrid, minOR(1189.7, pgrid), lwd = 2, col = green)
rect(0.02, 1.1, 0.50, 1.5, col = "#00000012", border = NA)
text(0.26, 1.28, "effect sizes of established\ncancer susceptibility loci", cex = 0.6)
points(0.085, 2.4, pch = 18, cex = 1.4)
text(0.085, 2.4, "simulated ABO effect", pos = 4, cex = 0.6, font = 3)
legend("topright", bty = "n", cex = 0.7, lwd = 2, col = c(orange, green),
       legend = c(expression(European~set~","~N[eff]==580),
                  expression(three~strata~combined~","~N[eff]==1190)))
mtext("c", 3, 0.3, adj = 0, font = 2, cex = 1.1)
dev.off()

cat("Written to", out_dir, ":\n")
cat("  Figure2_qc_trajectory.png\n  Figure3_population_structure.png\n  Figure4_association_and_power.png\n")
cat("\nlambda =", sprintf("%.4f", lam), " lambda for 80% power =", sprintf("%.2f", lam80), "\n")

# ==============================================================================
# Figure 5: fine mapping at ABO, and what ancestral diversity resolves
# ==============================================================================
fe <- read.delim("results/finemap/pdac_demo_06_pp_eur.tsv", stringsAsFactors = FALSE)
fm <- read.delim("results/finemap/pdac_demo_06_pp_meta.tsv", stringsAsFactors = FALSE)
TRUTH <- "9:133249045:A:G"
in_cs <- function(v) v == TRUE | v == "TRUE"

png(file.path(out_dir, "Figure5_fine_mapping.png"),
    width = 2400, height = 2000, res = 300)
layout(matrix(c(1, 1, 2, 3, 4, 4), 3, 2, byrow = TRUE), heights = c(1, 1, 0.95))
par(mar = c(4, 4.6, 2, 1), mgp = c(2.6, 0.7, 0))

# a: regional association
ce <- ifelse(in_cs(fe$in_cs), orange, grey)
plot(fe$POS / 1e6, -log10(fe$P), pch = 19, cex = 0.8, col = ce,
     xlab = "Position on chromosome 9 (Mb)", ylab = expression(-log[10](italic(P))))
abline(h = -log10(5e-8), col = orange, lty = 2)
tr <- fe[fe$ID == TRUTH, ]
points(tr$POS / 1e6, -log10(tr$P), pch = 5, cex = 2, lwd = 2)
legend("topleft", bty = "n", cex = 0.7, pch = c(19, 19, 5), col = c(orange, grey, "black"),
       legend = c("95% credible set", "outside the set", "simulated causal variant"))
mtext("a", 3, 0.3, adj = 0, font = 2, cex = 1.1)

# b: posterior in the European set
plot(fe$POS / 1e6, fe$PP, type = "h", lwd = 2, col = ce, ylim = c(0, 1),
     xlab = "Position (Mb)", ylab = "Posterior probability")
points(tr$POS / 1e6, tr$PP, pch = 5, cex = 1.8, lwd = 2)
mtext("b", 3, 0.3, adj = 0, font = 2, cex = 1.1)
mtext(sprintf("European set: %d variants", sum(in_cs(fe$in_cs))), 3, -1.2, cex = 0.65)

# c: posterior in the combined strata
cm <- ifelse(in_cs(fm$in_cs), green, grey)
plot(fm$POS / 1e6, fm$PP, type = "h", lwd = 2, col = cm, ylim = c(0, 1),
     xlab = "Position (Mb)", ylab = "Posterior probability")
trm <- fm[fm$ID == TRUTH, ]
points(trm$POS / 1e6, trm$PP, pch = 5, cex = 1.8, lwd = 2)
mtext("c", 3, 0.3, adj = 0, font = 2, cex = 1.1)
mtext(sprintf("three strata: %d variant", sum(in_cs(fm$in_cs))), 3, -1.2, cex = 0.65)

# d: how many variants are needed to reach 95%
par(mar = c(4.2, 4.6, 2, 1))
plot(cumsum(sort(fe$PP, decreasing = TRUE)), type = "s", lwd = 2, col = orange,
     xlim = c(0, 70), ylim = c(0, 1),
     xlab = "Variants, ranked by posterior probability",
     ylab = "Cumulative posterior probability")
lines(cumsum(sort(fm$PP, decreasing = TRUE)), type = "s", lwd = 2, col = green)
abline(h = 0.95, col = "grey50", lty = 3)
legend("bottomright", bty = "n", cex = 0.75, lwd = 2, col = c(orange, green),
       legend = c(sprintf("European set, %d variants to 95%%", sum(in_cs(fe$in_cs))),
                  sprintf("three strata, %d variant to 95%%", sum(in_cs(fm$in_cs)))))
mtext("d", 3, 0.3, adj = 0, font = 2, cex = 1.1)
dev.off()

cat("  Figure5_fine_mapping.png\n")
