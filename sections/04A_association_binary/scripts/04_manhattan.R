################################################################################
# Section 4A: Association testing — Step 04: Manhattan plot
#
# PURPOSE:
#   Display -log10(P) against genomic position, with the two thresholds that
#   structure interpretation.
#
#     5e-8   genome-wide significance. Approximately a Bonferroni correction for
#            one million independent tests in European-ancestry data. It is
#            somewhat conservative for a single-ancestry analysis with fewer
#            independent tests, and anticonservative for trans-ancestry work
#            where LD blocks are shorter.
#
#     1e-5   the conventional suggestive threshold. Variants here are not
#            findings and must never be presented as such. In a rare cancer,
#            however, discarding them entirely also discards the most plausible
#            candidates for the next consortium round. Report them transparently,
#            flagged as requiring replication, next to the power curve that
#            explains why the study could not resolve them.
#
#   Expect this plot to look sparse. On the demonstration data one locus rises
#   above the genome-wide line and the rest of the genome is flat. That is the
#   correct result, not a failed analysis: the phenotype was simulated from a
#   single causal locus, and at this effective sample size nothing weaker could
#   have been detected anyway.
#
# INPUT:
#   - results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid
#
# OUTPUT:
#   - results/assoc/pdac_demo_04A_manhattan.png
################################################################################

out_dir <- "results/assoc"

f <- Sys.glob(file.path(out_dir, "pdac_demo_04A_gwas.*.glm.logistic.hybrid"))
if (!length(f)) f <- Sys.glob(file.path(out_dir, "pdac_demo_04A_gwas.*.glm.logistic"))
if (!length(f)) stop("No association output found in ", out_dir)

d <- read.table(f[1], header = TRUE, comment.char = "", check.names = FALSE)
names(d)[1] <- sub("^#", "", names(d)[1])
if ("TEST" %in% names(d)) d <- d[d$TEST == "ADD", ]
d <- d[is.finite(d$P) & d$P > 0, ]

d$CHR <- suppressWarnings(as.integer(sub("^chr", "", d$CHROM)))
d <- d[!is.na(d$CHR) & d$CHR >= 1 & d$CHR <= 22, ]
d <- d[order(d$CHR, d$POS), ]

# Cumulative position, so that chromosomes sit side by side on one axis.
chr_len   <- tapply(d$POS, d$CHR, max)
chr_start <- c(0, cumsum(as.numeric(chr_len))[-length(chr_len)])
names(chr_start) <- names(chr_len)
d$cum <- d$POS + chr_start[as.character(d$CHR)]

axis_at <- tapply(d$cum, d$CHR, function(x) (min(x) + max(x)) / 2)

# Alternate two shades so that chromosome boundaries are visible.
cols <- ifelse(d$CHR %% 2 == 0, "#4B6584", "#A5B1C2")

# Thin the uninformative bulk: points below -log10(P) = 1 are drawn sparsely.
logp   <- -log10(d$P)
is_top <- logp > 1
keep   <- is_top | (seq_len(nrow(d)) %% 4 == 0)

png(file.path(out_dir, "pdac_demo_04A_manhattan.png"),
    width = 2400, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 3, 1))
plot(d$cum[keep], logp[keep],
     col = cols[keep], pch = 19, cex = 0.35,
     xaxt = "n", xlab = "Chromosome",
     ylab = expression(-log[10](italic(P))),
     ylim = c(0, max(3, max(logp) * 1.15)),
     main = "Association results, European analysis set")
axis(1, at = axis_at, labels = names(axis_at), cex.axis = 0.7, tick = FALSE)
abline(h = -log10(5e-8), col = "#D55E00", lty = 2)
abline(h = -log10(1e-5), col = "#0072B2", lty = 3)
legend("topright", bty = "n", cex = 0.8,
       legend = c("genome-wide, 5e-8", "suggestive, 1e-5"),
       col = c("#D55E00", "#0072B2"), lty = c(2, 3))

# Label the lead variant, since with a single peak the reader should be able to
# read its identity straight off the plot.
lead <- d[which.min(d$P), ]
text(lead$cum, -log10(lead$P), labels = lead$ID, pos = 2, cex = 0.65)
dev.off()

cat("Genome-wide significant (P < 5e-8):", sum(d$P < 5e-8), "\n")
cat("Suggestive (P < 1e-5):             ", sum(d$P < 1e-5), "\n")
cat("Lead variant:", lead$ID, " P =", format(lead$P, digits = 3), "\n")
cat("\nWritten:\n  ", file.path(out_dir, "pdac_demo_04A_manhattan.png"), "\n")
cat("\n=== NEXT STEP ===\n\n")
cat("  Rscript scripts/04A_association_binary/05_summary.R\n\n")
