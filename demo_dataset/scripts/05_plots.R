#!/usr/bin/env Rscript
# =============================================================================
# 05_plots.R
# -----------------------------------------------------------------------------
# Reads the PLINK2 GWAS summary statistics and produces the three standard
# diagnostics of a GWAS:
#   * a Manhattan plot      (where in the genome the signals are)
#   * a QQ plot             (are p-values inflated / well-calibrated?)
#   * the genomic inflation factor lambda_GC (a single-number summary)
#
# The known causal variants (causal_truth.tsv) are highlighted on the Manhattan
# plot so the figure doubles as a teaching aid: students can see which of the
# six simulated loci the GWAS actually recovered at n = 6,000.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

# ---- options -----------------------------------------------------------------
option_list <- list(
  make_option("--results_dir", type = "character", default = "./results",
              help = "directory holding the PLINK2 GWAS output [default %default]"),
  make_option("--truth", type = "character",
              default = "./data/subset/causal_truth.tsv",
              help = "table of simulated causal variants [default %default]"),
  make_option("--out_dir", type = "character", default = "./docs",
              help = "where to write the figures [default %default]"),
  make_option("--sig", type = "double", default = 5e-8,
              help = "genome-wide significance threshold [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- locate and read the summary statistics ---------------------------------
res_file <- list.files(opt$results_dir,
                       pattern = "gwas\\.PHENO\\.glm\\.logistic.*$",
                       full.names = TRUE)[1]
if (is.na(res_file)) stop("No GWAS results found in ", opt$results_dir)
cat("Reading:", res_file, "\n")

gw <- fread(res_file)
setnames(gw, old = names(gw), new = sub("^#", "", names(gw)))  # strip leading #

# PLINK2 column names: CHROM, POS, ID, ..., P  (and TEST when firth-fallback used)
if ("TEST" %in% names(gw)) gw <- gw[TEST == "ADD"]   # keep the additive test row
gw <- gw[!is.na(P) & P > 0]
gw[, CHROM := as.integer(CHROM)]
gw[, POS := as.integer(POS)]

# ---- genomic inflation factor lambda_GC --------------------------------------
# lambda = median(chi^2 observed) / median(chi^2 expected under the null)
chisq  <- qchisq(gw$P, df = 1, lower.tail = FALSE)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
cat(sprintf("\nGenomic inflation factor  lambda_GC = %.3f\n", lambda))
cat("(values near 1.0 indicate well-controlled stratification)\n")

# ---- cumulative genomic position for the Manhattan x-axis --------------------
setorder(gw, CHROM, POS)
chr_len <- gw[, .(maxpos = max(POS)), by = CHROM][order(CHROM)]
chr_len[, offset := cumsum(as.numeric(maxpos)) - maxpos]
gw <- merge(gw, chr_len[, .(CHROM, offset)], by = "CHROM")
gw[, pos_cum := POS + offset]
axis_df <- gw[, .(center = (min(pos_cum) + max(pos_cum)) / 2), by = CHROM]

# ---- flag the true causal variants -------------------------------------------
truth <- if (file.exists(opt$truth)) fread(opt$truth) else data.table(snp = character())
gw[, is_causal := ID %in% truth$snp]

# ---- Manhattan plot ----------------------------------------------------------
gw[, logp := -log10(P)]
manh <- ggplot(gw, aes(pos_cum, logp, colour = factor(CHROM %% 2))) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_point(data = gw[is_causal == TRUE],
             colour = "red", size = 2.2) +
  geom_hline(yintercept = -log10(opt$sig),
             linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c("0" = "#3b6fb6", "1" = "#9fb8d8"),
                      guide = "none") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHROM,
                     expand = c(0.01, 0.01)) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p))),
       title = "Simulated PDAC GWAS — Manhattan plot",
       subtitle = sprintf("red = simulated causal loci;  dashed = 5e-8;  lambda_GC = %.3f",
                          lambda)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(opt$out_dir, "manhattan.png"), manh,
       width = 10, height = 4.5, dpi = 300)

# ---- QQ plot -----------------------------------------------------------------
obs <- -log10(sort(gw$P))
exp <- -log10(ppoints(length(obs)))
qq_df <- data.table(exp = exp, obs = obs)
qq <- ggplot(qq_df, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") +
  geom_point(size = 0.6, colour = "#3b6fb6", alpha = 0.7) +
  labs(x = expression(Expected~-log[10](italic(p))),
       y = expression(Observed~-log[10](italic(p))),
       title = "QQ plot",
       subtitle = sprintf("lambda_GC = %.3f", lambda)) +
  theme_minimal(base_size = 12)

ggsave(file.path(opt$out_dir, "qqplot.png"), qq,
       width = 5, height = 5, dpi = 300)

# ---- short report on the causal loci -----------------------------------------
if (nrow(truth) > 0) {
  hit <- gw[is_causal == TRUE, .(CHROM, ID, POS, P)]
  hit <- merge(hit, truth[, .(snp, gene, OR)],
               by.x = "ID", by.y = "snp", all.x = TRUE)
  hit[, genome_wide := P < opt$sig]
  cat("\nRecovery of the simulated causal loci:\n")
  print(hit[order(P), .(gene, ID, CHROM, OR, P = signif(P, 3), genome_wide)])
  fwrite(hit, file.path(opt$out_dir, "causal_recovery.tsv"), sep = "\t")
}

cat("\nFigures written to", opt$out_dir,
    ":\n  manhattan.png\n  qqplot.png\n  causal_recovery.tsv\n")
cat("\nPipeline complete.\n")
