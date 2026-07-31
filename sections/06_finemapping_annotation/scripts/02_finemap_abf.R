#!/usr/bin/env Rscript

################################################################################
# Section 6: Fine mapping — Step 02: posterior probabilities and credible sets
#
# PURPOSE:
#   Turn a region of association statistics into a ranked shortlist of candidate
#   causal variants.
#
# THE METHOD USED HERE:
#   Wakefield's approximate Bayes factor. For each variant it compares the
#   evidence for an effect against the evidence for no effect, using only the
#   estimate and its standard error, then normalises across the region to give
#   each variant a posterior probability of being the causal one. The variants
#   are ranked by that probability and taken in order until their probabilities
#   sum to 0.95: that is the 95% credible set.
#
#   It is written out in full below because it is short enough to read, and
#   because seeing it makes clear what the method does and does not assume.
#
# THE ASSUMPTION THAT MATTERS:
#   This method assumes there is exactly **one** causal variant in the region.
#   Where a locus carries two independent signals, it will not find them: it
#   returns a single credible set that mixes both, or favours whichever is
#   stronger. Methods that allow several causal variants, such as SuSiE, FINEMAP
#   or CAVIAR, return one credible set per signal and are the right choice when
#   multiplicity is plausible. In PDAC that is not hypothetical: TERT-CLPTM1L
#   carries multiple independent signals in a narrow window.
#
#   The single causal variant assumption is used here because it can be
#   implemented transparently in a few lines, and because at the sample size of
#   this demonstration a second signal would not be detectable anyway. On real
#   data with adequate power, prefer a multi-signal method and report how many
#   sets it returns.
#
# THE PRIOR:
#   W is the prior variance of the effect size on the log odds scale. W = 0.04
#   corresponds to a prior that puts most weight on odds ratios between about
#   0.67 and 1.5, which is the range of genuine cancer susceptibility effects.
#   The choice matters and should be stated: a prior centred on larger effects
#   favours rarer variants, and conversely.
#
# INPUT:
#   - results/assoc/pdac_demo_04A_gwas.*.glm.logistic.hybrid
#   - results/meta/pdac_demo_05_meta.tsv.gz          (optional, for comparison)
#   - results/finemap/abo_eur.bim
#
# OUTPUT:
#   - results/finemap/pdac_demo_06_credible_set.tsv
#   - results/finemap/pdac_demo_06_finemap_summary.txt
#   - results/finemap/credible_set_for_VEP.txt      (default VEP input format)
################################################################################

out_dir <- "results/finemap"
W       <- 0.04          # prior variance of the log odds ratio
COVER   <- 0.95          # credible set coverage
# The generative truth is read from the released truth table rather than
# hardcoded. A constant written here goes stale the moment the phenotype is
# regenerated, and a stale constant of exactly this kind is why this dataset
# was rebuilt.
.truth <- read.delim("demo_data/truth.tsv", stringsAsFactors = FALSE)
.truth <- .truth[grepl(":", .truth$variant), ]
TRUTH   <- .truth$variant[1]   # the ABO variant, this section's target

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

reg <- read.delim(file.path(out_dir, "abo_eur.bim"), header = FALSE,
                  col.names = c("chr", "id", "cm", "pos", "a1", "a2"),
                  stringsAsFactors = FALSE)

read_assoc <- function() {
  f <- Sys.glob("results/assoc/pdac_demo_04A_gwas.*.glm.logistic.hybrid")
  d <- read.table(f[1], header = TRUE, comment.char = "", check.names = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1])
  d <- d[d$TEST == "ADD", ]
  d[is.finite(d$BETA) & is.finite(d$SE) & d$SE > 0, ]
}

# --- Wakefield's approximate Bayes factor ------------------------------------
# For each variant, with estimate b and standard error s:
#   shrinkage  r  = W / (W + s^2)
#   log ABF       = 0.5 * ( log(1 - r) + r * b^2 / s^2 )
# The posterior probability is the ABF normalised over the region, under a flat
# prior that any one variant in the region is equally likely to be causal.
abf_finemap <- function(b, s) {
  r      <- W / (W + s^2)
  lbf    <- 0.5 * (log(1 - r) + r * b^2 / s^2)
  lbf    <- lbf - max(lbf)                 # for numerical stability
  pp     <- exp(lbf) / sum(exp(lbf))
  list(lbf = lbf, pp = pp)
}

credible <- function(pp, coverage = COVER) {
  o   <- order(pp, decreasing = TRUE)
  cs  <- cumsum(pp[o])
  k   <- which(cs >= coverage)[1]
  o[seq_len(k)]
}

# --- European analysis set ---------------------------------------------------
a <- read_assoc()
e <- a[a$ID %in% reg$id, ]
e <- e[order(e$POS), ]
fm <- abf_finemap(e$BETA, e$SE)
e$PP <- fm$pp
cs_idx <- credible(e$PP)
e$in_cs <- seq_len(nrow(e)) %in% cs_idx

# --- the same region in the meta-analysis, if it exists ----------------------
meta_available <- file.exists("results/meta/pdac_demo_05_meta.tsv.gz")
if (meta_available) {
  m <- read.table(gzfile("results/meta/pdac_demo_05_meta.tsv.gz"),
                  header = TRUE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE)
  m <- m[m$ID %in% reg$id & is.finite(m$BETA) & is.finite(m$SE) & m$SE > 0, ]
  m <- m[order(m$POS), ]
  fm_m <- abf_finemap(m$BETA, m$SE)
  m$PP <- fm_m$pp
  cs_m <- credible(m$PP)
  m$in_cs <- seq_len(nrow(m)) %in% cs_m
}

# --- report ------------------------------------------------------------------
con <- file(file.path(out_dir, "pdac_demo_06_finemap_summary.txt"), open = "wt")
w <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = con) }

w("Fine mapping summary, approximate Bayes factor, single causal variant")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("Prior variance W = ", W, ";  coverage = ", COVER)
w("")
w("REGION")
w("  variants analysed: ", nrow(e))
w("  span: chr", e$CHROM[1], ":", min(e$POS), "-", max(e$POS))
w("")
w("EUROPEAN ANALYSIS SET")
w("  variants in the 95% credible set: ", sum(e$in_cs))
w("  top variant: ", e$ID[which.max(e$PP)],
  sprintf("  PP = %.3f", max(e$PP)))
w("  simulated causal variant in the set: ",
  if (TRUTH %in% e$ID[e$in_cs]) "YES" else "NO")
if (TRUTH %in% e$ID) {
  w("  its posterior probability: ", sprintf("%.3f", e$PP[e$ID == TRUTH]),
    "   rank ", which(order(e$PP, decreasing = TRUE) == which(e$ID == TRUTH)))
}
w("")

if (meta_available) {
  w("THREE STRATA COMBINED")
  w("  variants in the 95% credible set: ", sum(m$in_cs))
  w("  top variant: ", m$ID[which.max(m$PP)], sprintf("  PP = %.3f", max(m$PP)))
  w("  simulated causal variant in the set: ",
    if (TRUTH %in% m$ID[m$in_cs]) "YES" else "NO")
  if (TRUTH %in% m$ID) {
    w("  its posterior probability: ", sprintf("%.3f", m$PP[m$ID == TRUTH]),
      "   rank ", which(order(m$PP, decreasing = TRUE) == which(m$ID == TRUTH)))
  }
  w("")
  w("  The comparison is the point of this section. Combining the strata does")
  w("  not merely lower the P value; it changes how well the signal can be")
  w("  localised, because the causal variant sits on different linkage")
  w("  disequilibrium backgrounds in different populations.")
  w("")
}

w("CREDIBLE SET, EUROPEAN ANALYSIS SET")
close(con)

e$OR <- exp(e$BETA)
keep_cols <- intersect(c("CHROM", "POS", "ID", "A1", "A1_FREQ", "OR", "P", "PP"), names(e))
cs <- e[e$in_cs, keep_cols]
cs <- cs[order(-cs$PP), ]
cs$is_causal <- ifelse(cs$ID == TRUTH, "yes", "")
write.table(cs, file.path(out_dir, "pdac_demo_06_credible_set.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# The same set in the default VEP input format, so that anyone moving on to
# functional annotation starts from the credible set this run produced. It is
# written here, beside the set itself, because a file of variants maintained
# separately from the analysis that defines them will eventually describe a
# different analysis.
vep <- data.frame(chr = sub(":.*", "", cs$ID),
                  pos = cs$POS,
                  id  = ".",
                  ref = sapply(strsplit(cs$ID, ":"), `[`, 3),
                  alt = sapply(strsplit(cs$ID, ":"), `[`, 4),
                  q   = ".", filt = ".", info = ".")
write.table(vep, file.path(out_dir, "credible_set_for_VEP.txt"),
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(format(cs, digits = 3),
            file.path(out_dir, "pdac_demo_06_finemap_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)

# keep the per-variant posteriors for the plot
out <- e[, c("CHROM", "POS", "ID", "BETA", "SE", "P", "PP", "in_cs")]
write.table(out, file.path(out_dir, "pdac_demo_06_pp_eur.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (meta_available) {
  write.table(m[, c("CHROM", "POS", "ID", "BETA", "SE", "P", "PP", "in_cs")],
              file.path(out_dir, "pdac_demo_06_pp_meta.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("\nWritten to ", out_dir, "\n", sep = "")
cat("\n=== NEXT STEP ===\n\n")
cat("  Rscript scripts/06_finemapping_annotation/03_finemap_plots.R\n\n")
