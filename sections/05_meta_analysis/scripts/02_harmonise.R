################################################################################
# Section 5: Meta-analysis — Step 02: harmonising the summary statistics
#
# PURPOSE:
#   Bring the per-stratum results onto a common footing before combining them.
#   This step does no statistics. It exists because the commonest errors in
#   meta-analysis are bookkeeping errors, and they are silent: the numbers still
#   combine, the plot still looks reasonable, and the effect is in the wrong
#   direction.
#
# THE FOUR THINGS THAT GO WRONG:
#
#   1. Effect allele mismatch. Study A reports the effect of allele G, study B
#      the effect of allele A at the same variant. Combining them without
#      flipping the sign of one estimate cancels a real signal, or manufactures
#      one. This is the single most common error.
#
#   2. Strand ambiguity. At A/T and C/G variants, the effect allele cannot be
#      resolved by comparing alleles alone, because the complement of A is T and
#      of C is G. Frequency can sometimes resolve it, but not near 50%. The
#      conservative treatment, adopted here, is to flag them; many consortia
#      simply exclude ambiguous variants with frequency between 0.4 and 0.6.
#
#   3. Effect-size scale. A log odds ratio and an odds ratio cannot be averaged
#      together, and neither can a beta from a linear model of a binary trait.
#      Confirm what each contributing study actually reported.
#
#   4. Sample overlap. If the same individuals appear in two contributing
#      studies, the estimates are correlated and a standard fixed-effect
#      meta-analysis understates the standard error, inflating significance.
#      There is no way to detect this from summary statistics alone; it has to
#      be established from the study descriptions.
#
# IN THIS DEMONSTRATION:
#   All three strata were derived from a single genotype file, so they share one
#   variant set, one allele coding and one genome build, and no individual
#   appears twice. Every check below therefore passes trivially. That is exactly
#   why the checks are written out in full: in a real meta-analysis of
#   independently genotyped cohorts, this is the step where the errors are, and
#   a script that assumes agreement will not find them.
#
# INPUT:
#   - results/meta/<group>/pdac_demo_05_<group>_gwas.PHENO.glm.logistic.hybrid
#
# OUTPUT:
#   - results/meta/pdac_demo_05_harmonised.tsv.gz
#   - results/meta/pdac_demo_05_harmonisation_report.txt
################################################################################

out_dir <- "results/meta"
strata  <- c("eur", "afr", "eas")
ref     <- "eur"   # the stratum whose effect allele defines the reference

report <- file(file.path(out_dir, "pdac_demo_05_harmonisation_report.txt"), open = "wt")
w <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = report) }

read_stratum <- function(g) {
  f <- Sys.glob(file.path(out_dir, g, sprintf("pdac_demo_05_%s_gwas.*.glm.logistic.hybrid", g)))
  if (!length(f)) stop("No summary statistics found for stratum ", g)
  d <- read.table(f[1], header = TRUE, comment.char = "", check.names = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1])
  if ("TEST" %in% names(d)) d <- d[d$TEST == "ADD", ]
  d <- d[is.finite(d$BETA) & is.finite(d$SE) & d$SE > 0 & is.finite(d$P), ]
  d[, c("CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "OBS_CT", "BETA", "SE", "P")]
}

w("Harmonisation report")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("")

dat <- list()
for (g in strata) {
  d <- read_stratum(g)
  dat[[g]] <- d
  w(sprintf("%-4s  variants with a usable estimate: %d", g, nrow(d)))
}
w("")

# --- 1. variants shared by all strata ----------------------------------------
common <- Reduce(intersect, lapply(dat, function(d) d$ID))
w("Variants present in all ", length(strata), " strata: ", length(common))
w("  (variants tested in only some strata are dropped here for simplicity;")
w("   a real meta-analysis usually keeps them, meta-analysing each variant")
w("   across whatever studies carry it, and reports the study count per variant)")
w("")

base <- dat[[ref]][match(common, dat[[ref]]$ID), ]

# --- 2. allele consistency and orientation -----------------------------------
n_flip <- setNames(integer(length(strata)), strata)
n_mismatch <- setNames(integer(length(strata)), strata)

for (g in strata) {
  d <- dat[[g]][match(common, dat[[g]]$ID), ]

  # Same allele pair at the same variant?
  same_pair <- (d$REF == base$REF & d$ALT == base$ALT)
  n_mismatch[g] <- sum(!same_pair)

  # Is the effect allele the same as in the reference stratum?
  flip <- d$A1 != base$A1
  n_flip[g] <- sum(flip)

  # Align: flip the sign of the effect and complement the frequency.
  d$BETA[flip]    <- -d$BETA[flip]
  d$A1_FREQ[flip] <- 1 - d$A1_FREQ[flip]
  d$A1            <- base$A1

  dat[[g]] <- d
}

w("Allele checks, against the ", ref, " stratum as reference:")
for (g in strata) {
  w(sprintf("  %-4s  allele-pair mismatches: %d   effect-allele flips applied: %d",
            g, n_mismatch[g], n_flip[g]))
}
w("")

# --- 3. strand-ambiguous variants --------------------------------------------
amb <- (base$REF == "A" & base$ALT == "T") | (base$REF == "T" & base$ALT == "A") |
       (base$REF == "C" & base$ALT == "G") | (base$REF == "G" & base$ALT == "C")
amb_mid <- amb & base$A1_FREQ > 0.4 & base$A1_FREQ < 0.6

w("Strand-ambiguous variants (A/T or C/G): ", sum(amb),
  " of ", length(common), sprintf(" (%.1f%%)", 100 * sum(amb) / length(common)))
w("  of which frequency is between 0.4 and 0.6, where frequency cannot resolve")
w("  the strand: ", sum(amb_mid))
w("  These are flagged, not removed. Because all strata here come from one")
w("  genotype file the strand is known to agree; across independently")
w("  genotyped cohorts, many groups exclude the ambiguous variants outright.")
w("")

# --- 4. frequency agreement --------------------------------------------------
# Large frequency differences between studies of the same ancestry indicate a
# harmonisation problem. Between ancestries they are expected and are precisely
# what a trans-ancestry meta-analysis has to accommodate.
w("Allele-frequency correlation between strata (Pearson r):")
for (i in seq_along(strata)) {
  for (j in seq_along(strata)) {
    if (i < j) {
      r <- cor(dat[[strata[i]]]$A1_FREQ, dat[[strata[j]]]$A1_FREQ, use = "complete.obs")
      w(sprintf("  %s vs %s: %.3f", strata[i], strata[j], r))
    }
  }
}
w("  Values well below 1 are expected here: these are different continental")
w("  populations. Between cohorts of the same ancestry, anything below about")
w("  0.9 would indicate a problem rather than biology.")
w("")

# --- 5. assemble the harmonised table ----------------------------------------
h <- data.frame(
  ID    = common,
  CHROM = base$CHROM,
  POS   = base$POS,
  REF   = base$REF,
  ALT   = base$ALT,
  A1    = base$A1,
  strand_ambiguous = amb,
  stringsAsFactors = FALSE
)
for (g in strata) {
  h[[paste0("BETA_", g)]] <- dat[[g]]$BETA
  h[[paste0("SE_",   g)]] <- dat[[g]]$SE
  h[[paste0("P_",    g)]] <- dat[[g]]$P
  h[[paste0("FREQ_", g)]] <- dat[[g]]$A1_FREQ
  h[[paste0("N_",    g)]] <- dat[[g]]$OBS_CT
}

gz <- gzfile(file.path(out_dir, "pdac_demo_05_harmonised.tsv.gz"), "wt")
write.table(h, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

w("Harmonised table written: ", nrow(h), " variants x ", length(strata), " strata")
close(report)

cat("\nWritten:\n")
cat("  ", file.path(out_dir, "pdac_demo_05_harmonised.tsv.gz"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_05_harmonisation_report.txt"), "\n")
cat("\n=== NEXT STEP ===\n\n")
cat("  Rscript scripts/05_meta_analysis/03_meta_analysis.R\n\n")
