#!/usr/bin/env Rscript

################################################################################
# Section 4A, step 01b: choose the number of principal components
#
# PURPOSE
#   Replace the default of ten components with the pre-specified, three-part
#   justification the guide recommends: the scree plot, the association of
#   each component with case status, and the stability of lambda and of the
#   leading results across component counts. The number chosen here is used
#   by the primary scan (02_association.sh), by the European stratum of the
#   meta-analysis (05_meta_analysis/scripts/01_stratified_gwas.sh) and by the
#   sensitivity step (06_sensitivity.R).
#
# THE RULE (fixed before the canonical run; our operationalisation of the
# guide's procedure, not a published criterion)
#   Candidates: K = 2, 5, 10, 20 within-European components (02 step 05
#   extracts 20).
#   (a) Scree. e = number of components before the largest ratio between
#       consecutive eigenvalues among the first ten (lambda_j / lambda_j+1).
#       Within one ancestry group the scree is often flat and e small; it
#       sets a floor, not the answer.
#   (b) Association. a = highest index among components associated with case
#       status, each tested alone given SEX and AGE (logistic regression), at
#       a Bonferroni threshold over the 20 tested, P < 0.05 / 20. A nominal
#       P < 0.05 is not used because one in twenty components would pass by
#       chance. a = 0 if none passes.
#   (c) Stability. For each candidate k, the scan is run with SEX, AGE and
#       PC1-PCk. k is stable if, against every larger candidate k': adding
#       the extra components does not lower lambda by more than 0.01
#       (lambda(k') >= lambda(k) - 0.01), the lead variant is the same, and
#       at least 8 of the 10 variants with the smallest P are shared. The
#       criterion is one-sided on purpose: components that remove confounding
#       lower lambda; components that only add parameters leave it unchanged
#       or raise it (the logic of consortia that stopped adding components
#       when inflation no longer fell). A lambda that rises with k is
#       over-adjustment, not instability of k.
#   Choice: the smallest candidate k >= max(e, a) that is stable. If none
#   smaller than 20 is, 20.
#
# HISTORY OF THE RULE (stated because the rule was revised once)
#   The first version used a two-sided criterion, |lambda(k) - lambda(k')| <
#   0.01. Tested on the data of the 23 November run, before the final
#   canonical run, it chose 20 because lambda rose with every component
#   added (1.006, 1.010, 1.022, 1.038 at 2, 5, 10, 20), so only the most
#   adjusted model counted as stable. That inverts the procedure's logic and
#   was corrected to the one-sided criterion above before the final run.
#
# INPUT
#   results/pca/pdac_demo_02_hwe_filt.*, results/pca/pdac_demo_02_eur_keep.txt
#   results/pca/pdac_demo_02_pca_eur.eigenval
#   results/assoc/pdac_demo_04A_covar.txt   (01_make_covariates.R: SEX, AGE, PC1-PC20)
#   demo_data/phenotype.txt
#
# OUTPUT (results/assoc/, collected into sections/04A_association_binary/results/)
#   pdac_demo_04A_npc_scans.tsv       one row per candidate: lambda, lambda_1000,
#                                     genome-wide hits, lead, top-10 overlap
#   pdac_demo_04A_npc_criteria.tsv    e, a, stability per candidate, choice
#   pdac_demo_04A_npc_chosen.txt      one line: chosen_n_pcs<TAB>k (read downstream)
#   pdac_demo_04A_npc_choice.txt      plain-text report with the reasoning
#   pdac_demo_04A_npc_scree.png       within-European scree, 20 components
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
root <- script_dir
for (up in list(c("..", ".."), c("..", "..", ".."))) {
  cand <- normalizePath(do.call(file.path, c(list(script_dir), as.list(up))), mustWork = FALSE)
  if (dir.exists(file.path(cand, "scripts", "dev"))) { root <- cand; break }
}
if (!dir.exists(file.path(root, "scripts", "dev"))) root <- getwd()
setwd(root)

IN       <- "results/pca/pdac_demo_02_hwe_filt"
KEEP     <- "results/pca/pdac_demo_02_eur_keep.txt"
EIGENVAL <- "results/pca/pdac_demo_02_pca_eur.eigenval"
COVAR    <- "results/assoc/pdac_demo_04A_covar.txt"
PHENO    <- "demo_data/phenotype.txt"
OUT_DIR  <- "results/assoc"
CANDIDATES <- c(2L, 5L, 10L, 20L)
N_TESTED <- 20L
ALPHA_ASSOC <- 0.05 / N_TESTED
LAMBDA_TOL <- 0.01
TOP_N <- 10L; TOP_MIN_SHARED <- 8L
TMP <- file.path(tempdir(), "npc"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
for (f in c(paste0(IN, ".bed"), KEEP, EIGENVAL, COVAR, PHENO)) if (!file.exists(f)) stop("missing input: ", f)
if (Sys.which("plink2") == "") stop("plink2 not on PATH")

report <- character()
w <- function(...) { line <- paste0(...); report <<- c(report, line); cat(line, "\n", sep = "") }
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }

cov <- rd(COVAR); ph <- rd(PHENO)
pcs_avail <- grep("^PC[0-9]+$", names(cov), value = TRUE)
if (length(pcs_avail) < max(CANDIDATES))
  stop(COVAR, " has ", length(pcs_avail), " components; 02 step 05 must extract ", max(CANDIDATES))

## ---- (a) scree --------------------------------------------------------------------------
ev <- scan(EIGENVAL, quiet = TRUE)
if (length(ev) < N_TESTED) stop(EIGENVAL, " has ", length(ev), " eigenvalues, need ", N_TESTED)
ratio <- ev[1:10] / ev[2:11]
e <- which.max(ratio)
png(file.path(OUT_DIR, "pdac_demo_04A_npc_scree.png"), width = 1600, height = 1000, res = 200)
plot(seq_along(ev), 100 * ev / sum(ev), type = "b", pch = 19, col = "#0072B2",
     xlab = "Component", ylab = "% of variance captured by the extracted components",
     main = "Within-European scree, 20 components")
abline(v = CANDIDATES, lty = 3, col = "grey50")
invisible(dev.off())

## ---- (b) association of each component with case status -----------------------------------------
d <- merge(cov, ph[, c("IID", "PHENO")], by = "IID")
d$status <- ifelse(d$PHENO == 2, 1, ifelse(d$PHENO == 1, 0, NA))
d <- d[!is.na(d$status), ]
pa <- sapply(pcs_avail[1:N_TESTED], function(pc)
  summary(glm(as.formula(paste("status ~ SEX + AGE +", pc)), data = d, family = binomial()))$coefficients[pc, 4])
assoc_idx <- which(pa < ALPHA_ASSOC)
a <- if (length(assoc_idx)) max(assoc_idx) else 0L

## ---- (c) scans at each candidate count --------------------------------------------------------------
lam <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
scans <- list()
for (k in CANDIDATES) {
  out <- file.path(TMP, paste0("k", k))
  cmd <- paste("plink2 --bfile", IN, "--chr 1-22 --keep", KEEP, "--pheno", PHENO, "--covar", COVAR,
               "--covar-name", paste0("SEX,AGE,", paste0("PC", 1:k, collapse = ",")),
               "--covar-variance-standardize",
               "--glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err --out", out, "> /dev/null 2>&1")
  if (system(cmd) != 0) stop("plink2 failed: ", cmd)
  g <- rd(Sys.glob(paste0(out, ".*.glm.logistic.hybrid"))[1])
  if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
  g <- g[is.finite(g$P) & g$P > 0, ]
  g <- g[order(g$P), ]
  scans[[as.character(k)]] <- g
}
nca <- sum(d$status == 1 & d$IID %in% read.table(KEEP)[[2]]); nco <- sum(d$status == 0 & d$IID %in% read.table(KEEP)[[2]])
top <- function(k) head(scans[[as.character(k)]]$ID, TOP_N)
tab <- do.call(rbind, lapply(CANDIDATES, function(k) {
  g <- scans[[as.character(k)]]; l <- lam(g$P)
  data.frame(n_pcs = k, variants = nrow(g), lambda = round(l, 4),
             lambda_1000 = round(1 + (l - 1) * (1 / nca + 1 / nco) / (2 / 1000), 4),
             n_p_lt_5e8 = sum(g$P < 5e-8), n_p_lt_1e5 = sum(g$P < 1e-5),
             lead_id = g$ID[1], lead_p = signif(g$P[1], 3),
             top10_shared_with_20 = length(intersect(top(k), top(max(CANDIDATES)))))
}))
stable <- sapply(CANDIDATES, function(k) {
  larger <- CANDIDATES[CANDIDATES > k]
  if (!length(larger)) return(TRUE)
  all(sapply(larger, function(k2)
    tab$lambda[tab$n_pcs == k2] >= tab$lambda[tab$n_pcs == k] - LAMBDA_TOL &&
      tab$lead_id[tab$n_pcs == k] == tab$lead_id[tab$n_pcs == k2] &&
      length(intersect(top(k), top(k2))) >= TOP_MIN_SHARED))
})
floor_k <- max(e, a)
ok <- CANDIDATES[CANDIDATES >= floor_k & stable]
chosen <- if (length(ok)) min(ok) else max(CANDIDATES)

crit <- data.frame(n_pcs = CANDIDATES, at_or_above_floor = CANDIDATES >= floor_k, stable = stable,
                   chosen = CANDIDATES == chosen)
write.table(tab, file.path(OUT_DIR, "pdac_demo_04A_npc_scans.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(data.frame(quantity = c("scree_floor_e", "association_floor_a", "association_alpha",
                                          "lambda_tolerance", "top_n", "top_min_shared", "chosen_n_pcs"),
                             value = c(e, a, ALPHA_ASSOC, LAMBDA_TOL, TOP_N, TOP_MIN_SHARED, chosen)),
                  data.frame(quantity = paste0("pc", 1:N_TESTED, "_assoc_p"), value = signif(pa, 4))),
            file.path(OUT_DIR, "pdac_demo_04A_npc_criteria.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(paste0("chosen_n_pcs\t", chosen), file.path(OUT_DIR, "pdac_demo_04A_npc_chosen.txt"))

w("SECTION 4A: NUMBER OF PRINCIPAL COMPONENTS")
w("")
w(sprintf("European analysis set: %d cases, %d controls", nca, nco))
w(sprintf("(a) scree: largest ratio of consecutive eigenvalues among the first ten is after PC%d (ratio %.2f); floor e = %d",
          e, ratio[e], e))
w(sprintf("    variance shares of the 20 extracted components, PC1-PC6: %s %%",
          paste(sprintf("%.1f", 100 * ev[1:6] / sum(ev)), collapse = ", ")))
w(sprintf("(b) association with case status given SEX and AGE, threshold P < %.4f: %s; floor a = %d",
          ALPHA_ASSOC, if (length(assoc_idx)) paste(sprintf("PC%d (P %.2g)", assoc_idx, pa[assoc_idx]), collapse = ", ") else "none",
          a))
w(sprintf("    nominal P < 0.05, not used for the rule: %s",
          if (any(pa < 0.05)) paste(sprintf("PC%d", which(pa < 0.05)), collapse = ", ") else "none"))
w("(c) scans:")
for (i in seq_len(nrow(tab))) with(tab[i, ], w(sprintf(
  "    %2d PCs  lambda %.4f  lambda_1000 %.4f  P<5e-8 %d  P<1e-5 %d  lead %s (P %.3g)  top-10 shared with 20 PCs %d  stable %s",
  n_pcs, lambda, lambda_1000, n_p_lt_5e8, n_p_lt_1e5, lead_id, lead_p, top10_shared_with_20, stable[i])))
w("")
w(sprintf("CHOSEN: %d components (smallest candidate >= max(e, a) = %d that is stable: no larger candidate lowers lambda by more than %.2f or changes the lead or the top 10)",
          chosen, floor_k, LAMBDA_TOL))
writeLines(report, file.path(OUT_DIR, "pdac_demo_04A_npc_choice.txt"))
