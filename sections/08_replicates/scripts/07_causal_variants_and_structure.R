#!/usr/bin/env Rscript

################################################################################
# Section 8, step 07: do the causal genotypes follow the within-European
# components, and does structure explain the ABO x TERT/CLPTM1L correlation?
#
# PURPOSE
#   In the generator case status depends on ancestry group only through the
#   intercept, so within a group it is independent of structure GIVEN the
#   causal variants. A single-variant scan does not condition on them: if the
#   frequency of a causal variant follows a within-group axis, case status
#   follows that axis too, and the components are then confounders for every
#   other variant on it. This step measures how far that is the case here.
#
# WHAT IT DOES (European analysis set, non-missing genotypes, risk-allele count)
#   1. Each causal genotype regressed on PC1..PCk, k = 2, 5, 10, 20: R2 and
#      the P of the F test; and on each of the 20 components alone.
#   2. Correlation between the causal genotypes, and partial correlation
#      given PC1..PCk, k = 2, 10, 20 (residuals of each genotype on the
#      components; t test with n - 2 - k degrees of freedom).
#   Nothing here depends on the phenotype: one answer for all replicates.
#
# NOT PART OF THE CANONICAL RUN: reads its outputs in results/.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_causal_structure.tsv        locus x k: R2, P
#   pdac_demo_08_causal_structure_by_pc.tsv  locus x component: r, P
#   pdac_demo_08_causal_partial_r.tsv        pair x k: r, P
#   pdac_demo_08_causal_structure.txt        plain-text report
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(normalizePath(file.path(script_dir, "..", "..", "..")))
OUT <- "sections/08_replicates/results"
IN <- "results/pca/pdac_demo_02_hwe_filt"
CAUSAL <- c(G1 = "9:133273682:A:T", G2 = "5:1286401:C:A", G3 = "13:73362027:C:T")
LABEL <- c(G1 = "ABO", G2 = "TERT/CLPTM1L", G3 = "13q22.1")
K_R2 <- c(2L, 5L, 10L, 20L); K_PART <- c(0L, 2L, 10L, 20L)
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
tmp <- tempfile("cs"); dir.create(tmp)
writeLines(unname(CAUSAL), file.path(tmp, "snps"))
if (system2("plink2", c("--bfile", IN, "--extract", file.path(tmp, "snps"), "--export", "A", "--out", file.path(tmp, "c")),
            stdout = FALSE, stderr = FALSE) != 0) stop("plink2 export failed")
system2("plink2", c("--bfile", "demo_data/pdac_demo", "--extract", file.path(tmp, "snps"), "--recode", "A",
                    "--out", file.path(tmp, "gen")), stdout = FALSE, stderr = FALSE)
gen_hdr <- names(read.table(file.path(tmp, "gen.raw"), header = TRUE, check.names = FALSE, nrows = 1))
raw <- read.table(file.path(tmp, "c.raw"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
col_of <- function(h, id) grep(paste0("^", gsub("([.:])", "\\\\\\1", id), "_"), h, value = TRUE)
keep <- read.table("results/pca/pdac_demo_02_eur_keep.txt", stringsAsFactors = FALSE)[[2]]
cv <- rd("results/assoc/pdac_demo_04A_covar.txt")
G <- sapply(CAUSAL, function(id) { cn <- col_of(names(raw), id); x <- raw[[cn]]
  if (sub(".*_", "", cn) != sub(".*_", "", col_of(gen_hdr, id))) 2 - x else x })[match(keep, raw$IID), ]
PC <- as.matrix(cv[match(keep, cv$IID), paste0("PC", 1:20)])
if (anyNA(PC)) stop("components missing")

r2 <- do.call(rbind, lapply(names(CAUSAL), function(L) do.call(rbind, lapply(K_R2, function(k) {
  ok <- !is.na(G[, L]); m <- summary(lm(G[ok, L] ~ PC[ok, 1:k])); f <- m$fstatistic
  data.frame(locus = L, label = LABEL[[L]], n_pcs = k, n = sum(ok), risk_allele_freq = mean(G[ok, L]) / 2, r2 = m$r.squared,
             adj_r2 = m$adj.r.squared, f = f[[1]], p = pf(f[[1]], f[[2]], f[[3]], lower.tail = FALSE)) }))))
bypc <- do.call(rbind, lapply(names(CAUSAL), function(L) do.call(rbind, lapply(1:20, function(j) {
  ok <- !is.na(G[, L]); ct <- cor.test(G[ok, L], PC[ok, j])
  data.frame(locus = L, label = LABEL[[L]], component = j, r = unname(ct$estimate), p = ct$p.value) }))))
pairs <- list(c("G1", "G2"), c("G1", "G3"), c("G2", "G3"))
part <- do.call(rbind, lapply(pairs, function(pp) do.call(rbind, lapply(K_PART, function(k) {
  ok <- complete.cases(G[, pp]); a <- G[ok, pp[1]]; b <- G[ok, pp[2]]
  if (k > 0) { a <- resid(lm(a ~ PC[ok, 1:k])); b <- resid(lm(b ~ PC[ok, 1:k])) }
  r <- cor(a, b); df <- sum(ok) - 2 - k; t <- r * sqrt(df / (1 - r^2))
  data.frame(pair = paste0(LABEL[[pp[1]]], " x ", LABEL[[pp[2]]]), given_n_pcs = k, n = sum(ok), r = r, p = 2 * pt(-abs(t), df)) }))))
sg <- function(x) { num <- sapply(x, is.numeric); x[num] <- lapply(x[num], function(v) signif(v, 4)); x }
wr <- function(x, f) write.table(sg(x), file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
wr(r2, "pdac_demo_08_causal_structure.tsv"); wr(bypc, "pdac_demo_08_causal_structure_by_pc.tsv"); wr(part, "pdac_demo_08_causal_partial_r.tsv")

R <- c("SECTION 8: CAUSAL GENOTYPES AND THE WITHIN-EUROPEAN COMPONENTS",
       sprintf("European analysis set, %d individuals; risk-allele counts, non-missing genotypes", length(keep)), "",
       "1. GENOTYPE REGRESSED ON PC1..PCk: R2 (P of the F test)")
for (L in names(CAUSAL)) { x <- r2[r2$locus == L, ]; b <- bypc[bypc$locus == L, ]; s <- b[order(b$p), ][1:3, ]
  R <- c(R, sprintf("  %-13s n %d  %s", LABEL[[L]], x$n[1], paste(sprintf("k %d: %.4f (%.2g)", x$n_pcs, x$r2, x$p), collapse = "   ")),
         sprintf("  %-13s single components with the smallest P: %s; at P < 0.05: %d of 20", "",
                 paste(sprintf("PC%d r %+.3f (%.2g)", s$component, s$r, s$p), collapse = ", "), sum(b$p < 0.05))) }
R <- c(R, "", "2. CORRELATION BETWEEN CAUSAL GENOTYPES, PLAIN AND GIVEN PC1..PCk: r (P)")
for (pn in unique(part$pair)) { x <- part[part$pair == pn, ]
  R <- c(R, sprintf("  %-26s %s", pn, paste(sprintf("k %d: %+.4f (%.2g)", x$given_n_pcs, x$r, x$p), collapse = "   "))) }
R <- c(R, "", "HOW TO READ THIS",
  "  A causal genotype that follows the components makes case status follow",
  "  them, although the generator gives structure no effect of its own. The",
  "  number of components that is correct for the single-variant scan is then",
  "  not zero for variants on the same axes; how much it matters is measured",
  "  by step 06 (lambda with 0 against 2 components, variants following",
  "  PC1 and PC2). A correlation between two causal genotypes that the",
  "  components do not remove is not structure captured by them.")
writeLines(R, file.path(OUT, "pdac_demo_08_causal_structure.txt")); cat(R, sep = "\n")
