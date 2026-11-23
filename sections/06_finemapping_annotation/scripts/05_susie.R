#!/usr/bin/env Rscript

################################################################################
# Section 6, step 05: SuSiE with in-sample LD, as a comparison with the ABF
#
# PURPOSE
#   Step 02 fine-maps ABO with the single-causal-variant approximate Bayes
#   factor (ABF). SuSiE allows up to L causal signals and needs an LD matrix.
#   With individual genotypes the LD matrix can be computed in the analysis
#   set itself (in-sample LD), the setting in which SuSiE is best behaved.
#   Here L = 5: if SuSiE returns one credible set, the multi-signal model
#   agrees with the single-signal one; more than one would be a second signal
#   to explain (the simulation has one causal variant in this region).
#
# METHOD
#   z = BETA / SE of the primary scan (Section 4A) for every variant of the
#   fine-mapping region present in the scan. R = correlation of genotype
#   dosages in the European analysis set, counting the same allele (A1) as
#   the scan so that signs agree; missing calls set to the variant mean; and,
#   as susieR advises when the scan had covariates, the dosages are first
#   residualised on the scan's covariates (SEX, AGE and the components
#   chosen by 04A step 01b), so that R is the LD of the same adjusted data.
#   susie_rss(z, R, n, L = 5), coverage 0.95, residual variance estimated
#   (allowed because the LD is in-sample), other settings at susieR
#   defaults. n is the number of individuals in the scan.
#
# DEPENDENCY
#   susieR is not part of the pipeline toolchain. It is loaded from a separate
#   library, SUSIE_LIB (default ~/R/susie-lib-4.5), so that installing it
#   could not change any package the other steps use. If it cannot be loaded,
#   the step writes pdac_demo_06_susie.txt saying so and exits without error.
#
# INPUT
#   results/finemap/abo_eur.*                 (step 01: region, European set)
#   results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid
#   results/finemap/pdac_demo_06_pp_eur.tsv   (step 02, ABF, for comparison)
#   demo_data/truth.tsv
#
# OUTPUT (results/finemap/, collected into sections/06_finemapping_annotation/results/)
#   pdac_demo_06_susie_pip.tsv      PIP and credible-set membership per variant
#   pdac_demo_06_susie_cs.tsv       one row per credible set
#   pdac_demo_06_susie_summary.tsv  key quantities, side by side with the ABF
#   pdac_demo_06_susie.txt          plain-text report
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

OUT_DIR <- "results/finemap"; REG <- file.path(OUT_DIR, "abo_eur")
GWAS <- "results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid"
COVAR <- "results/assoc/pdac_demo_04A_covar.txt"; NPCF <- "results/assoc/pdac_demo_04A_npc_chosen.txt"
ABF <- file.path(OUT_DIR, "pdac_demo_06_pp_eur.tsv")
L <- 5L; COVERAGE <- 0.95
SUSIE_LIB <- Sys.getenv("SUSIE_LIB", path.expand("~/R/susie-lib-4.5"))
REPORT <- file.path(OUT_DIR, "pdac_demo_06_susie.txt")

if (dir.exists(SUSIE_LIB)) .libPaths(c(SUSIE_LIB, .libPaths()))
if (!suppressWarnings(requireNamespace("susieR", quietly = TRUE))) {
  msg <- c("SECTION 6: SuSiE COMPARISON NOT RUN",
           sprintf("susieR could not be loaded (looked in %s and the default library).", SUSIE_LIB),
           "The pipeline toolchain does not include it; the ABF results of step 02 stand alone.")
  writeLines(msg, REPORT); cat(msg, sep = "\n"); quit(save = "no", status = 0)
}
for (f in c(paste0(REG, ".bed"), GWAS, ABF, COVAR, NPCF, "demo_data/truth.tsv")) if (!file.exists(f)) stop("missing input: ", f)
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
TMP <- file.path(tempdir(), "susie"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

## ---- association statistics in the region ----------------------------------------------
bim <- read.table(paste0(REG, ".bim"), stringsAsFactors = FALSE)
g <- rd(GWAS); if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
g <- g[g$ID %in% bim$V2 & is.finite(g$BETA) & is.finite(g$SE) & g$SE > 0, ]
g <- g[order(g$POS), ]
n <- round(median(g$OBS_CT))

## ---- in-sample LD, counting the scan's A1 allele -------------------------------------------
writeLines(g$ID, file.path(TMP, "ids"))
write.table(g[, c("ID", "A1")], file.path(TMP, "a1"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cmd <- paste("plink2 --bfile", REG, "--extract", file.path(TMP, "ids"), "--export A",
             "--export-allele", file.path(TMP, "a1"), "--out", file.path(TMP, "geno"), "> /dev/null 2>&1")
if (system(cmd) != 0) stop("plink2 export failed: ", cmd)
raw <- read.table(file.path(TMP, "geno.raw"), header = TRUE, check.names = FALSE)
X <- as.matrix(raw[, -(1:6)])
colnames(X) <- sub("_[^_]+$", "", colnames(X))
X <- X[, g$ID]
for (j in seq_len(ncol(X))) { m <- is.na(X[, j]); if (any(m)) X[m, j] <- mean(X[, j], na.rm = TRUE) }
## residualise on the scan's covariates, individuals in the scan only
cv <- rd(COVAR); k <- as.integer(read.table(NPCF)[1, 2])
cv <- cv[match(raw$IID, cv$IID), c("SEX", "AGE", paste0("PC", seq_len(k)))]
inscan <- complete.cases(cv)
X <- X[inscan, , drop = FALSE]; cv <- cv[inscan, ]
X <- apply(X, 2, function(x) residuals(lm(x ~ ., data = cv)))
keep <- apply(X, 2, sd) > 0
X <- X[, keep]; g <- g[keep, ]
R <- cor(X)

## ---- SuSiE-RSS ------------------------------------------------------------------------------------
z <- g$BETA / g$SE
fit <- susieR::susie_rss(z = z, R = R, n = n, L = L, coverage = COVERAGE, estimate_residual_variance = TRUE)
cs <- susieR::susie_get_cs(fit, Xcorr = R, coverage = COVERAGE)
pip <- data.frame(ID = g$ID, POS = g$POS, z = round(z, 4), PIP = fit$pip, cs = NA_integer_)
for (i in seq_along(cs$cs)) pip$cs[cs$cs[[i]]] <- i
causal <- read.delim("demo_data/truth.tsv", stringsAsFactors = FALSE)$variant
causal <- causal[grepl("^9:", causal)][1]
pip <- pip[order(-pip$PIP), ]
cs_tab <- if (length(cs$cs)) data.frame(
  cs = seq_along(cs$cs), size = lengths(cs$cs),
  coverage = round(cs$coverage, 4),
  min_abs_corr = round(cs$purity[, "min.abs.corr"], 4),
  top_id = sapply(cs$cs, function(ix) g$ID[ix][which.max(fit$pip[ix])]),
  contains_causal = sapply(cs$cs, function(ix) causal %in% g$ID[ix])) else
  data.frame(cs = integer(), size = integer(), coverage = numeric(), min_abs_corr = numeric(),
             top_id = character(), contains_causal = logical())

abf <- rd(ABF); abf <- abf[order(-abf$PP), ]
summ <- data.frame(
  quantity = c("susieR_version", "L", "n", "variants", "credible_sets", "cs1_size", "causal_pip", "causal_rank_pip",
               "causal_in_a_cs", "top_pip_id", "top_pip", "converged",
               "abf_credible_set_size", "abf_causal_pp", "abf_causal_rank"),
  value = c(as.character(packageVersion("susieR")), L, n, nrow(g), length(cs$cs),
            if (length(cs$cs)) length(cs$cs[[1]]) else NA,
            signif(pip$PIP[pip$ID == causal], 4), which(pip$ID == causal),
            any(cs_tab$contains_causal), pip$ID[1], signif(pip$PIP[1], 4), isTRUE(fit$converged),
            sum(abf$in_cs), signif(abf$PP[abf$ID == causal], 4), which(abf$ID == causal)))
write.table(pip, file.path(OUT_DIR, "pdac_demo_06_susie_pip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cs_tab, file.path(OUT_DIR, "pdac_demo_06_susie_cs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ, file.path(OUT_DIR, "pdac_demo_06_susie_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rep <- c("SECTION 6: SuSiE WITH IN-SAMPLE LD (comparison with the single-signal ABF)",
         sprintf("susieR %s, L = %d, coverage %.2f, n = %d, %d variants in chr9:132.9-133.6 Mb",
                 packageVersion("susieR"), L, COVERAGE, n, nrow(g)), "",
         sprintf("credible sets: %d", length(cs$cs)),
         if (nrow(cs_tab)) sprintf("  set %d: %d variants, min |r| %.3f, top %s, contains the causal variant: %s",
                                   cs_tab$cs, cs_tab$size, cs_tab$min_abs_corr, cs_tab$top_id, cs_tab$contains_causal),
         sprintf("causal %s: PIP %.3f, rank %d", causal, pip$PIP[pip$ID == causal], which(pip$ID == causal)),
         sprintf("ABF (step 02): credible set %d variants, causal PP %.3f, rank %d",
                 sum(abf$in_cs), abf$PP[abf$ID == causal], which(abf$ID == causal)), "",
         "HOW TO READ THIS",
         "  One credible set means the multi-signal model finds what the single-signal",
         "  model assumed. PIPs and ABF posterior probabilities answer the same question",
         "  under different priors; agreement in the ranking matters more than the values.")
writeLines(rep, REPORT); cat(rep, sep = "\n")
