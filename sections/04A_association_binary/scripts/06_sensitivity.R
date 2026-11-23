#!/usr/bin/env Rscript

################################################################################
# Section 4A, step 06: two sensitivity analyses of the primary scan
#
# PURPOSE
#   Show whether two choices of the pipeline mattered, by changing one choice
#   at a time and repeating everything downstream of it.
#
#   1. RELATEDNESS. The pipeline prunes relatives with a phenotype-aware rule
#      (01B step 07: in a case-control pair the case is kept). The alternative
#      is the PLINK 2 default, --king-cutoff at the same KING kinship 0.0884,
#      which ignores phenotype. Reported: how many cases and controls each rule
#      removes, the resulting European analysis set, lambda, lambda_1000, and
#      the odds ratio and P of the simulated causal ABO variant.
#   2. PRUNING r2 FOR PCA. The pipeline prunes at r2 0.2 before principal
#      component analysis. The alternative is r2 0.1. Reported: variants kept,
#      variance share of the leading components (pooled cohort and within the
#      European set), and lambda and the ABO result of the scan adjusted for
#      the components computed at each threshold.
#
# HOW
#   Each scenario repeats, with the settings of the pipeline itself:
#     01B step 08   MAF 0.01 on the samples kept by the relatedness rule
#     02  step 06   Hardy-Weinberg in controls within ancestry, P < 1e-6,
#                   union of failures excluded
#     02  step 04   European analysis set
#     02  steps 01, 02, 05   pruning (autosomes, MAF 0.05, long-range LD
#                   regions excluded, window 50, step 5, r2 as set), pooled
#                   and within-European PCA, 10 components
#     04A steps 01-03   covariates SEX, AGE and the number of components
#                   chosen by 04A step 01b; PLINK 2 logistic
#                   regression with firth-fallback; lambda as in 03_qq_lambda.R
#   The baseline scenario (phenotype-aware, r2 0.2) must reproduce the
#   canonical outputs of the same run: lambda, the ABO P value and the
#   eigenvalues. If it does not, the script stops, because a sensitivity
#   analysis is only informative when its baseline is the analysis it varies.
#
# INPUT (written earlier in the canonical run)
#   results/qc/pdac_demo_06_filt.*                     before relatedness
#   results/qc/pdac_demo_07_removal_list.txt           phenotype-aware removals
#   results/qc/pdac_demo_07_king_cutoff.king.cutoff.out.id   PLINK 2 default
#   results/pca/pdac_demo_02_pca_all.eigenval, _02_pca_eur.eigenval
#   results/assoc/pdac_demo_04A_lambda.tsv, pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid
#   data_processed/highLD_b38.bed, demo_data/{phenotype,covariates}.txt,
#   demo_data/sample_ancestry.tsv
#
# OUTPUT (results/assoc/, collected into sections/04A_association_binary/results/)
#   pdac_demo_04A_sens_relatedness.tsv
#   pdac_demo_04A_sens_prune_r2.tsv
#   pdac_demo_04A_sens_validation.tsv
#   pdac_demo_04A_sens_summary.txt
#
# No random numbers are drawn; the output is deterministic.
################################################################################

## ---- project root, as in the other R scripts -------------------------------------
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

OUT_DIR   <- "results/assoc"
QC_IN     <- "results/qc/pdac_demo_06_filt"
RM_PHENO  <- "results/qc/pdac_demo_07_removal_list.txt"
RM_KING   <- "results/qc/pdac_demo_07_king_cutoff.king.cutoff.out.id"
LRLD      <- "data_processed/highLD_b38.bed"
PHENO     <- "demo_data/phenotype.txt"
COVAR     <- "demo_data/covariates.txt"
ANC       <- "demo_data/sample_ancestry.tsv"
CAUSAL    <- "9:133273682:A:T"   # demo_data/truth.tsv; risk allele A
MAF_QC    <- 0.01; HWE_P <- 1e-6; PRUNE_MAF <- 0.05
N_PCA_ALL <- 10L   # pooled PCA extracts 10 (02 step 02)
N_PCA_EUR <- 20L   # within-European PCA extracts 20 (02 step 05)
NPCF <- "results/assoc/pdac_demo_04A_npc_chosen.txt"
if (!file.exists(NPCF)) stop("missing input: ", NPCF, " (run 04A step 01b)")
N_PC <- as.integer(read.table(NPCF)[1, 2])   # components in the model (04A step 01b)
TMP <- file.path(tempdir(), "sens"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
for (f in c(paste0(QC_IN, ".bed"), RM_PHENO, RM_KING, LRLD, PHENO, COVAR, ANC,
            "results/pca/pdac_demo_02_pca_all.eigenval", "results/pca/pdac_demo_02_pca_eur.eigenval",
            "results/assoc/pdac_demo_04A_lambda.tsv"))
  if (!file.exists(f)) stop("missing input: ", f)
if (Sys.which("plink2") == "") stop("plink2 not on PATH")

report <- character()
w <- function(...) { line <- paste0(...); report <<- c(report, line); cat(line, "\n", sep = ""); flush.console() }
plink <- function(...) {
  cmd <- paste("plink2", ...)
  if (system(paste(cmd, "> /dev/null 2>&1")) != 0) stop("plink2 failed: ", cmd)
}
rd <- function(f, ...) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, ...)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
ph  <- rd(PHENO); cov <- rd(COVAR)
anc <- read.table(ANC, header = FALSE, col.names = c("IID", "anc"), stringsAsFactors = FALSE)
ids <- function(f) { x <- read.table(f, header = FALSE, comment.char = "#", stringsAsFactors = FALSE); x[[ncol(x)]] }
status_of <- function(iid) ph$PHENO[match(iid, ph$IID)]

## ---- one scenario: everything downstream of the relatedness rule -------------------------
run_chain <- function(tag, removal_file, r2) {
  d <- file.path(TMP, tag); dir.create(d, showWarnings = FALSE)
  p <- function(x) file.path(d, x)
  plink("--bfile", QC_IN, "--remove", removal_file, "--maf", MAF_QC, "--make-bed --out", p("08"))
  fam <- read.table(p("08.fam"), stringsAsFactors = FALSE)
  ## Hardy-Weinberg in controls within each ancestry group, union of failures
  ctrl <- ph$IID[ph$PHENO == 1]
  fail <- character()
  for (g in sort(unique(anc$anc))) {
    k <- anc$IID[anc$anc == g & anc$IID %in% ctrl]
    write.table(cbind(k, k), p(paste0("ctrl_", g)), quote = FALSE, row.names = FALSE, col.names = FALSE)
    plink("--bfile", p("08"), "--keep", p(paste0("ctrl_", g)), "--hardy --out", p(paste0("hwe_", g)))
    h <- rd(p(paste0("hwe_", g, ".hardy")))
    fail <- c(fail, h$ID[h[[ncol(h)]] < HWE_P])
  }
  fail <- sort(unique(fail)); writeLines(fail, p("hwe_exclude"))
  plink("--bfile", p("08"), "--exclude", p("hwe_exclude"), "--make-bed --out", p("hwe_filt"))
  ## European analysis set, from the samples left by this rule
  eur <- fam$V2[fam$V2 %in% anc$IID[tolower(anc$anc) == "eur"]]
  write.table(cbind(eur, eur), p("eur_keep"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  ## pruning and PCA, pooled and within Europeans
  plink("--bfile", p("08"), "--autosome --maf", PRUNE_MAF, "--exclude bed0", LRLD,
        "--indep-pairwise 50 5", r2, "--out", p("prune_all"))
  plink("--bfile", p("08"), "--extract", p("prune_all.prune.in"), "--pca", N_PCA_ALL, "--out", p("pca_all"))
  plink("--bfile", p("08"), "--keep", p("eur_keep"), "--autosome --maf", PRUNE_MAF, "--exclude bed0", LRLD,
        "--indep-pairwise 50 5", r2, "--out", p("prune_eur"))
  plink("--bfile", p("08"), "--keep", p("eur_keep"), "--extract", p("prune_eur.prune.in"),
        "--pca", N_PCA_EUR, "--out", p("pca_eur"))
  ## covariates and association, as 04A steps 01 and 02
  ev <- rd(p("pca_eur.eigenvec"))
  m <- merge(cov, ev, by = intersect(c("FID", "IID"), intersect(names(cov), names(ev))))
  names(m)[1] <- paste0("#", names(m)[1])
  write.table(m, p("covar.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  pcs <- paste0("PC", seq_len(N_PC), collapse = ",")
  plink("--bfile", p("hwe_filt"), "--chr 1-22 --keep", p("eur_keep"), "--pheno", PHENO,
        "--covar", p("covar.txt"), "--covar-name", paste0("SEX,AGE,", pcs),
        "--covar-variance-standardize",
        "--glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err --out", p("gwas"))
  g <- rd(Sys.glob(p("gwas.*.glm.logistic.hybrid"))[1])
  if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
  g <- g[is.finite(g$P) & g$P > 0, ]
  lam <- median(qchisq(1 - g$P, 1)) / qchisq(0.5, 1)
  st <- status_of(eur); nca <- sum(st == 2); nco <- sum(st == 1)
  c1 <- g[g$ID == CAUSAL, ]
  or_risk <- if (nrow(c1)) { or_a1 <- if ("OR" %in% names(c1)) c1$OR else exp(c1$BETA); if (c1$A1 == "A") or_a1 else 1 / or_a1 } else NA
  val_all <- scan(p("pca_all.eigenval"), quiet = TRUE); val_eur <- scan(p("pca_eur.eigenval"), quiet = TRUE)
  rem <- ids(removal_file)
  list(tag = tag, r2 = r2,
       removed = length(rem), removed_cases = sum(status_of(rem) == 2, na.rm = TRUE),
       removed_controls = sum(status_of(rem) == 1, na.rm = TRUE),
       eur_cases = nca, eur_controls = nco, n_eff = 4 / (1 / nca + 1 / nco),
       hwe_excluded = length(fail), variants_tested = nrow(g),
       pruned_all = length(readLines(p("prune_all.prune.in"))),
       pruned_eur = length(readLines(p("prune_eur.prune.in"))),
       lambda = lam, lambda_1000 = 1 + (lam - 1) * (1 / nca + 1 / nco) / (2 / 1000),
       causal_or_risk = or_risk, causal_p = if (nrow(c1)) c1$P else NA,
       pct_all = 100 * val_all / sum(val_all), pct_eur = 100 * val_eur / sum(val_eur),
       val_all = val_all, val_eur = val_eur)
}

w("SECTION 4A: SENSITIVITY OF THE PRIMARY SCAN TO TWO PIPELINE CHOICES")
w("")
base <- run_chain("base", RM_PHENO, 0.2)

## ---- validation: the baseline must be the canonical analysis ------------------------------
lt <- rd("results/assoc/pdac_demo_04A_lambda.tsv")
can_lambda <- as.numeric(lt$value[lt$quantity == "lambda_gc"])
cg <- rd(Sys.glob("results/assoc/pdac_demo_04A_gwas.*.glm.logistic.hybrid")[1])
can_p <- cg$P[cg$ID == CAUSAL & (if ("TEST" %in% names(cg)) cg$TEST == "ADD" else TRUE)]
can_all <- scan("results/pca/pdac_demo_02_pca_all.eigenval", quiet = TRUE)
can_eur <- scan("results/pca/pdac_demo_02_pca_eur.eigenval", quiet = TRUE)
val <- data.frame(
  quantity = c("lambda", "causal ABO P", "max |eigenvalue diff|, pooled", "max |eigenvalue diff|, European"),
  canonical = c(can_lambda, can_p, 0, 0),
  baseline_here = c(base$lambda, base$causal_p,
                    max(abs(base$val_all - can_all)), max(abs(base$val_eur - can_eur))))
val$ok <- c(abs(base$lambda - can_lambda) < 5e-5, abs(log10(base$causal_p) - log10(can_p)) < 1e-6,
            val$baseline_here[3] < 1e-6, val$baseline_here[4] < 1e-6)
write.table(val, file.path(OUT_DIR, "pdac_demo_04A_sens_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w("validation: baseline (phenotype-aware pruning, r2 0.2) against the canonical run")
w(sprintf("  lambda %.4f against %.4f; ABO P %.3g against %.3g; eigenvalues max diff %.2g pooled, %.2g European",
          base$lambda, can_lambda, base$causal_p, can_p, val$baseline_here[3], val$baseline_here[4]))
if (!all(val$ok)) stop("the baseline does not reproduce the canonical run; see pdac_demo_04A_sens_validation.tsv")
w("  reproduced")
w("")

## ---- 1. relatedness rule -------------------------------------------------------------------
king <- run_chain("king_default", RM_KING, 0.2)
row <- function(x, label) data.frame(rule = label, removed = x$removed, removed_cases = x$removed_cases,
  removed_controls = x$removed_controls, eur_cases = x$eur_cases, eur_controls = x$eur_controls,
  n_eff = round(x$n_eff, 1), variants_tested = x$variants_tested, lambda = round(x$lambda, 4),
  lambda_1000 = round(x$lambda_1000, 4), causal_or_risk_A = round(x$causal_or_risk, 3),
  causal_p = signif(x$causal_p, 3))
rel <- rbind(row(base, "phenotype-aware (pipeline)"), row(king, "PLINK 2 --king-cutoff 0.0884 (default)"))
write.table(rel, file.path(OUT_DIR, "pdac_demo_04A_sens_relatedness.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w("1. RELATEDNESS RULE, KING kinship > 0.0884 in both")
for (i in seq_len(nrow(rel))) with(rel[i, ], w(sprintf(
  "  %-40s removed %d (%d cases, %d controls); European set %d cases, %d controls, N_eff %.1f;\n  %-40s lambda %.4f, lambda_1000 %.4f; ABO OR %.3f per A, P %.3g",
  rule, removed, removed_cases, removed_controls, eur_cases, eur_controls, n_eff, "",
  lambda, lambda_1000, causal_or_risk_A, causal_p)))
w("")

## ---- 2. pruning r2 for PCA ------------------------------------------------------------------
r01 <- run_chain("r2_0.1", RM_PHENO, 0.1)
prow <- function(x) data.frame(r2 = x$r2, pruned_variants_pooled = x$pruned_all, pruned_variants_eur = x$pruned_eur,
  pooled_PC1 = round(x$pct_all[1], 1), pooled_PC2 = round(x$pct_all[2], 1), pooled_PC3 = round(x$pct_all[3], 1),
  eur_PC1 = round(x$pct_eur[1], 1), eur_PC2 = round(x$pct_eur[2], 1), eur_PC3 = round(x$pct_eur[3], 1),
  eur_PC4 = round(x$pct_eur[4], 1), eur_PC5 = round(x$pct_eur[5], 1), eur_PC6 = round(x$pct_eur[6], 1),
  lambda = round(x$lambda, 4), causal_or_risk_A = round(x$causal_or_risk, 3), causal_p = signif(x$causal_p, 3))
pr <- rbind(prow(base), prow(r01))
write.table(pr, file.path(OUT_DIR, "pdac_demo_04A_sens_prune_r2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w(sprintf("2. PRUNING r2 BEFORE PCA (variance shares are of the variance captured by the extracted components: %d pooled, %d European)", N_PCA_ALL, N_PCA_EUR))
for (i in seq_len(nrow(pr))) with(pr[i, ], w(sprintf(
  "  r2 %.1f  kept %d pooled, %d European; pooled PC1-3 %.1f/%.1f/%.1f%%; European PC1-6 %.1f/%.1f/%.1f/%.1f/%.1f/%.1f%%;\n          lambda %.4f; ABO OR %.3f per A, P %.3g",
  r2, pruned_variants_pooled, pruned_variants_eur, pooled_PC1, pooled_PC2, pooled_PC3,
  eur_PC1, eur_PC2, eur_PC3, eur_PC4, eur_PC5, eur_PC6, lambda, causal_or_risk_A, causal_p)))
w("")
w("HOW TO READ THIS")
w("  A choice mattered if changing it moves lambda or the ABO result by more than")
w("  the difference between defensible analyses elsewhere in the guide. Each row")
w("  changes one choice and repeats everything downstream of it.")
writeLines(report, file.path(OUT_DIR, "pdac_demo_04A_sens_summary.txt"))
