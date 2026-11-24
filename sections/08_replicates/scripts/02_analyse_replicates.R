#!/usr/bin/env Rscript

################################################################################
# Section 8, step 02: analyse every replicate as the canonical run would
#
# For each replicate phenotype of step 01:
#   - the European scan of Section 4A: logistic regression with firth-fallback,
#     SEX, AGE and the components chosen by 04A step 01b, European analysis set;
#   - the three ancestry strata of Section 5, each with its own components,
#     combined by fixed-effect inverse-variance meta-analysis aligned to the
#     European A1, with Cochran's Q, I-squared and DerSimonian-Laird random
#     effects, as 05_meta_analysis/scripts/03_meta_analysis.R computes them;
#   - the single-causal-variant approximate Bayes factor (W = 0.04, 95%
#     coverage, as 06 step 02) in a window of +/- 350 kb around each causal
#     variant, in the European scan and in the meta-analysis.
#
# WHAT IS HELD FIXED, AND WHY (simplifications, stated)
#   Samples, quality control, strata, components and the number of components
#   are those of the canonical run (results/ after scripts/dev/run_all.sh).
#   Two QC steps of the pipeline depend on the phenotype: relatedness pruning
#   prefers cases, and Hardy-Weinberg is tested in controls. They are not
#   repeated per replicate: the replicate samples and variants are the
#   canonical ones. Age is the replicate's own (the generator draws it).
#   Variants: the LD-pruned European set used for its components (enough for
#   lambda, as agreed) plus every variant within 350 kb of each causal
#   variant (for power, effect sizes and fine-mapping).
#
# HISTORY
#   Until 24 November 2026 (results of 01e7f98) the European scan passed all
#   20 component columns of pdac_demo_04A_covar.txt: k was read and not used.
#   Found by the manuscript session (checks/diag_replicates_20261124). Fixed:
#   the scan now uses PC1..PCk; the strata of the meta-analysis use the
#   components of their own covar.txt, as before.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_per_replicate.tsv   one row per replicate x locus x analysis
#   pdac_demo_08_lambda.tsv          one row per replicate
################################################################################

suppressPackageStartupMessages(library(parallel))
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "..", "..", ".."))
setwd(root)

OUT <- "sections/08_replicates/results"; SIM <- file.path(OUT, "sim")
IN <- "results/pca/pdac_demo_02_hwe_filt"
EUR_KEEP <- "results/pca/pdac_demo_02_eur_keep.txt"
EUR_COV <- "results/assoc/pdac_demo_04A_covar.txt"
NPCF <- "results/assoc/pdac_demo_04A_npc_chosen.txt"
PRUNE <- "results/pca/pdac_demo_02_prune_eur.prune.in"
STRATA <- c("eur", "afr", "eas")
CAUSAL <- c(G1 = "9:133273682:A:T", G2 = "5:1286401:C:A", G3 = "13:73362027:C:T")
OR_GEN <- c(G1 = 2.60, G2 = 1.50, G3 = 1.97)
WIN <- 350000L; W <- 0.04; COVERAGE <- 0.95; ALPHA <- 5e-8
N_PAR <- as.integer(Sys.getenv("GWAS_08_PARALLEL", "6"))
THREADS <- as.integer(Sys.getenv("GWAS_08_THREADS", "4"))
for (f in c(paste0(IN, ".bed"), EUR_KEEP, EUR_COV, NPCF, PRUNE,
            file.path("results/meta", STRATA, "keep.txt"), file.path("results/meta", STRATA, "covar.txt")))
  if (!file.exists(f)) stop("missing input from the canonical run: ", f)
reps <- sort(Sys.glob(file.path(SIM, "rep_*.pheno")))
if (!length(reps)) stop("no replicates in ", SIM, ": run step 01")
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
k <- as.integer(read.table(NPCF)[1, 2])

## ---- fixed inputs ------------------------------------------------------------------------
tmp <- file.path(tempdir(), "rep"); dir.create(tmp, showWarnings = FALSE)
bim <- read.table(paste0(IN, ".bim"), col.names = c("CHR", "ID", "CM", "POS", "A1", "A2"), stringsAsFactors = FALSE)
cpos <- bim[match(CAUSAL, bim$ID), c("CHR", "POS")]
if (anyNA(cpos$POS)) stop("a causal variant is missing from ", IN)
in_win <- Reduce(`|`, lapply(seq_along(CAUSAL), function(i) bim$CHR == cpos$CHR[i] & abs(bim$POS - cpos$POS[i]) <= WIN))
pruned <- readLines(PRUNE)
V <- unique(c(pruned, bim$ID[in_win]))
writeLines(V, file.path(tmp, "extract"))
# counted allele of the generator (plink2 --recode A on the released genotypes)
writeLines(unname(CAUSAL), file.path(tmp, "c.snps"))
system2("plink2", c("--bfile", "demo_data/pdac_demo", "--extract", file.path(tmp, "c.snps"), "--recode", "A",
                    "--out", file.path(tmp, "c")), stdout = FALSE, stderr = FALSE)
hdr <- names(read.table(file.path(tmp, "c.raw"), header = TRUE, check.names = FALSE, nrows = 1))
counted <- sapply(CAUSAL, function(id) sub(".*_", "", grep(paste0("^", gsub("([.:])", "\\\\\\1", id), "_"), hdr, value = TRUE)))
cov_e <- rd(EUR_COV)
# the European scan of Section 4A adjusts for SEX, AGE and PC1..PCk, k chosen by
# 04A step 01b; the covariate file carries all 20 extracted components
pc_all <- grep("^PC[0-9]+$", names(cov_e), value = TRUE)
if (k < 1 || k > length(pc_all)) stop("chosen number of components ", k, " outside 1..", length(pc_all))
cov_e <- cov_e[, c(setdiff(names(cov_e), pc_all), paste0("PC", seq_len(k)))]
cov_s <- lapply(setNames(STRATA, STRATA), function(g) rd(file.path("results/meta", g, "covar.txt")))

lam <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
abf_cs <- function(b, s, ids, causal) {
  ok <- is.finite(b) & is.finite(s) & s > 0; b <- b[ok]; s <- s[ok]; ids <- ids[ok]
  V <- s^2; r <- W / (V + W); z <- b / s
  lbf <- 0.5 * (log(1 - r) + r * z^2)
  pp <- exp(lbf - max(lbf)); pp <- pp / sum(pp)
  o <- order(-pp); n_cs <- which(cumsum(pp[o]) >= COVERAGE)[1]
  rk <- which(ids[o] == causal)
  c(cs_size = n_cs, causal_in_cs = length(rk) == 1 && rk <= n_cs, causal_rank = if (length(rk)) rk else NA,
    causal_pp = if (length(rk)) pp[o][rk] else NA, lead_is_causal = ids[o][1] == causal)
}
glm_run <- function(keep, cov, pheno, age, tag) {
  cv <- cov; cv$AGE <- age$AGE[match(cv$IID, age$IID)]
  pcs <- grep("^PC[0-9]+$", names(cv), value = TRUE)
  cf <- file.path(tmp, paste0(tag, ".cov")); names(cv)[1] <- paste0("#", names(cv)[1])
  write.table(cv, cf, sep = "\t", quote = FALSE, row.names = FALSE)
  out <- file.path(tmp, tag)
  st <- system2("plink2", c("--bfile", IN, "--extract", file.path(tmp, "extract"), "--chr", "1-22", "--keep", keep,
                            "--pheno", pheno, "--covar", cf, "--covar-name", paste(c("SEX", "AGE", pcs), collapse = ","),
                            "--covar-variance-standardize", "--threads", THREADS,
                            "--glm", "firth-fallback", "hide-covar", "cols=+a1freq,+beta,+orbeta,+nobs,+err",
                            "--out", out), stdout = FALSE, stderr = FALSE)
  if (st != 0) stop("plink2 failed for ", tag)
  r <- rd(Sys.glob(paste0(out, ".*.glm.logistic.hybrid"))[1]); if ("TEST" %in% names(r)) r <- r[r$TEST == "ADD", ]
  unlink(Sys.glob(paste0(out, "*")))
  r
}
cc <- function(pheno, keep) { p <- rd(pheno); k <- read.table(keep)[[2]]; s <- p$PHENO[p$IID %in% k]; c(sum(s == 2), sum(s == 1)) }

one <- function(f) {
  i <- as.integer(sub(".*rep_([0-9]+)\\.pheno$", "\\1", f)); tag <- sprintf("r%03d", i)
  age <- read.delim(sub("\\.pheno$", ".age", f))
  e <- glm_run(EUR_KEEP, cov_e, f, age, paste0(tag, "_4a"))
  s <- lapply(STRATA, function(g) glm_run(file.path("results/meta", g, "keep.txt"), cov_s[[g]], f, age, paste0(tag, "_", g)))
  names(s) <- STRATA
  ## meta-analysis on variants present in all strata, aligned to the European A1
  ids <- Reduce(intersect, lapply(s, function(d) d$ID[is.finite(d$BETA) & is.finite(d$SE)]))
  base <- s$eur[match(ids, s$eur$ID), ]
  B <- S <- matrix(NA_real_, length(ids), 3, dimnames = list(ids, STRATA))
  for (g in STRATA) { d <- s[[g]][match(ids, s[[g]]$ID), ]; fl <- d$A1 != base$A1
    B[, g] <- ifelse(fl, -d$BETA, d$BETA); S[, g] <- d$SE }
  Wt <- 1 / S^2; bF <- rowSums(Wt * B) / rowSums(Wt); sF <- sqrt(1 / rowSums(Wt))
  pF <- 2 * pnorm(-abs(bF / sF)); Q <- rowSums(Wt * (B - bF)^2); df <- 2
  I2 <- pmax(0, 100 * (Q - df) / Q); Ph <- pchisq(Q, df, lower.tail = FALSE)
  C <- rowSums(Wt) - rowSums(Wt^2) / rowSums(Wt); tau2 <- pmax(0, (Q - df) / C)
  Wr <- 1 / (S^2 + tau2); bR <- rowSums(Wr * B) / rowSums(Wr); sR <- sqrt(1 / rowSums(Wr)); pR <- 2 * pnorm(-abs(bR / sR))
  meta <- data.frame(ID = ids, CHR = base$CHROM, POS = base$POS, A1 = base$A1, BETA = bF, SE = sF, P = pF,
                     Q = Q, I2 = I2, P_het = Ph, P_re = pR, stringsAsFactors = FALSE)
  ## lambda on the pruned set
  lam_rows <- data.frame(replicate = i, lambda_eur_4a = lam(e$P[e$ID %in% pruned]),
                         lambda_meta = lam(meta$P[meta$ID %in% pruned]),
                         lambda_eur_stratum = lam(s$eur$P[s$eur$ID %in% pruned]),
                         lambda_afr = lam(s$afr$P[s$afr$ID %in% pruned]), lambda_eas = lam(s$eas$P[s$eas$ID %in% pruned]))
  ## per locus
  ce <- cc(f, EUR_KEEP); neff_e <- 4 / (1 / ce[1] + 1 / ce[2])
  neff_s <- sapply(STRATA, function(g) { x <- cc(f, file.path("results/meta", g, "keep.txt")); 4 / (1 / x[1] + 1 / x[2]) })
  loc <- do.call(rbind, lapply(names(CAUSAL), function(L) {
    id <- CAUSAL[[L]]; ch <- cpos$CHR[names(CAUSAL) == L]; ps <- cpos$POS[names(CAUSAL) == L]
    rows <- list()
    for (an in c("eur_4a", "meta")) {
      d <- if (an == "eur_4a") e else meta
      chrcol <- if (an == "eur_4a") d$CHROM else d$CHR
      w <- d[chrcol == ch & abs(d$POS - ps) <= WIN & is.finite(d$BETA) & is.finite(d$SE), ]
      r <- w[w$ID == id, ]
      if (nrow(r) != 1) { rows[[an]] <- data.frame(replicate = i, locus = L, analysis = an, present = FALSE); next }
      sg <- if (r$A1 == counted[[L]]) 1 else -1
      maf_e <- s$eur$A1_FREQ[s$eur$ID == id]; maf_e <- min(maf_e, 1 - maf_e)
      fr <- sapply(STRATA, function(g) { x <- s[[g]]$A1_FREQ[s[[g]]$ID == id]; min(x, 1 - x) })
      info <- if (an == "eur_4a") neff_e * 0.25 * 2 * maf_e * (1 - maf_e) else sum(neff_s * 0.25 * 2 * fr * (1 - fr))
      pw <- pchisq(qchisq(1 - ALPHA, 1), 1, ncp = info * log(OR_GEN[[L]])^2, lower.tail = FALSE)
      fm <- abf_cs(w$BETA, w$SE, w$ID, id)
      lead_p <- min(w$P, na.rm = TRUE)
      rows[[an]] <- data.frame(replicate = i, locus = L, analysis = an, present = TRUE,
        n_eff = if (an == "eur_4a") neff_e else sum(neff_s), predicted_power = pw,
        or_hat = exp(sg * r$BETA), log_or_hat = sg * r$BETA, se = r$SE, p = r$P, significant = r$P < ALPHA,
        region_lead_p = lead_p, region_significant = lead_p < ALPHA,
        I2 = if (an == "meta") r$I2 else NA, Q = if (an == "meta") r$Q else NA,
        P_het = if (an == "meta") r$P_het else NA, p_random = if (an == "meta") r$P_re else NA,
        cs_size = fm[["cs_size"]], causal_in_cs = as.logical(fm[["causal_in_cs"]]),
        causal_rank = fm[["causal_rank"]], causal_pp = fm[["causal_pp"]], lead_is_causal = as.logical(fm[["lead_is_causal"]]))
    }
    do.call(rbind, lapply(rows, function(x) { for (cn in setdiff(names(rows[[1]]), names(x))) x[[cn]] <- NA; x[names(rows[[1]])] }))
  }))
  list(loc = loc, lam = lam_rows)
}

t0 <- Sys.time()
cat(sprintf("European scan: SEX, AGE, PC1-PC%d\n", k)); cat(sprintf("%d replicates, %d in parallel, %d threads each, %d variants (%d pruned + windows)\n",
            length(reps), N_PAR, THREADS, length(V), length(pruned)))
res <- mclapply(reps, function(f) tryCatch(one(f), error = function(e) e), mc.cores = N_PAR, mc.preschedule = FALSE)
bad <- which(sapply(res, inherits, "error"))
if (length(bad)) stop("replicates failed: ", paste(basename(reps[bad]), collapse = ", "), ": ", conditionMessage(res[[bad[1]]]))
loc <- do.call(rbind, lapply(res, `[[`, "loc")); lamt <- do.call(rbind, lapply(res, `[[`, "lam"))
write.table(loc[order(loc$replicate, loc$locus, loc$analysis), ], file.path(OUT, "pdac_demo_08_per_replicate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lamt[order(lamt$replicate), ], file.path(OUT, "pdac_demo_08_lambda.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
