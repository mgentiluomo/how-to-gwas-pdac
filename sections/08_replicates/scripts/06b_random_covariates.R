#!/usr/bin/env Rscript

################################################################################
# Section 8, step 06b: random covariates in place of the components
#
# PURPOSE
#   Step 06 shows that the mean lambda of the European scan rises with every
#   component added, under a null that holds by construction. Components
#   differ from other covariates in two ways: they are parameters, and they
#   are computed from the genotypes under test. In the Cox scan of Section 4B
#   ten random covariates reproduce the whole excess (4B step 06). This step
#   asks the same of the logistic scan: is the rise in lambda with the
#   components explained by the number of parameters alone?
#
# WHAT IT DOES, per replicate phenotype of step 01, European analysis set,
# LD-pruned European set (the set and the lambda of steps 02 and 06)
#   Four scans, logistic regression with firth-fallback as step 06:
#     SEX, AGE, PC1..PC5        SEX, AGE, R1..R5
#     SEX, AGE, PC1..PC20       SEX, AGE, R1..R20
#   R1..R20 are independent N(0,1) draws, one set of 20 per replicate (seed
#   RAND_SEED_BASE + replicate number), drawn for the individuals of the
#   European set in the order of the covariate file; the model with five uses
#   the first five, as the model with five components uses the first five
#   components. They are independent of the genotypes and of the phenotype.
#   The component models are run again here, on the same extracted variants
#   and in the same command form as the random ones, so that the paired
#   difference is not touched by the dependence of PLINK 2 on the variants
#   loaded together (step 06 scans all variants). Their lambda is compared
#   with step 06 and the largest difference is reported.
#
# NOT PART OF THE CANONICAL RUN: reads its outputs in results/, as step 02.
# SETTINGS: GWAS_08_PARALLEL (default 6) replicates at a time, GWAS_08_THREADS (4).
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_randcov_per_replicate.tsv   replicate x k: lambda with components, with random covariates
#   pdac_demo_08_randcov_summary.tsv         per k: mean, SD, SE; paired difference; lambda without components (step 06)
#   pdac_demo_08_randcov_summary.txt         plain-text report
################################################################################

suppressPackageStartupMessages(library(parallel))
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(normalizePath(file.path(script_dir, "..", "..", "..")))

OUT <- "sections/08_replicates/results"; SIM <- file.path(OUT, "sim")
IN <- "results/pca/pdac_demo_02_hwe_filt"
EUR_KEEP <- "results/pca/pdac_demo_02_eur_keep.txt"
EUR_COV <- "results/assoc/pdac_demo_04A_covar.txt"
PRUNE <- "results/pca/pdac_demo_02_prune_eur.prune.in"
STEP06 <- file.path(OUT, "pdac_demo_08_npc_per_replicate.tsv")
KS <- c(5L, 20L); K_MAX <- max(KS); RAND_SEED_BASE <- 20261200L
N_PAR <- as.integer(Sys.getenv("GWAS_08_PARALLEL", "6"))
THREADS <- as.integer(Sys.getenv("GWAS_08_THREADS", "4"))
for (f in c(paste0(IN, ".bed"), EUR_KEEP, EUR_COV, PRUNE, STEP06)) if (!file.exists(f)) stop("missing input: ", f)
reps <- sort(Sys.glob(file.path(SIM, "rep_*.pheno")))
if (!length(reps)) stop("no replicates in ", SIM, ": run step 01")
rd <- function(f, ...) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, ...)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
lam <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
tmp <- file.path(tempdir(), "randcov"); dir.create(tmp, showWarnings = FALSE)
cov_e <- rd(EUR_COV); keep_ids <- read.table(EUR_KEEP, stringsAsFactors = FALSE)[[2]]
cov_e <- cov_e[cov_e$IID %in% keep_ids, ]          # European set, order of the covariate file
if (nrow(cov_e) != length(keep_ids)) stop("covariates missing for individuals of the European set")

one <- function(f) {
  i <- as.integer(sub(".*rep_([0-9]+)\\.pheno$", "\\1", f)); tag <- sprintf("r%03d", i)
  age <- read.delim(sub("\\.pheno$", ".age", f))
  cv <- cov_e; cv$AGE <- age$AGE[match(cv$IID, age$IID)]
  set.seed(RAND_SEED_BASE + i)
  Rm <- matrix(rnorm(nrow(cv) * K_MAX), nrow(cv), K_MAX, dimnames = list(NULL, paste0("R", seq_len(K_MAX))))
  cv <- cbind(cv, Rm); names(cv)[1] <- paste0("#", names(cv)[1])
  cf <- file.path(tmp, paste0(tag, ".cov")); write.table(cv, cf, sep = "\t", quote = FALSE, row.names = FALSE)
  run <- function(pre, k) {
    out <- file.path(tmp, paste0(tag, "_", pre, k))
    st <- system2("plink2", c("--bfile", IN, "--extract", PRUNE, "--chr", "1-22", "--keep", EUR_KEEP, "--pheno", f, "--covar", cf,
                              "--covar-name", paste(c("SEX", "AGE", paste0(pre, seq_len(k))), collapse = ","),
                              "--covar-variance-standardize", "--threads", THREADS,
                              "--glm", "firth-fallback", "hide-covar", "cols=+a1freq,+beta,+orbeta,+nobs,+err",
                              "--out", out), stdout = FALSE, stderr = FALSE)
    if (st != 0) stop("plink2 failed for ", tag, " ", pre, k)
    gf <- Sys.glob(paste0(out, ".*.glm.logistic.hybrid"))[1]
    h <- sub("^#", "", strsplit(readLines(gf, n = 1), "\t")[[1]])
    g <- rd(gf, sep = "\t", colClasses = ifelse(h %in% c("TEST", "P"), NA, "NULL"))
    unlink(Sys.glob(paste0(out, "*")))
    if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
    c(lambda = lam(g$P), n = sum(is.finite(g$P) & g$P > 0))
  }
  o <- do.call(rbind, lapply(KS, function(k) { a <- run("PC", k); b <- run("R", k)
    data.frame(replicate = i, k = k, seed = RAND_SEED_BASE + i, variants_pc = a[["n"]], variants_random = b[["n"]],
               lambda_pc = a[["lambda"]], lambda_random = b[["lambda"]]) }))
  unlink(cf); o
}
t0 <- Sys.time()
cat(sprintf("%d replicates, %d in parallel, %d threads each; k %s; seed of replicate i: %d + i\n", length(reps), N_PAR, THREADS,
            paste(KS, collapse = ", "), RAND_SEED_BASE))
res <- mclapply(reps, function(f) tryCatch(one(f), error = function(e) e), mc.cores = N_PAR, mc.preschedule = FALSE)
bad <- which(sapply(res, function(x) !is.data.frame(x)))
if (length(bad)) stop("replicates failed: ", paste(basename(reps[bad]), collapse = ", "), ": ", as.character(res[[bad[1]]]))
d <- do.call(rbind, res); d <- d[order(d$replicate, d$k), ]

s6 <- read.delim(STEP06, stringsAsFactors = FALSE)
d$lambda_pc_step06 <- s6$lambda_pruned[match(paste(d$replicate, d$k), paste(s6$replicate, s6$n_pcs))]
d$lambda_0pc_step06 <- s6$lambda_pruned[match(paste(d$replicate, 0), paste(s6$replicate, s6$n_pcs))]
if (anyNA(d$lambda_pc_step06) || anyNA(d$lambda_0pc_step06)) stop("step 06 lacks a replicate or a number of components")
dmax <- max(abs(d$lambda_pc - d$lambda_pc_step06))
chk_line <- sprintf("component models here (pruned variants extracted) against step 06 (all variants scanned): max |lambda difference| %.2g", dmax)
cat(chk_line, "\n"); if (dmax > 1e-3) stop("the component models do not reproduce step 06")

ms <- function(x) c(mean = mean(x), sd = sd(x), se = sd(x) / sqrt(length(x)))
summ <- do.call(rbind, lapply(KS, function(k) { x <- d[d$k == k, ]
  a <- ms(x$lambda_pc); b <- ms(x$lambda_random); z <- ms(x$lambda_0pc_step06)
  dd <- ms(x$lambda_pc - x$lambda_random); da <- ms(x$lambda_pc - x$lambda_0pc_step06); db <- ms(x$lambda_random - x$lambda_0pc_step06)
  data.frame(k = k, replicates = nrow(x), lambda_pc_mean = a[["mean"]], lambda_pc_sd = a[["sd"]], lambda_pc_se_of_mean = a[["se"]],
             lambda_random_mean = b[["mean"]], lambda_random_sd = b[["sd"]], lambda_random_se_of_mean = b[["se"]],
             pc_minus_random_mean = dd[["mean"]], pc_minus_random_sd = dd[["sd"]], pc_minus_random_se = dd[["se"]],
             pc_minus_random_z = dd[["mean"]] / dd[["se"]], n_replicates_pc_above_random = sum(x$lambda_pc > x$lambda_random),
             lambda_0pc_mean = z[["mean"]], pc_minus_0pc_mean = da[["mean"]], pc_minus_0pc_se = da[["se"]],
             random_minus_0pc_mean = db[["mean"]], random_minus_0pc_se = db[["se"]],
             share_of_rise_reproduced_by_random = db[["mean"]] / da[["mean"]]) }))
sg <- function(x, n) { num <- sapply(x, is.numeric); x[num] <- lapply(x[num], function(v) signif(v, n)); x }
wr <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
wr(sg(d, 6), "pdac_demo_08_randcov_per_replicate.tsv"); wr(sg(summ, 4), "pdac_demo_08_randcov_summary.tsv")
R <- c("SECTION 8: RANDOM COVARIATES IN PLACE OF THE COMPONENTS, EUROPEAN SCAN",
       sprintf("%d replicates; LD-pruned European set; SEX, AGE and k components, or k independent N(0,1) covariates (seed of replicate i: %d + i)",
               length(reps), RAND_SEED_BASE), chk_line, "",
       "LAMBDA (mean, SD over replicates; SE of the mean)")
for (j in seq_len(nrow(summ))) with(summ[j, ], R <<- c(R,
  sprintf("  k %2d  components %.4f, %.4f; %.4f   random covariates %.4f, %.4f; %.4f   no components (step 06) %.4f",
          k, lambda_pc_mean, lambda_pc_sd, lambda_pc_se_of_mean, lambda_random_mean, lambda_random_sd, lambda_random_se_of_mean, lambda_0pc_mean),
  sprintf("        paired, components minus random: %+.4f (SE %.4f, z %.1f); components above random in %d of %d replicates",
          pc_minus_random_mean, pc_minus_random_se, pc_minus_random_z, n_replicates_pc_above_random, replicates),
  sprintf("        rise over no components: %+.4f (SE %.4f) with components, %+.4f (SE %.4f) with random covariates: %.0f%% reproduced",
          pc_minus_0pc_mean, pc_minus_0pc_se, random_minus_0pc_mean, random_minus_0pc_se, 100 * share_of_rise_reproduced_by_random)))
writeLines(R, file.path(OUT, "pdac_demo_08_randcov_summary.txt")); cat(R, sep = "\n")
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
