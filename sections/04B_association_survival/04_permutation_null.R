#!/usr/bin/env Rscript

################################################################################
# Section 4B: permutation null of the Cox scan
#
# PURPOSE
#   Measure the genomic inflation factor that the Section 4B Cox scan shows
#   when the null holds by construction, so that the observed lambda can be
#   read against it rather than against 1.00.
#
# WHY
#   At this event count the Cox tests are not calibrated: lambda is above 1
#   even though the simulation gives no variant, and no covariate, an effect
#   on survival (demo_dataset/scripts/04_simulate_phenotype.R). Whether an
#   observed lambda reflects confounding, signal or the test itself cannot be
#   decided from the observed lambda alone. Permuting the outcome can decide
#   it: confounding and signal do not survive a permutation; miscalibration
#   of the test does.
#
# WHAT IT DOES
#   1. Rebuilds the survival set and the genotypes exactly as 01_cox_survival.R
#      does (same individuals, same within-set components, same MAF floor,
#      same fitter, which drops fits that warn), refits the primary model and
#      checks that its P values match the canonical scan.
#   2. Fits two models on the observed outcome:
#        primary   SEX, AGE, within-set PC1-PC10 (the model of 01_cox_survival.R)
#        no_pc     SEX, AGE only
#   3. Draws N_PERM permutations of the individuals and, for each, assigns the
#      (TIME, EVENT) pairs jointly to the permuted order. Genotypes, SEX, AGE
#      and the components stay together, so any link between outcome and
#      genotype or structure is destroyed while every other feature of the
#      data is kept. The same permutations are applied to both models, so
#      the two null distributions are paired.
#   4. Reports lambda overall and by minor-allele-frequency class for every
#      fit, and for each model the permutation mean, standard deviation and
#      range, and the empirical P of the observed lambda,
#      (number of permutations at or above it + 1) / (N_PERM + 1).
#
# SETTINGS (environment variables; run_all.sh sets and records them)
#   GWAS_4B_N_PERM      number of permutations, default 100
#   GWAS_4B_PERM_SEED   seed of the permutations, default 20261122
#   A run with values other than the defaults is not a canonical run;
#   run_all.sh says so in its log and in RUN_MANIFEST.txt.
#   GWAS_THREADS        cores for the per-variant fits (default: all)
#
# COST
#   One Cox fit per variant per model per permutation: about 2 x (N_PERM + 1)
#   genome-wide scans of about 240,000 variants. Expect hours, not minutes,
#   at the default N_PERM; elapsed time is printed after every permutation.
#
# INPUT (written earlier in the canonical run)
#   results/pca/pdac_demo_02_hwe_filt.bed/bim/fam
#   results/assoc/pdac_demo_04B_pca_within.eigenvec   (01_cox_survival.R)
#   results/assoc/pdac_demo_04B_cox.tsv               (01_cox_survival.R)
#   demo_data/survival.txt, demo_data/covariates.txt
#
# OUTPUT (results/assoc/, moved to sections/04B_association_survival/results/)
#   pdac_demo_04B_perm_lambda.tsv       lambda by fit and MAF class
#   pdac_demo_04B_perm_summary.tsv      per model: observed, permutation
#                                       mean/sd/min/max, n at or above, P
#   pdac_demo_04B_perm_validation.tsv   refitted primary scan vs canonical
#   pdac_demo_04B_perm_variants_by_maf.tsv
#   pdac_demo_04B_perm_summary.txt      plain-text report, with N_PERM and seed
#
# Derived from the diagnostics of 22 November 2026 (inflation_4B.R and
# perm_nopc_4B.R, 20 permutations each, outside the repository).
################################################################################

suppressPackageStartupMessages({ library(survival); library(parallel) })

## ---- project root, as in the other 4B scripts --------------------------------
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

## ---- settings ------------------------------------------------------------------
env_int <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (!nzchar(v)) return(default)
  x <- suppressWarnings(as.integer(v))
  if (is.na(x) || x < 1L) stop(name, " must be a positive integer, got '", v, "'")
  x
}
N_PERM_DEFAULT <- 100L
SEED_DEFAULT   <- 20261122L
N_PERM <- env_int("GWAS_4B_N_PERM", N_PERM_DEFAULT)
SEED   <- env_int("GWAS_4B_PERM_SEED", SEED_DEFAULT)
NCORE  <- env_int("GWAS_THREADS", max(1L, parallel::detectCores()))
MAF_MIN <- 0.05
N_PC <- 10L
CANONICAL <- N_PERM == N_PERM_DEFAULT && SEED == SEED_DEFAULT

IN       <- "results/pca/pdac_demo_02_hwe_filt"
EIGENVEC <- "results/assoc/pdac_demo_04B_pca_within.eigenvec"
COX      <- "results/assoc/pdac_demo_04B_cox.tsv"
OUT_DIR  <- "results/assoc"
TMP      <- file.path(tempdir(), "cox_perm"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
for (f in c(paste0(IN, ".bed"), EIGENVEC, COX, "demo_data/survival.txt", "demo_data/covariates.txt"))
  if (!file.exists(f)) stop("missing input: ", f, " (run 01_cox_survival.R first)")
if (Sys.which("plink2") == "") stop("plink2 not on PATH")

report <- character()
say <- function(...) { cat(..., "\n", sep = ""); flush.console() }
w   <- function(...) { line <- paste0(...); report <<- c(report, line); say(line) }
lam <- function(p) { p <- p[is.finite(p)]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }

say("permutations: ", N_PERM, "   seed: ", SEED, "   cores: ", NCORE,
    if (CANONICAL) "" else "   (NOT the canonical settings)")

## ---- the survival set, exactly as 01_cox_survival.R defined it --------------------
ev <- read.table(EIGENVEC, header = TRUE, comment.char = "", check.names = FALSE)
names(ev)[1:2] <- c("FID", "IID")
pc_names <- paste0("PC", seq_len(N_PC))
surv <- read.delim("demo_data/survival.txt", stringsAsFactors = FALSE)
cov  <- read.delim("demo_data/covariates.txt", stringsAsFactors = FALSE)
dat <- merge(merge(ev[, c("IID", pc_names)], surv[, c("IID", "TIME", "EVENT")], by = "IID"),
             cov[, c("IID", "SEX", "AGE")], by = "IID")
say("survival set: ", nrow(dat), " individuals, ", sum(dat$EVENT), " events")
keep <- file.path(TMP, "keep.txt")
write.table(data.frame(dat$IID, dat$IID), keep, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

## ---- genotypes -----------------------------------------------------------------------
pref <- file.path(TMP, "geno")
cmd <- sprintf("plink2 --bfile %s --keep %s --chr 1-22 --maf %s --threads %d --export A-transpose --out %s",
               IN, keep, MAF_MIN, NCORE, pref)
if (system(cmd, ignore.stdout = TRUE) != 0) stop("plink2 export failed: ", cmd)
say("reading genotypes ...")
hdr <- strsplit(readLines(paste0(pref, ".traw"), n = 1L), "\t", fixed = TRUE)[[1]]
cc <- c(rep("character", 2), "NULL", "integer", rep("character", 2), rep("numeric", length(hdr) - 6))
tr <- read.table(paste0(pref, ".traw"), sep = "\t", header = TRUE, colClasses = cc,
                 comment.char = "", check.names = FALSE)
ids <- sub("^[^_]*_", "", names(tr)[-(1:5)])
dat <- dat[match(ids, dat$IID), ]
if (anyNA(dat$IID)) stop("genotype samples do not match the survival set")
G <- as.matrix(tr[, -(1:5)]); vid <- tr[[2]]; rm(tr); invisible(gc())
unlink(paste0(pref, ".traw"))
freq <- rowMeans(G, na.rm = TRUE) / 2
maf  <- pmin(freq, 1 - freq)
say("variants: ", nrow(G))

## ---- fitter identical to 01_cox_survival.R ---------------------------------------------
scan <- function(time, event, covars) {
  so <- Surv(time, event)
  f <- function(i) {
    g <- G[i, ]; good <- is.finite(g)
    if (sum(good) < 20L || length(unique(g[good])) < 2L) return(NA_real_)
    d <- cbind(G = g[good], covars[good, , drop = FALSE])
    fit <- tryCatch(coxph(so[good] ~ ., data = d), error = function(e) NULL,
                    warning = function(w) NULL)
    if (is.null(fit) || !is.finite(coef(fit)[["G"]])) return(NA_real_)
    summary(fit)$coefficients["G", "Pr(>|z|)"]
  }
  unlist(mclapply(seq_len(nrow(G)), f, mc.cores = NCORE, mc.preschedule = TRUE))
}

models <- list(
  primary = as.data.frame(dat[, c("SEX", "AGE", pc_names)]),
  no_pc   = as.data.frame(dat[, c("SEX", "AGE")])
)
maf_class <- cut(maf, c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50), include.lowest = TRUE,
                 labels = c("0.05-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50"))
row_for <- function(model, fit, p) {
  data.frame(model = model, fit = fit, n_valid_P = sum(is.finite(p)),
             lambda = round(lam(p), 4),
             t(round(tapply(p, maf_class, lam), 4)), check.names = FALSE,
             stringsAsFactors = FALSE)
}

## ---- observed outcome ------------------------------------------------------------------
t0 <- Sys.time()
mins <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
rows <- list()
p_obs <- list()
for (m in names(models)) {
  p_obs[[m]] <- scan(dat$TIME, dat$EVENT, models[[m]])
  rows[[length(rows) + 1L]] <- row_for(m, "observed", p_obs[[m]])
  say(sprintf("observed  %-7s  lambda %.4f  (%.1f min)", m, lam(p_obs[[m]]), mins()))
}

## The refitted primary scan must reproduce the canonical one.
can <- read.delim(COX, stringsAsFactors = FALSE)
mm <- match(vid, can$ID)
valid <- data.frame(n_variants_here = length(vid), n_matched_canonical = sum(!is.na(mm)),
                    max_abs_diff_P = signif(max(abs(p_obs$primary - can$P[mm]), na.rm = TRUE), 3),
                    lambda_here = round(lam(p_obs$primary), 4),
                    lambda_canonical = round(lam(can$P), 4))
if (valid$n_matched_canonical != valid$n_variants_here || valid$max_abs_diff_P > 1e-8)
  stop("the refitted primary scan does not reproduce 01_cox_survival.R: ",
       paste(names(valid), unlist(valid), collapse = ", "))
say("validation: ", valid$n_matched_canonical, " of ", valid$n_variants_here,
    " variants matched, max |P difference| ", valid$max_abs_diff_P)

## ---- permutations ----------------------------------------------------------------------
## All permutations are drawn before any fit, so that each one depends only on
## the seed and its index, not on how many fits ran before it.
set.seed(SEED)
perms <- lapply(seq_len(N_PERM), function(k) sample.int(nrow(dat)))
for (k in seq_len(N_PERM)) {
  o <- perms[[k]]
  for (m in names(models)) {
    p <- scan(dat$TIME[o], dat$EVENT[o], models[[m]])
    rows[[length(rows) + 1L]] <- row_for(m, sprintf("perm%03d", k), p)
  }
  last <- tail(rows, 2L)
  say(sprintf("permutation %3d/%d  primary %.4f  no_pc %.4f  (%.1f min elapsed)",
              k, N_PERM, last[[1]]$lambda, last[[2]]$lambda, mins()))
}
tab <- do.call(rbind, rows)

## ---- summaries --------------------------------------------------------------------------
summ <- do.call(rbind, lapply(names(models), function(m) {
  obs  <- tab$lambda[tab$model == m & tab$fit == "observed"]
  perm <- tab$lambda[tab$model == m & tab$fit != "observed"]
  n_ge <- sum(perm >= obs)
  data.frame(model = m, n_perm = N_PERM, seed = SEED,
             lambda_observed = obs,
             perm_mean = round(mean(perm), 4), perm_sd = round(sd(perm), 4),
             perm_min = min(perm), perm_max = max(perm),
             n_perm_at_or_above_observed = n_ge,
             empirical_P = round((n_ge + 1) / (N_PERM + 1), 4),
             stringsAsFactors = FALSE)
}))
nvar <- data.frame(maf_class = levels(maf_class), n_variants = as.integer(table(maf_class)))

out <- function(x, name) write.table(x, file.path(OUT_DIR, name), sep = "\t",
                                     quote = FALSE, row.names = FALSE)
out(tab,   "pdac_demo_04B_perm_lambda.tsv")
out(summ,  "pdac_demo_04B_perm_summary.tsv")
out(valid, "pdac_demo_04B_perm_validation.tsv")
out(nvar,  "pdac_demo_04B_perm_variants_by_maf.tsv")

w("SECTION 4B: PERMUTATION NULL OF THE COX SCAN")
w("  permutations: ", N_PERM, "   seed: ", SEED,
  if (CANONICAL) "   (canonical settings)" else "   (NOT the canonical settings)")
w("  survival set: ", nrow(dat), " individuals, ", sum(dat$EVENT), " events; ",
  nrow(G), " variants (MAF >= ", MAF_MIN, " within the set)")
w("  validation: refitted primary scan matches 01_cox_survival.R on ",
  valid$n_matched_canonical, " of ", valid$n_variants_here,
  " variants, max |P difference| ", valid$max_abs_diff_P)
w("")
for (i in seq_len(nrow(summ))) {
  s <- summ[i, ]
  w(sprintf("  %-7s observed lambda %.4f | permutation null: mean %.4f, sd %.4f, range %.4f to %.4f",
            s$model, s$lambda_observed, s$perm_mean, s$perm_sd, s$perm_min, s$perm_max))
  w(sprintf("          %d of %d permutations at or above the observed value; empirical P %.4f",
            s$n_perm_at_or_above_observed, N_PERM, s$empirical_P))
}
w("")
w("HOW TO READ THIS")
w("  The permuted outcome carries no confounding and no signal, so the")
w("  permutation mean is the lambda this test gives under the null with this")
w("  covariate set. Read the observed lambda against it, not against 1.00.")
w("  The standard deviation across permutations is the noise of a single")
w("  lambda from these correlated variants.")
w("")
w("  elapsed: ", sprintf("%.1f min", mins()), " on ", NCORE, " cores")
writeLines(report, file.path(OUT_DIR, "pdac_demo_04B_perm_summary.txt"))
unlink(TMP, recursive = TRUE)
