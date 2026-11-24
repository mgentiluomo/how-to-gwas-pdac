#!/usr/bin/env Rscript

################################################################################
# Section 4B: permutation null of the Cox scan with ten random covariates
#
# PURPOSE
#   04_permutation_null.R shows that the Cox scan with SEX, AGE and ten
#   within-set components has a null lambda near 1.10, against 1.02 with SEX
#   and AGE only. The ten components differ from SEX and AGE in two ways at
#   once: they are ten more parameters, and they are computed from the very
#   genotypes under test in 219 individuals, so each tested variant is partly
#   predicted by them. This step separates the two: it replaces the ten
#   components with ten covariates drawn from N(0, 1), independent of
#   everything, and repeats the permutation null. Same number of parameters,
#   no dependence on the genotypes.
#
# WHAT IT DOES
#   1. Rebuilds the survival set and the genotypes exactly as
#      04_permutation_null.R does, with the same fitter.
#   2. Draws the same N_PERM permutations as 04_permutation_null.R (same seed,
#      same call), so the three nulls (SEX and AGE; ten components; ten random
#      covariates) are paired permutation by permutation. Checked: the SEX and
#      AGE model is refitted on the first permutation and must reproduce the
#      lambda that step 04 recorded for it.
#   3. Draws N_PERM + 1 independent sets of ten N(0, 1) covariates (seed
#      GWAS_4B_RANDCOV_SEED): set 0 for the observed outcome, set k for
#      permutation k. A new set per permutation, so that the null mean is an
#      average over covariate draws and does not depend on one of them. All
#      sets are drawn before any fit.
#   4. Fits SEX, AGE and the ten random covariates on the observed outcome and
#      on every permutation, and reports lambda overall and by MAF class, the
#      permutation mean, SD and range, and the paired differences against the
#      two models of step 04.
#   5. Counts tied survival times in the survival set (the times are recorded
#      to 0.1 month; coxph handles ties with the Efron approximation).
#
# NOT PART OF THE CANONICAL RUN
#   scripts/dev/run_all.sh does not call this step. It reads the outputs of
#   the canonical run in results/ (as Section 8 does) and takes about half
#   the time of 04_permutation_null.R. scripts/dev/run_experiments.sh runs
#   it after run_all.sh, with Section 8, copies its output to
#   sections/04B_association_survival/results/ and records it in
#   EXPERIMENTS_MANIFEST.txt.
#
# SETTINGS (environment variables)
#   GWAS_4B_N_PERM, GWAS_4B_PERM_SEED   as 04_permutation_null.R (100, 20261122);
#                                       they must be those of the run whose
#                                       pdac_demo_04B_perm_lambda.tsv is read
#   GWAS_4B_RANDCOV_SEED                seed of the random covariates, default 20261125
#   GWAS_THREADS                        cores for the per-variant fits
#
# INPUT (written by the canonical run)
#   results/pca/pdac_demo_02_hwe_filt.bed/bim/fam
#   results/assoc/pdac_demo_04B_pca_within.eigenvec   (defines the survival set)
#   results/assoc/pdac_demo_04B_perm_lambda.tsv       (04_permutation_null.R)
#   demo_data/survival.txt, demo_data/covariates.txt
#
# OUTPUT (results/assoc/)
#   pdac_demo_04B_randcov_lambda.tsv    lambda by fit and MAF class
#   pdac_demo_04B_randcov_summary.tsv   the three nulls side by side, paired differences
#   pdac_demo_04B_ties.tsv              tied survival times in the survival set
#   pdac_demo_04B_randcov_summary.txt   plain-text report
################################################################################

suppressPackageStartupMessages({ library(survival); library(parallel) })

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

env_int <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (!nzchar(v)) return(default)
  x <- suppressWarnings(as.integer(v))
  if (is.na(x) || x < 1L) stop(name, " must be a positive integer, got '", v, "'")
  x
}
N_PERM   <- env_int("GWAS_4B_N_PERM", 100L)
SEED     <- env_int("GWAS_4B_PERM_SEED", 20261122L)
COV_SEED <- env_int("GWAS_4B_RANDCOV_SEED", 20261125L)
NCORE    <- env_int("GWAS_THREADS", max(1L, parallel::detectCores()))
MAF_MIN <- 0.05
N_RAND  <- 10L

IN       <- "results/pca/pdac_demo_02_hwe_filt"
EIGENVEC <- "results/assoc/pdac_demo_04B_pca_within.eigenvec"
PERM     <- "results/assoc/pdac_demo_04B_perm_lambda.tsv"
OUT_DIR  <- "results/assoc"
TMP      <- file.path(tempdir(), "cox_randcov"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
for (f in c(paste0(IN, ".bed"), EIGENVEC, PERM, "demo_data/survival.txt", "demo_data/covariates.txt"))
  if (!file.exists(f)) stop("missing input: ", f, " (run scripts/dev/run_all.sh first)")
if (Sys.which("plink2") == "") stop("plink2 not on PATH")

report <- character()
say <- function(...) { cat(..., "\n", sep = ""); flush.console() }
w   <- function(...) { line <- paste0(...); report <<- c(report, line); say(line) }
lam <- function(p) { p <- p[is.finite(p)]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
say("permutations: ", N_PERM, "   permutation seed: ", SEED, "   covariate seed: ", COV_SEED, "   cores: ", NCORE)

## ---- the survival set and the genotypes, exactly as 04_permutation_null.R ---------------
ev <- read.table(EIGENVEC, header = TRUE, comment.char = "", check.names = FALSE)
names(ev)[1:2] <- c("FID", "IID")
surv <- read.delim("demo_data/survival.txt", stringsAsFactors = FALSE)
cov  <- read.delim("demo_data/covariates.txt", stringsAsFactors = FALSE)
dat <- merge(merge(ev[, "IID", drop = FALSE], surv[, c("IID", "TIME", "EVENT")], by = "IID"),
             cov[, c("IID", "SEX", "AGE")], by = "IID")
say("survival set: ", nrow(dat), " individuals, ", sum(dat$EVENT), " events")
keep <- file.path(TMP, "keep.txt")
write.table(data.frame(dat$IID, dat$IID), keep, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
pref <- file.path(TMP, "geno")
cmd <- sprintf("plink2 --bfile %s --keep %s --chr 1-22 --maf %s --threads %d --export A-transpose --out %s",
               IN, keep, MAF_MIN, NCORE, pref)
if (system(cmd, ignore.stdout = TRUE) != 0) stop("plink2 export failed: ", cmd)
say("reading genotypes ...")
hdr <- strsplit(readLines(paste0(pref, ".traw"), n = 1L), "\t", fixed = TRUE)[[1]]
cc <- c(rep("character", 2), "NULL", "integer", rep("character", 2), rep("numeric", length(hdr) - 6))
tr <- read.table(paste0(pref, ".traw"), sep = "\t", header = TRUE, colClasses = cc, comment.char = "", check.names = FALSE)
ids <- sub("^[^_]*_", "", names(tr)[-(1:5)])
dat <- dat[match(ids, dat$IID), ]
if (anyNA(dat$IID)) stop("genotype samples do not match the survival set")
G <- as.matrix(tr[, -(1:5)]); rm(tr); invisible(gc())
unlink(paste0(pref, ".traw"))
freq <- rowMeans(G, na.rm = TRUE) / 2
maf  <- pmin(freq, 1 - freq)
say("variants: ", nrow(G))

scan <- function(time, event, covars) {
  so <- Surv(time, event)
  f <- function(i) {
    g <- G[i, ]; good <- is.finite(g)
    if (sum(good) < 20L || length(unique(g[good])) < 2L) return(NA_real_)
    d <- cbind(G = g[good], covars[good, , drop = FALSE])
    fit <- tryCatch(coxph(so[good] ~ ., data = d), error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit) || !is.finite(coef(fit)[["G"]])) return(NA_real_)
    summary(fit)$coefficients["G", "Pr(>|z|)"]
  }
  unlist(mclapply(seq_len(nrow(G)), f, mc.cores = NCORE, mc.preschedule = TRUE))
}
maf_class <- cut(maf, c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50), include.lowest = TRUE,
                 labels = c("0.05-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50"))
row_for <- function(model, fit, p) data.frame(model = model, fit = fit, n_valid_P = sum(is.finite(p)),
  lambda = round(lam(p), 4), t(round(tapply(p, maf_class, lam), 4)), check.names = FALSE, stringsAsFactors = FALSE)

## ---- permutations of step 04, and the random covariates ---------------------------------
set.seed(SEED)
perms <- lapply(seq_len(N_PERM), function(k) sample.int(nrow(dat)))
set.seed(COV_SEED)
rand <- lapply(0:N_PERM, function(k) {
  m <- matrix(rnorm(nrow(dat) * N_RAND), nrow(dat), N_RAND, dimnames = list(NULL, paste0("R", seq_len(N_RAND)))); m })
base <- as.data.frame(dat[, c("SEX", "AGE")])

t0 <- Sys.time(); mins <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
old <- read.delim(PERM, check.names = FALSE, stringsAsFactors = FALSE)
if (sum(old$model == "no_pc" & old$fit != "observed") != N_PERM)
  stop(PERM, " does not hold ", N_PERM, " permutations: settings differ from the run that wrote it")
## the permutations here must be those of step 04
p_chk <- scan(dat$TIME[perms[[1]]], dat$EVENT[perms[[1]]], base)
chk <- c(here = round(lam(p_chk), 4), step04 = old$lambda[old$model == "no_pc" & old$fit == "perm001"])
if (!isTRUE(all.equal(chk[["here"]], chk[["step04"]])))
  stop("first permutation, SEX and AGE: lambda ", chk[["here"]], " here against ", chk[["step04"]], " in ", PERM)
say(sprintf("check: SEX and AGE on permutation 1 gives lambda %.4f, as step 04 (%.1f min)", chk[["here"]], mins()))

rows <- list()
p <- scan(dat$TIME, dat$EVENT, cbind(base, rand[[1]]))
rows[[1]] <- row_for("random10", "observed", p)
say(sprintf("observed  random10  lambda %.4f  (%.1f min)", rows[[1]]$lambda, mins()))
for (k in seq_len(N_PERM)) {
  o <- perms[[k]]
  p <- scan(dat$TIME[o], dat$EVENT[o], cbind(base, rand[[k + 1L]]))
  rows[[length(rows) + 1L]] <- row_for("random10", sprintf("perm%03d", k), p)
  say(sprintf("permutation %3d/%d  random10 %.4f  (%.1f min elapsed)", k, N_PERM, rows[[length(rows)]]$lambda, mins()))
}
tab <- do.call(rbind, rows)

## ---- summaries ---------------------------------------------------------------------------
pl <- function(d, m) { x <- d[d$model == m & d$fit != "observed", ]; setNames(x$lambda, x$fit) }
nulls <- list(no_pc = pl(old, "no_pc"), random10 = pl(tab, "random10"), primary = pl(old, "primary"))
fits <- names(nulls$random10)
if (!all(sapply(nulls, function(x) identical(names(x), fits)))) stop("permutation labels differ between the tables")
desc <- c(no_pc = "SEX, AGE", random10 = "SEX, AGE, ten N(0,1) covariates", primary = "SEX, AGE, within-set PC1-PC10")
npar <- c(no_pc = 3L, random10 = 13L, primary = 13L)
obs <- c(no_pc = old$lambda[old$model == "no_pc" & old$fit == "observed"],
         random10 = tab$lambda[tab$fit == "observed"],
         primary = old$lambda[old$model == "primary" & old$fit == "observed"])
n_ev <- sum(dat$EVENT)
summ <- do.call(rbind, lapply(names(nulls), function(m) { x <- nulls[[m]]
  data.frame(model = m, covariates = desc[[m]], parameters_with_variant = npar[[m]],
             events_per_parameter = round(n_ev / npar[[m]], 1), n_perm = N_PERM,
             source = if (m == "random10") "this step" else "04_permutation_null.R",
             lambda_observed = obs[[m]], perm_mean = round(mean(x), 4), perm_sd = round(sd(x), 4),
             perm_se_of_mean = round(sd(x) / sqrt(length(x)), 4), perm_min = min(x), perm_max = max(x),
             stringsAsFactors = FALSE) }))
pd <- function(a, b) { d <- nulls[[a]] - nulls[[b]]
  data.frame(contrast = paste(a, "minus", b), mean_diff = round(mean(d), 4), sd_diff = round(sd(d), 4),
             se_diff = round(sd(d) / sqrt(length(d)), 4), n_perm_positive = sum(d > 0), stringsAsFactors = FALSE) }
paired <- rbind(pd("random10", "no_pc"), pd("primary", "random10"), pd("primary", "no_pc"))
excess_total <- mean(nulls$primary) - mean(nulls$no_pc)
paired$share_of_primary_minus_no_pc <- round(paired$mean_diff / excess_total, 3)

## ---- tied times ----------------------------------------------------------------------------
et <- dat$TIME[dat$EVENT == 1]; tt <- table(et)
ties <- data.frame(quantity = c("individuals", "events", "distinct_times_all", "distinct_event_times",
                                "events_at_a_tied_event_time", "event_times_with_ties", "largest_number_of_events_at_one_time",
                                "censored_at_an_event_time", "time_resolution_months"),
                   value = c(nrow(dat), n_ev, length(unique(dat$TIME)), length(tt), sum(tt[tt > 1]), sum(tt > 1), max(tt),
                             sum(dat$EVENT == 0 & dat$TIME %in% et), signif(min(diff(sort(unique(dat$TIME)))), 4)))

out <- function(x, name) write.table(x, file.path(OUT_DIR, name), sep = "\t", quote = FALSE, row.names = FALSE)
out(tab, "pdac_demo_04B_randcov_lambda.tsv")
summ_out <- merge(summ, data.frame(model = "random10", randcov_seed = COV_SEED, perm_seed = SEED), all.x = TRUE)
summ_out <- summ_out[match(names(nulls), summ_out$model), ]
out(summ_out, "pdac_demo_04B_randcov_summary.tsv")
out(paired, "pdac_demo_04B_randcov_paired.tsv")
out(ties, "pdac_demo_04B_ties.tsv")

w("SECTION 4B: PERMUTATION NULL OF THE COX SCAN WITH TEN RANDOM COVARIATES")
w("  permutations: ", N_PERM, " (seed ", SEED, ", those of 04_permutation_null.R)   random covariates: seed ", COV_SEED,
  ", a new set of ten for every fit")
w("  survival set: ", nrow(dat), " individuals, ", n_ev, " events; ", nrow(G), " variants (MAF >= ", MAF_MIN, " within the set)")
w("")
w("  null lambda over the same permutations (mean, SD; SE of the mean)")
for (i in seq_len(nrow(summ))) with(summ[i, ], w(sprintf("  %-9s %-32s %2d parameters, %4.1f events per parameter: %.4f, %.4f; %.4f   observed %.4f",
  model, covariates, parameters_with_variant, events_per_parameter, perm_mean, perm_sd, perm_se_of_mean, lambda_observed)))
w("")
w("  paired differences, permutation by permutation (mean, SD; SE; share of primary minus no_pc)")
for (i in seq_len(nrow(paired))) with(paired[i, ], w(sprintf("  %-26s %+.4f, %.4f; %.4f; %.2f   positive in %d of %d",
  contrast, mean_diff, sd_diff, se_diff, share_of_primary_minus_no_pc, n_perm_positive, N_PERM)))
w("")
w("  tied survival times (recorded to ", ties$value[ties$quantity == "time_resolution_months"], " month): ",
  ties$value[ties$quantity == "events_at_a_tied_event_time"], " of ", n_ev, " events share their time with another event; ",
  ties$value[ties$quantity == "distinct_event_times"], " distinct event times; at most ",
  ties$value[ties$quantity == "largest_number_of_events_at_one_time"], " events at one time")
w("")
w("HOW TO READ THIS")
w("  random10 minus no_pc is what ten parameters cost when they carry no")
w("  information about the genotypes. primary minus random10 is what is added")
w("  by the components being computed from the tested genotypes.")
w("")
w("  elapsed: ", sprintf("%.1f min", mins()), " on ", NCORE, " cores")
writeLines(report, file.path(OUT_DIR, "pdac_demo_04B_randcov_summary.txt"))
unlink(TMP, recursive = TRUE)
