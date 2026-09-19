#!/usr/bin/env Rscript

################################################################################
# Section 4B, diagnostic: is the inflation in the survival scan real?
#
# THE QUESTION
#   The Cox scan returns lambda = 1.12 on 238,868 variants, where the binary
#   scan on the same individuals returns 1.02. With that many tests the median
#   estimator has a standard error near 0.004, so 1.12 is not sampling noise
#   and needs an explanation before anything from the scan is reported.
#
#   Two explanations are worth separating.
#
#   (a) The Wald statistic is not calibrated at this event count. The scan uses
#       the per-coefficient Wald test, which is anti-conservative when events
#       per parameter are low and becomes more so as the coefficient grows.
#       With 166 events and 13 parameters there are about 12.8 events per
#       parameter, at the conventional floor. If this is the cause, the
#       likelihood-ratio and score tests, which do not share the failure mode,
#       will return lambda near 1.00 on the same variants.
#
#   (b) Residual population structure among the cases. The principal components
#       in the model were computed across all 610 Europeans, cases and controls
#       together. The survival analysis uses the 217 cases only, and structure
#       specific to that subset need not be captured by components estimated in
#       the larger set. If this is the cause, all three test statistics will be
#       inflated together, because the problem is in the model rather than in
#       the statistic.
#
#   The two have different fixes: (a) means reporting the scan with a different
#   statistic, (b) means recomputing the components within the cases.
#
# METHOD
#   Draw a random subset of variants, refit each one, and compute all three
#   statistics for the genotype term. The likelihood-ratio and score tests are
#   computed against a null model fitted on the same individuals as the full
#   model, so that missing genotypes cannot make the comparison unfair.
#
# USAGE
#   Rscript 02_test_calibration.R [n_variants] [n_pcs]
#   defaults: 5000 variants, 10 principal components
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("survival", quietly = TRUE))
    stop("the 'survival' package is required")
  library(survival); library(parallel)
})

args_all  <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg))
  normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
root <- script_dir
for (up in list(c("..", ".."), c("..", "..", ".."))) {
  cand <- normalizePath(do.call(file.path, c(list(script_dir), as.list(up))),
                        mustWork = FALSE)
  if (dir.exists(file.path(cand, "scripts", "dev"))) { root <- cand; break }
}
if (!dir.exists(file.path(root, "scripts", "dev"))) root <- getwd()
setwd(root)

args  <- commandArgs(trailingOnly = TRUE)
N_VAR <- if (length(args) >= 1) as.integer(args[1]) else 5000L
N_PC  <- if (length(args) >= 2) as.integer(args[2]) else 10L

IN       <- "results/pca/pdac_demo_02_hwe_filt"
KEEP_EUR <- "results/pca/pdac_demo_02_eur_keep.txt"
SURV     <- "demo_data/survival.txt"
COVAR    <- "results/assoc/pdac_demo_04A_covar.txt"
SCAN     <- "results/assoc/pdac_demo_04B_cox.tsv"
OUT_DIR  <- "results/assoc"
TMP      <- file.path(tempdir(), "cox_calib")
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

for (f in c(paste0(IN, ".bed"), KEEP_EUR, SURV, COVAR, SCAN))
  if (!file.exists(f)) stop("required input missing: ", f)
if (Sys.which("plink2") == "") stop("plink2 not on PATH")

first_env_int <- function(vars) {
  for (v in vars) {
    x <- suppressWarnings(as.integer(Sys.getenv(v, "")))
    if (!is.na(x) && x > 0) return(x)
  }
  NULL
}
NCORE <- first_env_int(c("GWAS_THREADS", "SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS"))
if (is.null(NCORE)) NCORE <- max(1L, suppressWarnings(detectCores(logical = TRUE)) - 1L)
if (.Platform$OS.type != "unix") NCORE <- 1L

lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1, lower.tail = FALSE)
}

## --- phenotype and covariates, exactly as the scan built them ---------------
surv <- read.delim(SURV, stringsAsFactors = FALSE, check.names = FALSE)
names(surv) <- sub("^#", "", names(surv))
surv <- data.frame(IID = as.character(surv$IID),
                   TIME = as.numeric(surv$TIME),
                   EVENT = as.integer(surv$EVENT), stringsAsFactors = FALSE)
surv <- surv[is.finite(surv$TIME) & surv$TIME > 0 & surv$EVENT %in% c(0L, 1L), ]

cov <- read.delim(COVAR, stringsAsFactors = FALSE, check.names = FALSE)
names(cov)[1] <- sub("^#", "", names(cov)[1])
cov$IID <- as.character(cov$IID)
cov_cols <- c("SEX", "AGE", paste0("PC", seq_len(N_PC)))

keep_ids <- as.character(read.table(KEEP_EUR, stringsAsFactors = FALSE)[[2]])
dat <- merge(surv, cov[, c("IID", cov_cols)], by = "IID")
dat <- dat[dat$IID %in% keep_ids, ]
dat <- dat[complete.cases(dat), ]

cat("individuals: ", nrow(dat), "   events: ", sum(dat$EVENT), "\n", sep = "")
cat("parameters in the model: ", length(cov_cols) + 1L,
    "   events per parameter: ",
    sprintf("%.1f", sum(dat$EVENT) / (length(cov_cols) + 1L)), "\n\n", sep = "")

## --- a random subset of the variants the scan actually tested ---------------
scan <- read.delim(SCAN, stringsAsFactors = FALSE)
scan <- scan[is.finite(scan$P), ]
set.seed(2026)
sel <- scan[sample(nrow(scan), min(N_VAR, nrow(scan))), ]
snplist <- file.path(TMP, "sel.snplist")
writeLines(sel$ID, snplist)

keep_file <- file.path(TMP, "keep.txt")
write.table(data.frame(dat$IID, dat$IID), keep_file, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("exporting ", nrow(sel), " variants ...\n", sep = "")
cmd <- sprintf("plink2 --bfile %s --keep %s --extract %s --threads %d --export A-transpose --out %s",
               shQuote(IN), shQuote(keep_file), shQuote(snplist),
               NCORE, shQuote(file.path(TMP, "g")))
if (system(cmd, ignore.stdout = TRUE) != 0) stop("plink2 export failed:\n  ", cmd)

traw <- file.path(TMP, "g.traw")
blk  <- read.table(traw, sep = "\t", header = TRUE, comment.char = "",
                   check.names = FALSE, stringsAsFactors = FALSE)
meta_n <- 6L
ids  <- blk[[2]]
G    <- as.matrix(blk[, -seq_len(meta_n), drop = FALSE])
samp <- sub("^[^_]*_", "", colnames(blk)[-seq_len(meta_n)])
idx  <- match(samp, dat$IID)
if (anyNA(idx)) idx <- match(colnames(blk)[-seq_len(meta_n)], dat$IID)
if (anyNA(idx)) stop("sample mismatch")
dat <- dat[idx, ]

surv_obj <- Surv(dat$TIME, dat$EVENT)
base_df  <- as.data.frame(dat[, cov_cols])

## --- three statistics per variant -------------------------------------------
## The null model is refitted on the same individuals as the full model, so a
## variant with missing genotypes does not get a null fitted on more people
## than its alternative.
one <- function(i) {
  g <- G[i, ]
  good <- is.finite(g)
  if (sum(good) < 20L || length(unique(g[good])) < 2L) return(c(NA, NA, NA))
  d0 <- base_df[good, , drop = FALSE]
  d1 <- cbind(G = g[good], d0)
  f1 <- tryCatch(coxph(surv_obj[good] ~ ., data = d1),
                 error = function(e) NULL, warning = function(w) NULL)
  f0 <- tryCatch(coxph(surv_obj[good] ~ ., data = d0),
                 error = function(e) NULL, warning = function(w) NULL)
  if (is.null(f1) || is.null(f0) || !is.finite(coef(f1)[["G"]])) return(c(NA, NA, NA))
  z    <- summary(f1)$coefficients["G", "z"]
  wald <- pchisq(z^2, 1, lower.tail = FALSE)
  lr   <- pchisq(2 * (f1$loglik[2] - f0$loglik[2]), 1, lower.tail = FALSE)
  # score test for the genotype term: the full model's score statistic
  # evaluated at the null fit, obtained by starting the fit at zero
  fs <- tryCatch(coxph(surv_obj[good] ~ ., data = d1,
                       init = c(0, coef(f0)), control = coxph.control(iter.max = 0)),
                 error = function(e) NULL)
  sc <- if (is.null(fs)) NA_real_ else
    pchisq(fs$score - f0$score, 1, lower.tail = FALSE)
  c(wald, lr, sc)
}

cat("refitting on ", NCORE, " core(s) ...\n", sep = "")
res <- if (NCORE > 1L)
  mclapply(seq_len(nrow(G)), one, mc.cores = NCORE) else lapply(seq_len(nrow(G)), one)
M <- do.call(rbind, lapply(res, function(x)
  if (is.numeric(x) && length(x) == 3L) x else c(NA, NA, NA)))
colnames(M) <- c("P_wald", "P_lr", "P_score")

out <- data.frame(ID = ids, M, stringsAsFactors = FALSE)
write.table(out, file.path(OUT_DIR, "pdac_demo_04B_calibration.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

l_w <- lam(out$P_wald); l_l <- lam(out$P_lr); l_s <- lam(out$P_score)

cat("\n")
cat("GENOMIC INFLATION BY TEST STATISTIC, ", sum(is.finite(out$P_wald)),
    " variants\n", sep = "")
cat(sprintf("  Wald              lambda = %.4f\n", l_w))
cat(sprintf("  likelihood ratio  lambda = %.4f\n", l_l))
cat(sprintf("  score             lambda = %.4f\n", l_s))
cat("\n")
cat("HOW TO READ THIS\n")
if (is.finite(l_l) && l_l < 1.05 && l_w > 1.08) {
  cat("  The likelihood-ratio test is close to calibrated where the Wald test is\n")
  cat("  not. The inflation is a property of the statistic at this event count,\n")
  cat("  not of the sample. Report the scan on the likelihood-ratio test, and\n")
  cat("  state the event count and the events-per-parameter ratio alongside it.\n")
} else if (is.finite(l_l) && l_l > 1.08) {
  cat("  All three statistics are inflated together, so the problem is in the\n")
  cat("  model rather than in the test. The most likely cause is residual\n")
  cat("  structure among the cases: the components were estimated across cases\n")
  cat("  and controls, and this analysis uses cases only. Recompute the\n")
  cat("  components within the survival set and re-run before reporting.\n")
} else {
  cat("  The result is between the two clean cases. Inspect the QQ plot written\n")
  cat("  below: uniform departure across the range points to structure, whereas\n")
  cat("  departure concentrated in the tail points to the statistic.\n")
}

png(file.path(OUT_DIR, "pdac_demo_04B_calibration_qq.png"),
    width = 1800, height = 700, res = 170)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (nm in c("P_wald", "P_lr", "P_score")) {
  p <- out[[nm]]; p <- sort(p[is.finite(p) & p > 0 & p <= 1])
  if (!length(p)) { plot.new(); next }
  plot(-log10(ppoints(length(p))), -log10(p), pch = 20, cex = 0.4, col = "#0072B2",
       xlab = expression(Expected~-log[10](P)), ylab = expression(Observed~-log[10](P)),
       main = sprintf("%s  (lambda = %.3f)", sub("^P_", "", nm), lam(out[[nm]])))
  abline(0, 1, col = "grey40")
}
par(op); invisible(dev.off())

cat("\nWritten:\n")
cat("   ", file.path(OUT_DIR, "pdac_demo_04B_calibration.tsv"), "\n")
cat("   ", file.path(OUT_DIR, "pdac_demo_04B_calibration_qq.png"), "\n")
unlink(TMP, recursive = TRUE)