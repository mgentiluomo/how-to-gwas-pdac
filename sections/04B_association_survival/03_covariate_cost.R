#!/usr/bin/env Rscript

################################################################################
# Section 4B, diagnostic: what does the covariate set cost in power?
#
# THE GAP THIS FILLS
#   The power calculation in Section 7 computes the detectable hazard ratio for
#   a marginal test: genotype against outcome, no covariates. The scan actually
#   fits thirteen parameters. Those are not the same quantity, and the
#   difference is not negligible here, because 166 events against 13 parameters
#   is 12.8 events per parameter, at the conventional floor.
#
#   Adjusting for covariates costs precision in two ways. Any covariate
#   correlated with genotype removes genotype variation available to estimate
#   the genetic effect. And independently of correlation, each parameter
#   consumes information, which matters when events are few.
#
#   Rather than assume a correction factor, this script measures the cost
#   directly: it refits a sample of variants with the full covariate set, with
#   a reduced set, and with none, and reports the ratio of standard errors.
#   That ratio translates into a power penalty, because the non-centrality
#   parameter scales as (beta/SE)^2.
#
# WHAT IT REPORTS
#   For each covariate set: the median standard error of the genotype
#   coefficient, the inflation relative to the unadjusted fit, the implied
#   multiplier on the detectable hazard ratio, and the genomic inflation factor
#   so that the precision cost can be read against the calibration.
#
# USAGE
#   Rscript 03_covariate_cost.R [n_variants]
#   default 5000 variants; resources detected as in the other scripts
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

IN       <- "results/pca/pdac_demo_02_hwe_filt"
KEEP_EUR <- "results/pca/pdac_demo_02_eur_keep.txt"
SURV     <- "demo_data/survival.txt"
COVAR_IN <- "demo_data/covariates.txt"
PCA_IN   <- "results/assoc/pdac_demo_04B_pca_within.eigenvec"
SCAN     <- "results/assoc/pdac_demo_04B_cox.tsv"
OUT_DIR  <- "results/assoc"
TMP      <- file.path(tempdir(), "cov_cost")
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

for (f in c(paste0(IN, ".bed"), KEEP_EUR, SURV, COVAR_IN, PCA_IN, SCAN))
  if (!file.exists(f)) stop("required input missing: ", f,
                            "\n  run 01_cox_survival.R first")
if (Sys.which("plink2") == "") stop("plink2 not on PATH")

env_int <- function(v) { x <- suppressWarnings(as.integer(Sys.getenv(v, ""))); if (is.na(x)) NULL else x }
NCORE <- env_int("GWAS_THREADS")
if (is.null(NCORE)) NCORE <- max(1L, suppressWarnings(detectCores(logical = TRUE)) - 1L)
if (.Platform$OS.type != "unix") NCORE <- 1L

lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1, lower.tail = FALSE)
}

## --- assemble the survival set, exactly as the scan does --------------------
sv <- read.delim(SURV, stringsAsFactors = FALSE, check.names = FALSE)
names(sv) <- sub("^#", "", names(sv))
sv <- data.frame(IID = as.character(sv$IID), TIME = as.numeric(sv$TIME),
                 EVENT = as.integer(sv$EVENT), stringsAsFactors = FALSE)
sv <- sv[is.finite(sv$TIME) & sv$TIME > 0 & sv$EVENT %in% c(0L, 1L), ]

cv <- read.delim(COVAR_IN, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
names(cv)[1] <- sub("^#", "", names(cv)[1]); cv$IID <- as.character(cv$IID)

pc <- read.delim(PCA_IN, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
names(pc)[1] <- sub("^#", "", names(pc)[1]); pc$IID <- as.character(pc$IID)
pc_names <- grep("^PC[0-9]+$", names(pc), value = TRUE)

keep <- as.character(read.table(KEEP_EUR, stringsAsFactors = FALSE)[[2]])
dat <- merge(merge(sv, cv[, c("IID", "SEX", "AGE")], by = "IID"),
             pc[, c("IID", pc_names)], by = "IID")
dat <- dat[dat$IID %in% keep, ]
dat <- dat[complete.cases(dat), ]

n_ind <- nrow(dat); n_ev <- sum(dat$EVENT == 1L)
cat("\nindividuals: ", n_ind, "   events: ", n_ev, "\n", sep = "")

## --- genotypes for a random sample of the variants the scan tested ----------
scan <- read.delim(SCAN, stringsAsFactors = FALSE)
scan <- scan[is.finite(scan$P), ]
set.seed(2026)
sel <- scan[sample(nrow(scan), min(N_VAR, nrow(scan))), ]
writeLines(sel$ID, file.path(TMP, "sel.snplist"))
write.table(data.frame(dat$IID, dat$IID), file.path(TMP, "keep.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("exporting ", nrow(sel), " variants ...\n", sep = "")
cmd <- sprintf("plink2 --bfile %s --keep %s --extract %s --threads %d --export A-transpose --out %s",
               shQuote(IN), shQuote(file.path(TMP, "keep.txt")),
               shQuote(file.path(TMP, "sel.snplist")), NCORE, shQuote(file.path(TMP, "g")))
if (system(cmd, ignore.stdout = TRUE) != 0) stop("plink2 export failed:\n  ", cmd)

blk <- read.table(file.path(TMP, "g.traw"), sep = "\t", header = TRUE,
                  comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
G <- as.matrix(blk[, -(1:6), drop = FALSE])
samp <- sub("^[^_]*_", "", colnames(blk)[-(1:6)])
idx <- match(samp, dat$IID)
if (anyNA(idx)) idx <- match(colnames(blk)[-(1:6)], dat$IID)
if (anyNA(idx)) stop("sample mismatch")
dat <- dat[idx, ]

surv_obj <- Surv(dat$TIME, dat$EVENT)

## --- the covariate sets to compare -----------------------------------------
## Unadjusted is the model the Section 7 power calculation actually describes.
## The others are what the scan fits. The gap between them is the answer.
sets <- list(
  "unadjusted"        = character(0),
  "SEX + AGE"         = c("SEX", "AGE"),
  "SEX + AGE + 5 PCs" = c("SEX", "AGE", paste0("PC", 1:5)),
  "SEX + AGE + 10 PCs"= c("SEX", "AGE", paste0("PC", 1:10))
)

fit_set <- function(cols) {
  covars <- if (length(cols)) as.data.frame(dat[, cols, drop = FALSE]) else NULL
  one <- function(i) {
    g <- G[i, ]; good <- is.finite(g)
    if (sum(good) < 20L || length(unique(g[good])) < 2L) return(c(NA, NA))
    d <- if (is.null(covars)) data.frame(G = g[good])
         else cbind(G = g[good], covars[good, , drop = FALSE])
    f <- tryCatch(coxph(surv_obj[good] ~ ., data = d),
                  error = function(e) NULL, warning = function(w) NULL)
    if (is.null(f) || !is.finite(coef(f)[["G"]])) return(c(NA, NA))
    s <- summary(f)$coefficients
    c(s["G", "se(coef)"], s["G", "Pr(>|z|)"])
  }
  r <- if (NCORE > 1L) mclapply(seq_len(nrow(G)), one, mc.cores = NCORE)
       else lapply(seq_len(nrow(G)), one)
  do.call(rbind, lapply(r, function(x) if (is.numeric(x) && length(x) == 2L) x else c(NA, NA)))
}

cat("refitting ", length(sets), " covariate sets on ", NCORE, " core(s) ...\n", sep = "")
res <- lapply(sets, fit_set)

se0 <- median(res[["unadjusted"]][, 1], na.rm = TRUE)

tab <- do.call(rbind, lapply(names(sets), function(nm) {
  m <- res[[nm]]
  se <- median(m[, 1], na.rm = TRUE)
  data.frame(
    covariates       = nm,
    n_parameters     = length(sets[[nm]]) + 1L,
    events_per_par   = round(n_ev / (length(sets[[nm]]) + 1L), 1),
    median_SE        = round(se, 5),
    SE_inflation     = round(se / se0, 4),
    # power scales as (beta/SE)^2, so the detectable effect on the log hazard
    # scale scales linearly with SE
    HR_multiplier    = round(se / se0, 4),
    lambda           = round(lam(m[, 2]), 4),
    stringsAsFactors = FALSE)
}))

cat("\n")
print(tab, row.names = FALSE)

cat("\nHOW TO READ THIS\n")
cat("  median_SE is the precision of the genotype coefficient under each model.\n")
cat("  SE_inflation is that precision relative to the unadjusted fit, which is\n")
cat("  the model the Section 7 power calculation describes. Because power\n")
cat("  depends on beta divided by its standard error, the same figure is the\n")
cat("  multiplier on the smallest detectable log hazard ratio: a value of 1.05\n")
cat("  means the covariate set raises the detectable effect by 5% over what the\n")
cat("  marginal calculation reports.\n\n")
cat("  Read lambda alongside it. A covariate set that costs precision without\n")
cat("  improving calibration is buying nothing, and at ",
    round(n_ev / 13, 1), " events per parameter\n", sep = "")
cat("  for the full model that trade needs to be made deliberately.\n")

adj <- tab$SE_inflation[tab$covariates == "SEX + AGE + 10 PCs"]
if (is.finite(adj)) {
  for (h in c(1.3, 1.5, 2.0)) {
    cat(sprintf("  detectable HR %.1f unadjusted becomes %.2f with ten components\n",
                h, exp(log(h) * adj)))
  }
}

write.table(tab, file.path(OUT_DIR, "pdac_demo_04B_covariate_cost.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWritten: ", file.path(OUT_DIR, "pdac_demo_04B_covariate_cost.tsv"), "\n", sep = "")
unlink(TMP, recursive = TRUE)