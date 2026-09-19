#!/usr/bin/env Rscript

################################################################################
# Section 4B: Association testing — survival (Cox proportional hazards)
#
# PURPOSE
#   Run a genome-wide time-to-event scan in the European case series, using Cox
#   proportional hazards regression, and establish that the scan is calibrated
#   before reporting anything from it.
#
# WHY THIS SECTION EXISTS
#   Section 03B computes the power available for a survival scan and finds it
#   very low: with 166 events, the smallest hazard ratio detectable at 80%
#   power and genome-wide significance is 2.12 at a minor allele frequency of
#   0.30. Running the scan anyway is not futile, because the point is to show
#   what an honest underpowered survival analysis looks like and how it should
#   be reported. A null scan is the expected outcome; an *uncontrolled* null
#   scan is not the same thing and must not be reported as one.
#
# WHAT CHANGED IN THIS REVISION, AND WHY
#   The first version of this script reused the principal components built for
#   Section 4A. Those were estimated across all 610 Europeans, cases and
#   controls together. This analysis uses the 217 cases only, and the structure
#   within that subset does not have to lie along the axes estimated in the
#   larger sample. The consequence was measurable: lambda = 1.12 across the
#   whole scan.
#
#   Three diagnostics ruled out the alternatives. The Wald and likelihood-ratio
#   tests agreed to three decimal places (1.2168 against 1.2154 on a 5,000
#   variant subsample), so the inflation was not an artefact of the Wald
#   statistic at low event counts. Lambda was flat across the frequency
#   spectrum (1.121, 1.121, 1.122, 1.110, 1.130 from MAF 0.05 to 0.50), so it
#   was not the asymptotic approximation failing where information is thin;
#   that failure mode rises steeply as frequency falls. Flat inflation across
#   frequency is the signature of confounding.
#
#   The components are therefore recomputed here, inside the survival set. That
#   has to happen in this script rather than upstream, because the survival set
#   is defined by the intersection of European assignment, case status and
#   available follow-up, and does not exist as a group before this point.
#
# WHAT GOVERNS POWER
#   The number of events, not the number of individuals. A survival scan on
#   217 people with 166 events carries roughly the information of a
#   case-control scan with 166 cases. It also constrains the covariate count:
#   166 events against 13 parameters is 12.8 events per parameter, at the
#   conventional floor, which is why the script fits a second model with fewer
#   components and reports the comparison rather than assuming ten are safe.
#
# COMPARISONS PERFORMED
#   1. Components estimated within the 217 cases, against those estimated in
#      the 610-person Section 4A sample. Same individuals, same variants, same
#      outcome; only the covariate axes differ. Answers whether the inflation
#      was caused by borrowing components from a different sample.
#   2. Ten components against a smaller number, both within-set. Answers
#      whether ten can be afforded at 12.8 events per parameter.
#   The primary scan is the within-set PCA at N_PC components. The other two
#   fits are reported alongside it, not instead of it.
#
# INPUT
#   results/pca/pdac_demo_02_hwe_filt.bed/bim/fam
#   results/pca/pdac_demo_02_eur_keep.txt
#   demo_data/survival.txt                       (FID IID TIME EVENT)
#   demo_data/covariates.txt                     (SEX, AGE)
#   results/assoc/pdac_demo_04A_covar.txt        (the Section 4A components)
#   data_processed/highLD_b38.bed
#
# OUTPUT
#   results/assoc/pdac_demo_04B_cox.tsv          primary scan
#   results/assoc/pdac_demo_04B_lambda.tsv
#   results/assoc/pdac_demo_04B_top_hits.tsv
#   results/assoc/pdac_demo_04B_ph_check.tsv     proportional hazards, top hits
#   results/assoc/pdac_demo_04B_pc_comparison.tsv
#   results/assoc/pdac_demo_04B_pca_within.eigenvec / .eigenval
#   results/assoc/pdac_demo_04B_qq.png
#   results/assoc/pdac_demo_04B_summary.txt
#
# USAGE
#   Rscript 01_cox_survival.R [n_pcs] [maf_floor] [n_pcs_reduced]
#   defaults: 10 components, MAF 0.05 within the survival set, 5 for the
#   reduced-covariate comparison
#   Resources are detected at run time; override with GWAS_THREADS, GWAS_MEM_MB.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("survival", quietly = TRUE))
    stop("the 'survival' package is required (it ships with r-recommended)")
  library(survival); library(parallel)
})

## ---------------------------------------------------------------------------
## Locate the project root, the same way the shell scripts do
## ---------------------------------------------------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) {
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
} else getwd()

root <- script_dir
for (up in list(c("..", ".."), c("..", "..", ".."))) {
  cand <- normalizePath(do.call(file.path, c(list(script_dir), as.list(up))),
                        mustWork = FALSE)
  if (dir.exists(file.path(cand, "scripts", "dev"))) { root <- cand; break }
}
if (!dir.exists(file.path(root, "scripts", "dev"))) root <- getwd()
setwd(root)

args      <- commandArgs(trailingOnly = TRUE)
N_PC      <- if (length(args) >= 1) as.integer(args[1]) else 10L
MAF_MIN   <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
N_PC_RED  <- if (length(args) >= 3) as.integer(args[3]) else 5L
PRUNE_MAF <- 0.05

IN        <- "results/pca/pdac_demo_02_hwe_filt"
KEEP_EUR  <- "results/pca/pdac_demo_02_eur_keep.txt"
SURV      <- "demo_data/survival.txt"
COVAR_IN  <- "demo_data/covariates.txt"
COVAR_4A  <- "results/assoc/pdac_demo_04A_covar.txt"
LRLD      <- "data_processed/highLD_b38.bed"
OUT_DIR   <- "results/assoc"
TMP_DIR   <- file.path(tempdir(), "cox_scan")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

for (f in c(paste0(IN, ".bed"), KEEP_EUR, SURV, COVAR_IN, LRLD))
  if (!file.exists(f)) stop("required input missing: ", f)
if (Sys.which("plink2") == "")
  stop("plink2 not on PATH. export PATH=\"$HOME/gwas_tutorial/tools/bin:$PATH\"")

say <- function(...) cat(..., "\n", sep = "")
log_lines <- character(0)
w <- function(...) { s <- paste0(...); log_lines <<- c(log_lines, s); cat(s, "\n", sep = "") }

lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1, lower.tail = FALSE)
}

## ===========================================================================
## Resource detection
##
## A scheduler's allocation always wins over what the hardware reports: on a
## shared node the two differ and only the first is yours to use.
## ===========================================================================
first_env_int <- function(vars) {
  for (v in vars) {
    x <- suppressWarnings(as.integer(Sys.getenv(v, "")))
    if (!is.na(x) && x > 0) return(list(n = x, src = v))
  }
  NULL
}

cgroup_cpu_quota <- function() {
  if (file.exists("/sys/fs/cgroup/cpu.max")) {
    p <- strsplit(trimws(readLines("/sys/fs/cgroup/cpu.max", warn = FALSE)[1]), "\\s+")[[1]]
    if (length(p) == 2L && p[1] != "max") {
      q <- suppressWarnings(as.numeric(p[1])); pr <- suppressWarnings(as.numeric(p[2]))
      if (is.finite(q) && is.finite(pr) && pr > 0) return(q / pr)
    }
  }
  qf <- "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"; pf <- "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
  if (file.exists(qf) && file.exists(pf)) {
    q  <- suppressWarnings(as.numeric(readLines(qf, warn = FALSE)[1]))
    pr <- suppressWarnings(as.numeric(readLines(pf, warn = FALSE)[1]))
    if (is.finite(q) && q > 0 && is.finite(pr) && pr > 0) return(q / pr)
  }
  NA_real_
}

detect_cores <- function() {
  hit <- first_env_int(c("GWAS_THREADS", "SLURM_CPUS_PER_TASK", "NSLOTS",
                         "PBS_NUM_PPN", "OMP_NUM_THREADS"))
  if (!is.null(hit)) return(list(n = hit$n, src = hit$src, physical = NA_integer_))
  phys <- suppressWarnings(detectCores(logical = TRUE))
  if (!is.finite(phys) || phys < 1L) phys <- 1L
  n <- phys; src <- "detectCores"
  q <- cgroup_cpu_quota()
  if (is.finite(q) && floor(q) >= 1 && floor(q) < n) { n <- as.integer(floor(q)); src <- "cgroup quota" }
  list(n = max(1L, n - 1L), src = src, physical = phys)
}

cgroup_mem_limit_mb <- function() {
  for (f in c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
    if (file.exists(f)) {
      v <- trimws(readLines(f, warn = FALSE)[1])
      if (v != "max") {
        x <- suppressWarnings(as.numeric(v))
        if (is.finite(x) && x > 0 && x < 1e15) return(x / 1024^2)
      }
    }
  }
  NA_real_
}

detect_mem_mb <- function() {
  hit <- first_env_int(c("GWAS_MEM_MB", "SLURM_MEM_PER_NODE"))
  if (!is.null(hit)) return(list(mb = hit$n, src = hit$src))
  avail <- NA_real_
  if (file.exists("/proc/meminfo")) {
    mi <- readLines("/proc/meminfo", warn = FALSE)
    ln <- grep("^MemAvailable:", mi, value = TRUE)
    if (!length(ln)) ln <- grep("^MemFree:", mi, value = TRUE)
    if (length(ln)) avail <- as.numeric(gsub("[^0-9]", "", ln[1])) / 1024
  }
  limv <- cgroup_mem_limit_mb()
  if (is.finite(limv)) avail <- if (is.finite(avail)) min(avail, limv) else limv
  if (!is.finite(avail) || avail <= 0) return(list(mb = 2048, src = "fallback"))
  list(mb = avail, src = "system")
}

CORES <- detect_cores(); MEM <- detect_mem_mb()
NCORE <- CORES$n
if (.Platform$OS.type != "unix") NCORE <- 1L

MEM_BUDGET_MB  <- max(128, MEM$mb * 0.20)
# plink2 rejects --memory below 640 MB. Below that, omit the flag rather than
# claim memory the machine does not have, and let plink2 detect for itself.
PLINK_MEM_MB   <- as.integer(min(MEM$mb * 0.50, 16384))
PLINK_MEM_FLAG <- if (PLINK_MEM_MB >= 640L) sprintf("--memory %d ", PLINK_MEM_MB) else ""

say("")
say("=== Resources ===")
say("  cores:  ", NCORE, "  (source: ", CORES$src,
    if (is.finite(CORES$physical)) paste0(", ", CORES$physical, " detected") else "", ")")
say("  memory: ", sprintf("%.0f MB available", MEM$mb), "  (source: ", MEM$src, ")")
say("  chunk budget: ", sprintf("%.0f MB", MEM_BUDGET_MB), " | plink2 ",
    if (nzchar(PLINK_MEM_FLAG)) PLINK_MEM_FLAG else "--memory auto ",
    "--threads ", max(1L, NCORE))
say("  override with GWAS_THREADS and GWAS_MEM_MB")

## ---------------------------------------------------------------------------
## Survival phenotype
##
## Column names vary between datasets, so the time and status columns are
## detected rather than assumed. A hard failure here is better than silently
## fitting the wrong variable.
## ---------------------------------------------------------------------------
surv_raw <- read.delim(SURV, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
names(surv_raw) <- sub("^#", "", names(surv_raw))

pick <- function(df, candidates, what) {
  hit <- names(df)[toupper(names(df)) %in% toupper(candidates)]
  if (!length(hit)) stop("could not identify the ", what, " column in ", SURV,
                         "\n  columns present: ", paste(names(df), collapse = ", "),
                         "\n  expected one of: ", paste(candidates, collapse = ", "))
  hit[1]
}
col_iid   <- pick(surv_raw, c("IID", "ID", "SAMPLE", "SAMPLE_ID"), "sample identifier")
col_time  <- pick(surv_raw, c("TIME", "OS_TIME", "SURV_TIME", "FOLLOWUP", "FOLLOW_UP",
                              "MONTHS", "DAYS", "OS_MONTHS"), "follow-up time")
col_event <- pick(surv_raw, c("EVENT", "STATUS", "OS_EVENT", "DEAD", "DEATH",
                              "OS_STATUS", "CENSOR"), "event indicator")

say("")
say("survival file: ", SURV)
say("  identifier: ", col_iid, " | time: ", col_time, " | event: ", col_event)

surv <- data.frame(IID   = as.character(surv_raw[[col_iid]]),
                   TIME  = as.numeric(surv_raw[[col_time]]),
                   EVENT = as.integer(surv_raw[[col_event]]),
                   stringsAsFactors = FALSE)
surv <- surv[is.finite(surv$TIME) & surv$TIME > 0 & surv$EVENT %in% c(0L, 1L), ]

## ---------------------------------------------------------------------------
## Define the survival set: European, has follow-up, has sex and age
## ---------------------------------------------------------------------------
base_cov <- read.delim(COVAR_IN, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
names(base_cov)[1] <- sub("^#", "", names(base_cov)[1])
base_cov$IID <- as.character(base_cov$IID)
if (!all(c("SEX", "AGE") %in% names(base_cov)))
  stop("expected SEX and AGE in ", COVAR_IN)

keep_eur <- read.table(KEEP_EUR, stringsAsFactors = FALSE)
keep_ids <- as.character(keep_eur[[ncol(keep_eur)]])

dat <- merge(surv, base_cov[, c("IID", "SEX", "AGE")], by = "IID")
dat <- dat[dat$IID %in% keep_ids, ]
dat <- dat[complete.cases(dat), ]
if (nrow(dat) < 50) stop("fewer than 50 individuals with survival data and covariates")

n_ind <- nrow(dat); n_event <- sum(dat$EVENT == 1L); n_censor <- n_ind - n_event

## ---------------------------------------------------------------------------
## Principal components WITHIN the survival set
##
## Same pruning parameters as everywhere else in the pipeline, but --keep is
## applied first so that the frequency floor and the linkage disequilibrium
## estimates both use this set's own genotypes.
## ---------------------------------------------------------------------------
say("")
say("computing components within the survival set ...")

keep_file <- file.path(TMP_DIR, "keep.txt")
write.table(data.frame(FID = dat$IID, IID = dat$IID), keep_file,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

prune_pref <- file.path(TMP_DIR, "prune")
cmd <- sprintf(paste0("plink2 --bfile %s --keep %s --autosome --maf %s ",
                      "--exclude bed0 %s --indep-pairwise 50 5 0.2 ",
                      "--threads %d ", PLINK_MEM_FLAG, "--out %s"),
               shQuote(IN), shQuote(keep_file), format(PRUNE_MAF), shQuote(LRLD),
               max(1L, NCORE), shQuote(prune_pref))
if (system(cmd, ignore.stdout = TRUE) != 0)
  stop("pruning failed; run by hand to see why:\n  ", cmd)
n_prune_var <- length(readLines(paste0(prune_pref, ".prune.in")))

pca_pref <- file.path(TMP_DIR, "pca")
cmd <- sprintf(paste0("plink2 --bfile %s --keep %s --extract %s --pca %d ",
                      "--threads %d ", PLINK_MEM_FLAG, "--out %s"),
               shQuote(IN), shQuote(keep_file), shQuote(paste0(prune_pref, ".prune.in")),
               N_PC, max(1L, NCORE), shQuote(pca_pref))
if (system(cmd, ignore.stdout = TRUE) != 0)
  stop("within-set PCA failed; run by hand to see why:\n  ", cmd)

pcs_within <- read.delim(paste0(pca_pref, ".eigenvec"), check.names = FALSE, comment.char = "")
names(pcs_within)[1] <- sub("^#", "", names(pcs_within)[1])
pcs_within$IID <- as.character(pcs_within$IID)
file.copy(paste0(pca_pref, ".eigenvec"),
          file.path(OUT_DIR, "pdac_demo_04B_pca_within.eigenvec"), overwrite = TRUE)
file.copy(paste0(pca_pref, ".eigenval"),
          file.path(OUT_DIR, "pdac_demo_04B_pca_within.eigenval"), overwrite = TRUE)

pc_names <- paste0("PC", seq_len(N_PC))
dat <- merge(dat, pcs_within[, c("IID", pc_names)], by = "IID")
dat <- dat[complete.cases(dat), ]
n_ind <- nrow(dat); n_event <- sum(dat$EVENT == 1L); n_censor <- n_ind - n_event

## The Section 4A components, for the comparison fit. Missing is not fatal:
## the comparison is informative but not required.
pcs_4A <- NULL
if (file.exists(COVAR_4A)) {
  tmp <- read.delim(COVAR_4A, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
  names(tmp)[1] <- sub("^#", "", names(tmp)[1])
  tmp$IID <- as.character(tmp$IID)
  if (all(pc_names %in% names(tmp))) {
    pcs_4A <- tmp[match(dat$IID, tmp$IID), pc_names, drop = FALSE]
    if (anyNA(pcs_4A[[1]])) pcs_4A <- NULL
  }
}

n_par <- length(c("SEX", "AGE", pc_names)) + 1L

w("Cox proportional hazards survival scan")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("")
w("ANALYSIS SET")
w("  European cases with follow-up, sex and age: ", n_ind)
w("  events: ", n_event, "   censored: ", n_censor)
w("  minor allele frequency floor, evaluated in this set: ", MAF_MIN)
w("")
w("COVARIATES")
w("  components computed within this set, on ", n_prune_var, " pruned variants")
w("  primary model: SEX, AGE, PC1-PC", N_PC, "  (", n_par, " parameters)")
w("  events per parameter: ", sprintf("%.1f", n_event / n_par))
w("")
w("  The number governing power is the event count, not the sample size.")
w("")

## ---------------------------------------------------------------------------
## Export genotypes, variant-major, restricted to this set
## ---------------------------------------------------------------------------
traw_prefix <- file.path(TMP_DIR, "geno")
cmd <- sprintf(paste0("plink2 --bfile %s --keep %s --chr 1-22 --maf %s ",
                      "--threads %d ", PLINK_MEM_FLAG, "--export A-transpose --out %s"),
               shQuote(IN), shQuote(keep_file), format(MAF_MIN),
               max(1L, NCORE), shQuote(traw_prefix))
say("")
say("exporting genotypes ...")
if (system(cmd, ignore.stdout = TRUE) != 0)
  stop("plink2 export failed; run by hand to see why:\n  ", cmd)
traw <- paste0(traw_prefix, ".traw")
if (!file.exists(traw)) stop("expected export not found: ", traw)

con <- file(traw, "r")
hdr <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
meta_n <- 6L
sample_cols <- hdr[-seq_len(meta_n)]
idx <- match(sub("^[^_]*_", "", sample_cols), dat$IID)
if (anyNA(idx)) idx <- match(sample_cols, dat$IID)
if (anyNA(idx)) { close(con); stop("sample mismatch between the export and the analysis set") }
dat <- dat[idx, ]
if (!is.null(pcs_4A)) pcs_4A <- pcs_4A[idx, , drop = FALSE]

surv_obj <- Surv(dat$TIME, dat$EVENT)
base_df  <- as.data.frame(dat[, c("SEX", "AGE", pc_names)])
ev_flag  <- dat$EVENT == 1L
n_samp   <- nrow(dat)

bytes_per_variant <- n_samp * 8 * 4
CHUNK <- as.integer(MEM_BUDGET_MB * 1024^2 / (bytes_per_variant * max(1L, NCORE)))
CHUNK <- max(1000L, min(CHUNK, 25000L))
say("  chunk size: ", CHUNK, " variants (", n_samp, " samples per row)")

col_classes <- c(rep("character", 2), "NULL", "integer",
                 rep("character", 2), rep("numeric", length(sample_cols)))

NA_ROW <- c(freq = NA_real_, n = NA_real_, ev = NA_real_, beta = NA_real_,
            se = NA_real_, z = NA_real_, p = NA_real_, code = 1)

make_fitter <- function(covars) {
  force(covars)
  function(g) {
    good <- is.finite(g)
    if (sum(good) < 20L || length(unique(g[good])) < 2L) {
      r <- NA_ROW; r["code"] <- 2; return(r)
    }
    d <- cbind(G = g[good], covars[good, , drop = FALSE])
    fit <- tryCatch(coxph(surv_obj[good] ~ ., data = d),
                    error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit) || !is.finite(coef(fit)[["G"]])) {
      r <- NA_ROW
      r["freq"] <- mean(g[good]) / 2; r["n"] <- sum(good); r["ev"] <- sum(ev_flag[good])
      r["code"] <- 3; return(r)
    }
    s <- summary(fit)$coefficients
    c(freq = mean(g[good]) / 2, n = sum(good), ev = sum(ev_flag[good]),
      beta = s["G", "coef"], se = s["G", "se(coef)"],
      z = s["G", "z"], p = s["G", "Pr(>|z|)"], code = 0)
  }
}

## Models fitted in the same pass over the genotypes, so that all of them see
## identical variants and the comparison cannot be confounded by a differing
## variant set.
models <- list(primary = make_fitter(base_df))
models$reduced <- make_fitter(base_df[, c("SEX", "AGE", paste0("PC", seq_len(N_PC_RED)))])
if (!is.null(pcs_4A)) {
  df_4A <- as.data.frame(cbind(dat[, c("SEX", "AGE")], pcs_4A))
  models$pc_from_4A <- make_fitter(df_4A)
}

apply_chunk <- if (NCORE > 1L) {
  function(rows, f) mclapply(rows, f, mc.cores = NCORE, mc.preschedule = TRUE)
} else function(rows, f) lapply(rows, f)

## ---------------------------------------------------------------------------
## Stream and fit
## ---------------------------------------------------------------------------
res <- list(); alt_p <- list(); k <- 0L; n_read <- 0L
t0 <- Sys.time()

repeat {
  block <- tryCatch(read.table(con, sep = "\t", nrows = CHUNK, header = FALSE,
                               stringsAsFactors = FALSE, comment.char = "",
                               colClasses = col_classes),
                    error = function(e) NULL)
  if (is.null(block) || !nrow(block)) break

  G <- as.matrix(block[, -(1:5), drop = FALSE])
  rows <- lapply(seq_len(nrow(G)), function(i) G[i, ])

  fits <- apply_chunk(rows, models$primary)
  bad <- !vapply(fits, function(x) is.numeric(x) && length(x) == 8L, logical(1))
  if (any(bad)) fits[bad] <- list({ r <- NA_ROW; r["code"] <- 4; r })
  M <- do.call(rbind, fits)

  k <- k + 1L
  res[[k]] <- data.frame(
    CHROM = block[[1]], POS = block[[3]], ID = block[[2]],
    COUNTED = block[[4]], ALT = block[[5]],
    A1_FREQ = M[, "freq"], N = M[, "n"], EVENTS = M[, "ev"],
    BETA = M[, "beta"], SE = M[, "se"], HR = exp(M[, "beta"]),
    Z = M[, "z"], P = M[, "p"],
    ERRCODE = c(".", "MONOMORPHIC_OR_TOO_FEW", "FIT_FAILED",
                "WORKER_ERROR")[M[, "code"] + 1L],
    stringsAsFactors = FALSE)

  alt <- list()
  for (nm in setdiff(names(models), "primary")) {
    f2 <- apply_chunk(rows, models[[nm]])
    alt[[nm]] <- vapply(f2, function(x)
      if (is.numeric(x) && length(x) == 8L) x[["p"]] else NA_real_, numeric(1))
  }
  alt_p[[k]] <- as.data.frame(alt)

  n_read <- n_read + nrow(block)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\r  %d variants | %.1f min | %.0f variants/s",
              n_read, el, n_read / max(el * 60, 1e-9)))
  flush.console()
}
close(con); cat("\n")

if (!length(res)) stop("no variants were read from ", traw)
r <- do.call(rbind, res)
alt_df <- do.call(rbind, alt_p)
ord <- order(as.integer(r$CHROM), r$POS)
r <- r[ord, ]; alt_df <- alt_df[ord, , drop = FALSE]
n_skip <- sum(r$ERRCODE != ".")

## ---------------------------------------------------------------------------
## Calibration, top hits, proportional hazards
## ---------------------------------------------------------------------------
p_ok <- r$P[is.finite(r$P) & r$P > 0 & r$P <= 1]
lambda <- lam(r$P)
# No standard lambda_1000 exists for a survival scan. Inflation scales with the
# information in the study, which is carried by the events, so this is lambda
# rescaled to a nominal 1000 events. It is not the case-control lambda_1000 and
# the two are not comparable.
lambda_ev <- 1 + (lambda - 1) * 1000 / n_event

w("SCAN")
w("  variants tested: ", nrow(r))
w("  variants with a valid P: ", length(p_ok))
w("  skipped (monomorphic in this set, or the fit failed): ", n_skip)
w("  genomic inflation factor lambda: ", sprintf("%.4f", lambda))
w("  lambda per 1000 events: ", sprintf("%.4f", lambda_ev))
w("")

sig <- sum(p_ok < 5e-8); sugg <- sum(p_ok < 1e-5)
w("  genome-wide significant (P < 5e-8): ", sig)
w("  suggestive (P < 1e-5): ", sugg)
w("  expected under the null at 1e-5: ", sprintf("%.1f", length(p_ok) * 1e-5))

top <- r[is.finite(r$P), ]
top <- top[order(top$P), ][seq_len(min(20L, nrow(top))), ]
if (nrow(top)) w("  lead variant: ", top$ID[1], "  HR ", sprintf("%.3f", top$HR[1]),
                 "  P ", format(top$P[1], digits = 3))
w("")

## Covariate comparison. Same individuals, same variants, same outcome; only
## the covariate axes differ, so any difference in lambda is attributable to
## the covariates alone.
cmp <- data.frame(
  model = "primary: within-set PCs, PC1-PC" ,
  stringsAsFactors = FALSE)
cmp <- data.frame(
  model  = c(paste0("within-set PCs, PC1-PC", N_PC),
             paste0("within-set PCs, PC1-PC", N_PC_RED),
             if ("pc_from_4A" %in% names(alt_df)) paste0("Section 4A PCs, PC1-PC", N_PC)),
  n_par  = c(n_par, N_PC_RED + 3L,
             if ("pc_from_4A" %in% names(alt_df)) n_par),
  events_per_par = NA_real_,
  lambda = c(lambda, lam(alt_df$reduced),
             if ("pc_from_4A" %in% names(alt_df)) lam(alt_df$pc_from_4A)),
  suggestive = c(sugg, sum(alt_df$reduced < 1e-5, na.rm = TRUE),
                 if ("pc_from_4A" %in% names(alt_df))
                   sum(alt_df$pc_from_4A < 1e-5, na.rm = TRUE)),
  stringsAsFactors = FALSE)
cmp$events_per_par <- round(n_event / cmp$n_par, 1)
cmp$lambda <- round(cmp$lambda, 4)

w("COVARIATE COMPARISON")
w("  Same ", n_ind, " individuals, same ", nrow(r), " variants, same outcome.")
w("  Only the covariate axes differ, so any difference in lambda is theirs.")
w("")
for (i in seq_len(nrow(cmp)))
  w(sprintf("  %-34s  %2d par  %5.1f ev/par  lambda %.4f  suggestive %d",
            cmp$model[i], cmp$n_par[i], cmp$events_per_par[i],
            cmp$lambda[i], cmp$suggestive[i]))
w("")
if ("pc_from_4A" %in% names(alt_df)) {
  d <- lam(alt_df$pc_from_4A) - lambda
  if (is.finite(d) && d > 0.03) {
    w("  The components borrowed from Section 4A leave more inflation than those")
    w("  computed here. They were estimated across cases and controls together,")
    w("  and describe axes of variation in that larger sample rather than in")
    w("  this one. Structure specific to the case series is not on those axes.")
  } else if (is.finite(d)) {
    w("  The two sets of components leave similar inflation, so borrowing them")
    w("  from Section 4A was not the source of the problem. If lambda remains")
    w("  above about 1.05, look elsewhere: a frequency floor set by events")
    w("  rather than by sample size, or cryptic relatedness within the cases.")
  }
}
w("")

ph <- NULL
if (nrow(top)) {
  ph_rows <- lapply(seq_len(nrow(top)), function(j) {
    line <- suppressWarnings(system(
      sprintf("grep -m1 -F %s %s", shQuote(paste0("\t", top$ID[j], "\t")), shQuote(traw)),
      intern = TRUE, ignore.stderr = TRUE))
    if (!length(line)) return(NULL)
    v <- strsplit(line[1], "\t", fixed = TRUE)[[1]]
    g <- suppressWarnings(as.numeric(v[-seq_len(meta_n)]))
    good <- is.finite(g)
    fit <- tryCatch(coxph(surv_obj[good] ~ .,
                          data = cbind(G = g[good], base_df[good, , drop = FALSE])),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    z <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if (is.null(z)) return(NULL)
    data.frame(ID = top$ID[j], P_assoc = top$P[j],
               PH_P_genotype = z$table["G", "p"],
               PH_P_global   = z$table["GLOBAL", "p"], stringsAsFactors = FALSE)
  })
  ph <- do.call(rbind, Filter(Negate(is.null), ph_rows))
}

w("HOW TO READ THIS SCAN")
w("  Section 03B gives the detectable effect sizes for ", n_event, " events. A null")
w("  result at this event count is the expected outcome, not a failure of the")
w("  analysis. But a null scan and an uncontrolled scan are different things:")
w("  an excess of suggestive results is only interpretable once lambda is near")
w("  one. Report the event count, the covariates, lambda and the minimum")
w("  detectable hazard ratio together, so that a reader can tell an")
w("  uninformative scan from a negative one.")
w("")
if (is.finite(lambda) && lambda > 1.05) {
  w("  WARNING: lambda is ", sprintf("%.3f", lambda), ". The suggestive count above is")
  w("  not yet interpretable as signal. Resolve the inflation before reporting.")
}
w("")
w("RESOURCES USED")
w("  cores: ", NCORE, " (", CORES$src, ")   chunk: ", CHUNK, " variants")
w("  elapsed: ", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---------------------------------------------------------------------------
## Write
## ---------------------------------------------------------------------------
f_all <- file.path(OUT_DIR, "pdac_demo_04B_cox.tsv")
write.table(r, f_all, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top, file.path(OUT_DIR, "pdac_demo_04B_top_hits.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cmp, file.path(OUT_DIR, "pdac_demo_04B_pc_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(
  metric = c("n_individuals", "n_events", "n_censored", "n_parameters",
             "events_per_parameter", "pca_variants", "variants_tested",
             "variants_skipped", "lambda", "lambda_per_1000_events",
             "genome_wide_sig", "suggestive", "cores_used", "chunk_size"),
  value  = c(n_ind, n_event, n_censor, n_par, round(n_event / n_par, 1),
             n_prune_var, nrow(r), n_skip, round(lambda, 4),
             round(lambda_ev, 4), sig, sugg, NCORE, CHUNK)),
  file.path(OUT_DIR, "pdac_demo_04B_lambda.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(ph)) write.table(ph, file.path(OUT_DIR, "pdac_demo_04B_ph_check.tsv"),
                              sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(log_lines, file.path(OUT_DIR, "pdac_demo_04B_summary.txt"))

png(file.path(OUT_DIR, "pdac_demo_04B_qq.png"), width = 1600, height = 1600, res = 200)
o <- -log10(sort(p_ok)); e <- -log10(ppoints(length(o)))
plot(e, o, pch = 20, cex = 0.4, col = "#0072B2",
     xlab = expression(Expected~-log[10](P)), ylab = expression(Observed~-log[10](P)),
     main = sprintf("Cox survival scan, %d events (lambda = %.3f)", n_event, lambda))
abline(0, 1, col = "grey40")
invisible(dev.off())

say("")
say("Written:")
for (f in c(f_all,
            file.path(OUT_DIR, "pdac_demo_04B_top_hits.tsv"),
            file.path(OUT_DIR, "pdac_demo_04B_pc_comparison.tsv"),
            file.path(OUT_DIR, "pdac_demo_04B_lambda.tsv"),
            if (!is.null(ph)) file.path(OUT_DIR, "pdac_demo_04B_ph_check.tsv"),
            file.path(OUT_DIR, "pdac_demo_04B_pca_within.eigenvec"),
            file.path(OUT_DIR, "pdac_demo_04B_qq.png"),
            file.path(OUT_DIR, "pdac_demo_04B_summary.txt"))) say("   ", f)
unlink(TMP_DIR, recursive = TRUE)