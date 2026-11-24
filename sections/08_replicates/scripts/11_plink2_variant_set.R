#!/usr/bin/env Rscript

################################################################################
# Section 8, step 11: PLINK 2 logistic regression and the variants loaded
# beside the one tested (one documented case)
#
# PURPOSE
#   Steps 05 and 06 agree with step 02 to about four significant digits and
#   not to the last one, although the model is the same. This step records
#   the clearest case found, so that it can be quoted and rerun: the same
#   variant, samples, phenotype, covariates and options give a different
#   standard error when other variants are in the same PLINK 2 command.
#
# THE CASE
#   ABO causal variant 9:133273682:A:T, African stratum, replicate 96,
#   SEX, AGE and the components of the stratum, --glm firth-fallback
#   --covar-variance-standardize (the command of step 02):
#     alone    --extract with this variant only
#     window   every variant within 100 kb on either side (--chr/--from-bp/--to-bp)
#     r        glm() in R, double precision, convergence tolerance 1e-12
#   The estimate is reported per copy of the PLINK A1 allele in all three.
#   The cause is not established (our reading: an iterative fit in single
#   precision started from the previous variant; not verified in the PLINK
#   source).
#
# NOT PART OF THE CANONICAL RUN: reads its outputs in results/ and the
# replicates of step 01.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_plink2_variant_set.tsv, pdac_demo_08_plink2_variant_set.txt
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(normalizePath(file.path(script_dir, "..", "..", "..")))
OUT <- "sections/08_replicates/results"; IN <- "results/pca/pdac_demo_02_hwe_filt"
VARIANT <- "9:133273682:A:T"; CHR <- 9L; POS <- 133273682L; HALF <- 100000L; STRATUM <- "afr"; REP <- 96L
ph_f <- sprintf("%s/sim/rep_%03d.pheno", OUT, REP); age_f <- sub("\\.pheno$", ".age", ph_f)
keep_f <- file.path("results/meta", STRATUM, "keep.txt"); cov_f <- file.path("results/meta", STRATUM, "covar.txt")
for (f in c(ph_f, age_f, keep_f, cov_f, paste0(IN, ".bed"))) if (!file.exists(f)) stop("missing input: ", f)
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
tmp <- tempfile("pv"); dir.create(tmp)
cv <- rd(cov_f); age <- read.delim(age_f); cv$AGE <- age$AGE[match(cv$IID, age$IID)]
pcs <- grep("^PC[0-9]+$", names(cv), value = TRUE)
cw <- cv; names(cw)[1] <- paste0("#", names(cw)[1]); cf <- file.path(tmp, "cov"); write.table(cw, cf, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(VARIANT, file.path(tmp, "one"))
plink <- function(tag, sel) {
  out <- file.path(tmp, tag)
  st <- system2("plink2", c("--bfile", IN, sel, "--keep", keep_f, "--pheno", ph_f, "--covar", cf,
                            "--covar-name", paste(c("SEX", "AGE", pcs), collapse = ","), "--covar-variance-standardize", "--threads", "4",
                            "--glm", "firth-fallback", "hide-covar", "cols=+a1freq,+beta,+orbeta,+nobs,+err", "--out", out),
                stdout = FALSE, stderr = FALSE)
  if (st != 0) stop("plink2 failed: ", tag)
  # read as text: an allele column holding only "T" would otherwise become logical
  g <- read.table(Sys.glob(paste0(out, ".*.glm.logistic.hybrid"))[1], header = TRUE, comment.char = "", check.names = FALSE,
                  colClasses = "character"); names(g)[1] <- sub("^#", "", names(g)[1])
  for (cn in c("OBS_CT", "BETA", "SE", "P")) g[[cn]] <- as.numeric(g[[cn]])
  if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
  r <- g[g$ID == VARIANT, ]; if (nrow(r) != 1) stop("variant not in the output of ", tag)
  data.frame(fit = tag, variants_in_command = nrow(g), a1 = r$A1, n = r$OBS_CT, beta = r$BETA, se = r$SE, p = r$P,
             firth = r[[grep("^FIRTH", names(r))[1]]], errcode = r$ERRCODE, stringsAsFactors = FALSE)
}
a <- plink("alone", c("--extract", file.path(tmp, "one")))
w <- plink("window", c("--chr", CHR, "--from-bp", POS - HALF, "--to-bp", POS + HALF))
## R, per copy of the same A1 allele
if (system2("plink2", c("--bfile", IN, "--extract", file.path(tmp, "one"), "--export", "A", "--out", file.path(tmp, "g")),
            stdout = FALSE, stderr = FALSE) != 0) stop("plink2 export failed")
raw <- read.table(file.path(tmp, "g.raw"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gcol <- raw[[7]]; if (sub(".*_", "", names(raw)[7]) != a$a1) gcol <- 2 - gcol     # count the A1 allele of the scan
ids <- read.table(keep_f, stringsAsFactors = FALSE)[[2]]; ph <- read.delim(ph_f)
d <- data.frame(y = as.integer(ph$PHENO[match(ids, ph$IID)] == 2), g = gcol[match(ids, raw$IID)],
                cv[match(ids, cv$IID), c("SEX", "AGE", pcs)])
m <- summary(glm(y ~ ., binomial, d, control = glm.control(epsilon = 1e-12, maxit = 100)))$coefficients["g", ]
r <- data.frame(fit = "r", variants_in_command = 1L, a1 = a$a1, n = sum(!is.na(d$g)), beta = m[[1]], se = m[[2]], p = m[[4]],
                firth = "N", errcode = ".", stringsAsFactors = FALSE)
tab <- rbind(a, w, r)
tab$se_relative_to_r_pct <- 100 * (tab$se / r$se - 1); tab$p_relative_to_r_pct <- 100 * (tab$p / r$p - 1)
tab$se_window_relative_to_alone_pct <- 100 * (w$se / a$se - 1); tab$p_window_relative_to_alone_pct <- 100 * (w$p / a$p - 1)
info <- data.frame(variant = VARIANT, stratum = STRATUM, replicate = REP, window_half_width_bp = HALF, covariates = paste(c("SEX", "AGE", pcs), collapse = ","))
num <- sapply(tab, is.numeric); tw <- tab; tw[num] <- lapply(tw[num], function(v) signif(v, 7))
write.table(cbind(tw, info), file.path(OUT, "pdac_demo_08_plink2_variant_set.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ver <- system2("plink2", "--version", stdout = TRUE)[1]
R <- c("SECTION 8: PLINK 2 LOGISTIC REGRESSION AND THE VARIANTS LOADED BESIDE THE ONE TESTED",
       sprintf("%s; %s; %s, %s stratum, replicate %d; covariates %s; estimate per copy of %s", ver, R.version.string, VARIANT, STRATUM, REP, info$covariates, a$a1), "")
for (j in seq_len(nrow(tab))) with(tab[j, ], R <<- c(R, sprintf("  %-7s variants in the command %4d  n %d  beta %.6f  SE %.6f  P %.4g  Firth %s  error code %s",
  fit, variants_in_command, n, beta, se, p, firth, errcode)))
R <- c(R, "", sprintf("  window against alone: SE %+.2f%%, P %+.1f%%; alone against R: SE %+.4f%%; window against R: SE %+.2f%%",
                      tab$se_window_relative_to_alone_pct[1], tab$p_window_relative_to_alone_pct[1], tab$se_relative_to_r_pct[1], tab$se_relative_to_r_pct[2]),
       "", "  Identical commands give identical results. The standard error, and with it P, depends",
       "  on the other variants in the command; the cause was not established.")
writeLines(R, file.path(OUT, "pdac_demo_08_plink2_variant_set.txt")); cat(R, sep = "\n")
