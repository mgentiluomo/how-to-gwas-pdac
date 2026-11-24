#!/usr/bin/env Rscript

################################################################################
# Section 8, step 05: Wald, score, likelihood-ratio and Firth tests at the
# causal variants, and their calibration under the null
#
# PURPOSE
#   Step 03 finds realised power below the prediction of the Section 3B
#   formula, and step 04 traces it to the Wald statistic: its standard error
#   grows with the estimate, so z varies less than the estimate does. The
#   Section 3B formula uses the information at the null, which is the
#   information of the score test. This step asks whether another test on the
#   same 200 replicates reaches the predicted power, that is whether the
#   shortfall belongs to the formula or to the Wald test PLINK 2 applies.
#
# WHAT IT DOES
#   1. CAUSAL VARIANTS. For every replicate and each of the three causal
#      variants, in R, with the models of step 02:
#        European scan   status ~ G + SEX + AGE + PC1..PCk  (k of 04A step 01b)
#        meta-analysis   the three strata, each with SEX, AGE and its own
#                        components (results/meta/<stratum>/covar.txt)
#      four tests of G:
#        wald    estimate / standard error (what PLINK 2 --glm reports); in
#                the meta-analysis the fixed-effect inverse-variance
#                combination of the three stratum estimates, as step 02
#        score   U^2 / V from the model without G; in the meta-analysis the
#                sum of the stratum scores over the sum of their variances
#        lrt     twice the difference in log likelihood; in the meta-analysis
#                the stratified model with one common effect of G and all
#                covariate effects stratum-specific, against the same without G
#        firth   penalised likelihood ratio (Firth 1993; Heinze and Schemper
#                2002), same models as lrt
#      Missing genotypes: complete cases per variant, as PLINK 2.
#      CHECKS: the Wald chi-square must reproduce step 02 (PLINK 2) within
#      2% in the European scan and in the meta-analysis. It is not exact:
#      at a large effect the standard error PLINK 2 reports moves with the
#      variants loaded beside the one tested, which we attribute to its
#      iterative fit in single precision started from the previous variant
#      (not verified in the PLINK source) (ABO, African stratum,
#      replicate 96: the same estimate, standard error 0.2379 when the variant
#      is tested alone and 0.2368 within a 200 kb window; R in double
#      precision gives 0.2379). The Firth function is compared with logistf
#      when that package is installed (difference printed; logistf is not
#      needed to run the step).
#   2. NULL CALIBRATION, European scan. N_NULL variants drawn once (seed
#      NULL_SEED) from the LD-pruned European set, MAF at least 0.05, at least
#      2 Mb from every causal variant, tested with wald, score and lrt in
#      every replicate. They carry no effect, so lambda and the tail
#      fractions of each test are read against their nominal values.
#      Firth is left out of this part.
#
# NOT PART OF THE CANONICAL RUN: reads its outputs in results/, as step 02.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_tests_per_replicate.tsv   replicate x locus x analysis: the four tests
#   pdac_demo_08_tests_power.tsv           realised power per test, with 95% CI
#   pdac_demo_08_tests_null.tsv            null calibration per test
#   pdac_demo_08_tests_summary.txt         plain-text report
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
STRATA <- c("eur", "afr", "eas")
CAUSAL <- c(G1 = "9:133273682:A:T", G2 = "5:1286401:C:A", G3 = "13:73362027:C:T")
LABEL <- c(G1 = "ABO", G2 = "TERT/CLPTM1L", G3 = "13q22.1")
OR_GEN <- c(G1 = 2.60, G2 = 1.50, G3 = 1.97)
TESTS <- c("wald", "score", "lrt", "firth")
ALPHA <- 5e-8
N_NULL <- 2000L; NULL_SEED <- 20261124L; NULL_DIST <- 2e6; NULL_MAF <- 0.05
N_PAR <- as.integer(Sys.getenv("GWAS_08_PARALLEL", "6"))
reps <- sort(Sys.glob(file.path(SIM, "rep_*.pheno")))
if (!length(reps)) stop("no replicates in ", SIM, ": run step 01")
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
k <- as.integer(read.table("results/assoc/pdac_demo_04A_npc_chosen.txt")[1, 2])
per <- read.delim(file.path(OUT, "pdac_demo_08_per_replicate.tsv"), stringsAsFactors = FALSE)
tmp <- file.path(tempdir(), "tests"); dir.create(tmp, showWarnings = FALSE)

## ---- genotypes ------------------------------------------------------------------------
## causal variants: allele counted by the generator (the risk allele), missing kept
writeLines(unname(CAUSAL), file.path(tmp, "c.snps"))
if (system2("plink2", c("--bfile", IN, "--extract", file.path(tmp, "c.snps"), "--export", "A", "--out", file.path(tmp, "c")),
            stdout = FALSE, stderr = FALSE) != 0) stop("plink2 export failed")
system2("plink2", c("--bfile", "demo_data/pdac_demo", "--extract", file.path(tmp, "c.snps"), "--recode", "A",
                    "--out", file.path(tmp, "gen")), stdout = FALSE, stderr = FALSE)
gen_hdr <- names(read.table(file.path(tmp, "gen.raw"), header = TRUE, check.names = FALSE, nrows = 1))
raw <- read.table(file.path(tmp, "c.raw"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
col_of <- function(h, id) grep(paste0("^", gsub("([.:])", "\\\\\\1", id), "_"), h, value = TRUE)
Gc <- sapply(CAUSAL, function(id) { cn <- col_of(names(raw), id); x <- raw[[cn]]
  if (sub(".*_", "", cn) != sub(".*_", "", col_of(gen_hdr, id))) 2 - x else x })
rownames(Gc) <- raw$IID

## ---- designs ---------------------------------------------------------------------------
cov_e <- rd(EUR_COV)
eur_ids <- read.table(EUR_KEEP, stringsAsFactors = FALSE)[[2]]
design <- list(eur_4a = list(ids = eur_ids, cov = cov_e[match(eur_ids, cov_e$IID), c("IID", "SEX", paste0("PC", seq_len(k)))]))
for (g in STRATA) { ids <- read.table(file.path("results/meta", g, "keep.txt"), stringsAsFactors = FALSE)[[2]]
  cv <- rd(file.path("results/meta", g, "covar.txt"))
  design[[g]] <- list(ids = ids, cov = cv[match(ids, cv$IID), c("IID", "SEX", grep("^PC[0-9]+$", names(cv), value = TRUE))]) }
for (d in design) if (anyNA(d$cov)) stop("covariates missing for individuals of a keep file")

## ---- the tests --------------------------------------------------------------------------
ll <- function(X, y, b) { eta <- drop(X %*% b); sum(y * eta - log1p(exp(eta))) }
fit_ml <- function(X, y) { m <- suppressWarnings(glm.fit(X, y, family = binomial(), control = glm.control(epsilon = 1e-10, maxit = 50)))
  mu <- m$fitted.values; list(b = m$coefficients, mu = mu, w = mu * (1 - mu), dev = m$deviance) }
firth <- function(X, y, fix = integer()) {
  p <- ncol(X); b <- numeric(p); free <- setdiff(seq_len(p), fix)
  pl <- function(b) { eta <- drop(X %*% b); mu <- plogis(eta)
    sum(y * eta - log1p(exp(eta))) + 0.5 * as.numeric(determinant(crossprod(X * sqrt(mu * (1 - mu))))$modulus) }
  cur <- pl(b)
  for (it in 1:200) {
    mu <- plogis(drop(X %*% b)); XW <- X * sqrt(mu * (1 - mu)); I <- crossprod(XW)
    h <- rowSums((XW %*% solve(I)) * XW)
    U <- drop(crossprod(X, y - mu + h * (0.5 - mu)))
    step <- solve(I[free, free, drop = FALSE], U[free])
    nb <- b; nb[free] <- b[free] + step; new <- pl(nb); hv <- 0
    while (new < cur && hv < 20) { step <- step / 2; nb[free] <- b[free] + step; new <- pl(nb); hv <- hv + 1 }
    b <- nb; conv <- max(abs(step)) < 1e-9 && abs(new - cur) < 1e-10; cur <- new
    if (conv) break
  }
  list(b = b, pl = cur)
}
## all four tests of the first column of X (G) given the other columns
tests <- function(X, y, do_firth = TRUE) {
  f1 <- fit_ml(X, y); f0 <- fit_ml(X[, -1, drop = FALSE], y)
  se <- sqrt(solve(crossprod(X * sqrt(f1$w)))[1, 1]); g <- X[, 1]; X0 <- X[, -1, drop = FALSE]
  U <- sum(g * (y - f0$mu)); A <- crossprod(X0 * f0$w, g)
  V <- sum(f0$w * g^2) - drop(crossprod(A, solve(crossprod(X0 * sqrt(f0$w)), A)))
  o <- c(b_wald = unname(f1$b[1]), se_wald = se, chi_wald = unname((f1$b[1] / se)^2),
         score_u = U, score_v = V, chi_score = U^2 / V, chi_lrt = f0$dev - f1$dev)
  if (do_firth) { a <- firth(X, y); n0 <- firth(X, y, fix = 1L)
    o <- c(o, b_firth = a$b[1], chi_firth = max(0, 2 * (a$pl - n0$pl))) }
  o
}
block <- function(Xs) { nc <- sapply(Xs, ncol); nr <- sapply(Xs, nrow)
  M <- matrix(0, sum(nr), sum(nc)); r0 <- 0; c0 <- 0
  for (i in seq_along(Xs)) { M[r0 + seq_len(nr[i]), c0 + seq_len(nc[i])] <- Xs[[i]]; r0 <- r0 + nr[i]; c0 <- c0 + nc[i] }; M }

one <- function(f) {
  i <- as.integer(sub(".*rep_([0-9]+)\\.pheno$", "\\1", f))
  ph <- read.delim(f); age <- read.delim(sub("\\.pheno$", ".age", f))
  dat <- lapply(design, function(d) list(y = as.integer(ph$PHENO[match(d$ids, ph$IID)] == 2),
    X0 = cbind(1, as.matrix(d$cov[, -1]), AGE = age$AGE[match(d$ids, age$IID)]), ids = d$ids))
  do.call(rbind, lapply(names(CAUSAL), function(L) {
    part <- lapply(dat, function(d) { g <- Gc[d$ids, L]; ok <- !is.na(g); list(X = cbind(g[ok], d$X0[ok, , drop = FALSE]), y = d$y[ok]) })
    e <- tests(part$eur_4a$X, part$eur_4a$y)
    s <- lapply(part[STRATA], function(p) tests(p$X, p$y, do_firth = FALSE))
    bs <- sapply(s, `[[`, "b_wald"); ws <- 1 / sapply(s, `[[`, "se_wald")^2
    b_iv <- sum(ws * bs) / sum(ws); se_iv <- sqrt(1 / sum(ws))
    Xst <- cbind(unlist(lapply(part[STRATA], function(p) p$X[, 1])), block(lapply(part[STRATA], function(p) p$X[, -1, drop = FALSE])))
    yst <- unlist(lapply(part[STRATA], `[[`, "y"))
    f1 <- fit_ml(Xst, yst); a <- firth(Xst, yst); n0 <- firth(Xst, yst, fix = 1L)
    m <- c(b_wald = b_iv, se_wald = se_iv, chi_wald = (b_iv / se_iv)^2,
           score_u = sum(sapply(s, `[[`, "score_u")), score_v = sum(sapply(s, `[[`, "score_v")),
           chi_score = sum(sapply(s, `[[`, "score_u"))^2 / sum(sapply(s, `[[`, "score_v")),
           chi_lrt = sum(sapply(part[STRATA], function(p) fit_ml(p$X[, -1, drop = FALSE], p$y)$dev)) - f1$dev,
           b_firth = a$b[1], chi_firth = max(0, 2 * (a$pl - n0$pl)))
    rbind(data.frame(replicate = i, locus = L, analysis = "eur_4a", n = length(part$eur_4a$y), t(e)),
          data.frame(replicate = i, locus = L, analysis = "meta", n = length(yst), t(m)))
  }))
}
t0 <- Sys.time()
res <- mclapply(reps, function(f) tryCatch(one(f), error = function(e) e), mc.cores = N_PAR, mc.preschedule = FALSE)
bad <- which(sapply(res, function(x) inherits(x, "error") || inherits(x, "try-error")))
if (length(bad)) stop("replicates failed: ", paste(basename(reps[bad]), collapse = ", "), ": ", as.character(res[[bad[1]]]))
d <- do.call(rbind, res); d <- d[order(d$replicate, d$locus, d$analysis), ]
for (t in TESTS) d[[paste0("p_", t)]] <- pchisq(d[[paste0("chi_", t)]], 1, lower.tail = FALSE)
## Wald against step 02 (PLINK 2)
mm <- match(paste(d$replicate, d$locus, d$analysis), paste(per$replicate, per$locus, per$analysis))
d$p_step02 <- per$p[mm]; d$predicted_power <- per$predicted_power[mm]
d$chi_step02 <- qchisq(d$p_step02, 1, lower.tail = FALSE)
dchk <- tapply(abs(d$chi_wald / d$chi_step02 - 1), d$analysis, max); dmed <- tapply(abs(d$chi_wald / d$chi_step02 - 1), d$analysis, median)
chk_line <- sprintf("Wald chi-square here (R, double precision) against step 02 (PLINK 2): relative difference, max %.2g (median %.2g) European scan, max %.2g (median %.2g) meta-analysis",
                    dchk[["eur_4a"]], dmed[["eur_4a"]], dchk[["meta"]], dmed[["meta"]])
cat(chk_line, "\n"); if (any(dchk > 0.02)) stop("the Wald test does not reproduce step 02")
## Firth against logistf, first replicate, European scan
firth_line <- "Firth function not compared with logistf (package not installed)"
if (requireNamespace("logistf", quietly = TRUE)) {
  ph <- read.delim(reps[1]); age <- read.delim(sub("\\.pheno$", ".age", reps[1])); dd <- design$eur_4a
  df <- data.frame(y = as.integer(ph$PHENO[match(dd$ids, ph$IID)] == 2), dd$cov[, -1], AGE = age$AGE[match(dd$ids, age$IID)])
  dl <- sapply(names(CAUSAL), function(L) { x <- cbind(g = Gc[dd$ids, L], df); x <- x[!is.na(x$g), ]
    m <- logistf::logistf(y ~ ., data = x, pl = TRUE, control = logistf::logistf.control(maxit = 200, lconv = 1e-10, gconv = 1e-10, xconv = 1e-10),
                          plcontrol = logistf::logistpl.control(maxit = 200))
    mine <- d[d$replicate == 1 & d$locus == L & d$analysis == "eur_4a", ]
    c(abs(coef(m)[["g"]] - mine$b_firth), abs(log10(m$prob[["g"]]) - log10(mine$p_firth))) })
  firth_line <- sprintf("Firth function against logistf %s, replicate 1, European scan: max |estimate difference| %.2g, max |log10 P difference| %.2g",
                        as.character(packageVersion("logistf")), max(dl[1, ]), max(dl[2, ]))
}
cat(firth_line, "\n")

wilson <- function(x, n) { p <- x / n; z <- 1.959964
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
crit <- qchisq(1 - ALPHA, 1)
pow <- do.call(rbind, lapply(split(d, list(d$locus, d$analysis)), function(x) do.call(rbind, lapply(TESTS, function(t) {
  chi <- x[[paste0("chi_", t)]]; sig <- chi > crit; ci <- wilson(sum(sig), nrow(x)); z <- sqrt(chi) * sign(if (t == "score") x$score_u else x$b_wald)
  data.frame(locus = x$locus[1], label = LABEL[[x$locus[1]]], or_generative = OR_GEN[[x$locus[1]]], analysis = x$analysis[1], test = t,
             replicates = nrow(x), predicted_power_formula = mean(x$predicted_power), realised_power = mean(sig),
             realised_l95 = ci[1], realised_u95 = ci[2], formula_within_ci = mean(x$predicted_power) >= ci[1] & mean(x$predicted_power) <= ci[2],
             mean_chisq = mean(chi), median_chisq = median(chi), mean_z = mean(z), sd_z = sd(z),
             n_significant = sum(sig), n_significant_not_by_wald = sum(sig & x$chi_wald <= crit),
             n_wald_significant_not_by_this = sum(!sig & x$chi_wald > crit)) }))))
pow <- pow[order(pow$locus, pow$analysis, match(pow$test, TESTS)), ]

## ---- null calibration, European scan ---------------------------------------------------------
bim <- read.table(paste0(IN, ".bim"), col.names = c("CHR", "ID", "CM", "POS", "A1", "A2"), stringsAsFactors = FALSE)
cpos <- bim[match(CAUSAL, bim$ID), ]
near <- Reduce(`|`, lapply(seq_len(nrow(cpos)), function(j) bim$CHR == cpos$CHR[j] & abs(bim$POS - cpos$POS[j]) <= NULL_DIST))
cand <- intersect(readLines(PRUNE), bim$ID[!near & bim$CHR %in% 1:22])
writeLines(cand, file.path(tmp, "cand"));
system2("plink2", c("--bfile", IN, "--keep", EUR_KEEP, "--extract", file.path(tmp, "cand"), "--freq", "--out", file.path(tmp, "f")), stdout = FALSE, stderr = FALSE)
fr <- rd(file.path(tmp, "f.afreq")); cand <- sort(fr$ID[pmin(fr$ALT_FREQS, 1 - fr$ALT_FREQS) >= NULL_MAF])
set.seed(NULL_SEED); nullv <- sort(sample(cand, N_NULL))
writeLines(nullv, file.path(tmp, "nullv"))
if (system2("plink2", c("--bfile", IN, "--keep", EUR_KEEP, "--extract", file.path(tmp, "nullv"), "--export", "A", "--out", file.path(tmp, "n")),
            stdout = FALSE, stderr = FALSE) != 0) stop("plink2 export failed")
nr <- read.table(file.path(tmp, "n.raw"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
Gn <- as.matrix(nr[match(eur_ids, nr$IID), -(1:6)])
null_one <- function(f) {
  ph <- read.delim(f); age <- read.delim(sub("\\.pheno$", ".age", f)); dd <- design$eur_4a
  y <- as.integer(ph$PHENO[match(dd$ids, ph$IID)] == 2); X0 <- cbind(1, as.matrix(dd$cov[, -1]), AGE = age$AGE[match(dd$ids, age$IID)])
  t(sapply(seq_len(ncol(Gn)), function(j) { g <- Gn[, j]; ok <- !is.na(g)
    tests(cbind(g[ok], X0[ok, , drop = FALSE]), y[ok], do_firth = FALSE)[c("chi_wald", "chi_score", "chi_lrt")] }))
}
nres <- mclapply(reps, function(f) tryCatch(null_one(f), error = function(e) e), mc.cores = N_PAR, mc.preschedule = FALSE)
if (any(sapply(nres, function(x) inherits(x, "error") || inherits(x, "try-error")))) stop("null calibration failed: ", as.character(nres[[which(sapply(nres, function(x) !is.matrix(x)))[1]]]))
lam <- function(chi) median(chi, na.rm = TRUE) / qchisq(0.5, 1)
tails <- c(0.05, 0.01, 0.001, 1e-4)
nullt <- do.call(rbind, lapply(c("wald", "score", "lrt"), function(t) { cn <- paste0("chi_", t)
  per_rep <- sapply(nres, function(m) lam(m[, cn])); all <- unlist(lapply(nres, function(m) m[, cn])); p <- pchisq(all, 1, lower.tail = FALSE)
  o <- data.frame(test = t, variants = N_NULL, replicates = length(nres), tests = length(all), lambda_pooled = lam(all),
                  lambda_mean_over_replicates = mean(per_rep), lambda_sd_over_replicates = sd(per_rep),
                  lambda_se_of_mean = sd(per_rep) / sqrt(length(per_rep)), mean_chisq = mean(all))
  for (a in tails) o[[paste0("ratio_observed_expected_p_lt_", format(a, scientific = FALSE))]] <- mean(p < a) / a
  o }))

num <- sapply(pow, is.numeric); pow[num] <- lapply(pow[num], function(v) signif(v, 4))
num <- sapply(nullt, is.numeric); nullt[num] <- lapply(nullt[num], function(v) signif(v, 5))
dw <- d; num <- sapply(dw, is.numeric); dw[num] <- lapply(dw[num], function(v) signif(v, 6))
wr <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
wr(dw, "pdac_demo_08_tests_per_replicate.tsv"); wr(pow, "pdac_demo_08_tests_power.tsv"); wr(nullt, "pdac_demo_08_tests_null.tsv")

R <- c("SECTION 8: FOUR TESTS AT THE CAUSAL VARIANTS OVER THE REPLICATES",
       sprintf("%d replicates; European scan with SEX, AGE, PC1-PC%d; meta-analysis of the three strata; P < %g", length(reps), k, ALPHA),
       chk_line, firth_line, "",
       "REALISED POWER (95% CI) AGAINST THE SECTION 3B FORMULA; mean and SD of the signed statistic z")
for (j in seq_len(nrow(pow))) with(pow[j, ], R <<- c(R, sprintf("  %-13s OR %.2f  %-7s %-5s formula %5.1f%%  realised %5.1f%% (%.1f to %.1f)  z %.2f, %.2f%s",
  label, or_generative, analysis, test, 100 * predicted_power_formula, 100 * realised_power, 100 * realised_l95, 100 * realised_u95, mean_z, sd_z,
  if (test != "wald") sprintf("  significant where Wald is not: %d; the reverse: %d", n_significant_not_by_wald, n_wald_significant_not_by_this) else "")))
R <- c(R, "", sprintf("NULL CALIBRATION, European scan: %d pruned variants (MAF >= %g, more than %g Mb from the causal variants) x %d replicates",
                      N_NULL, NULL_MAF, NULL_DIST / 1e6, length(nres)))
for (j in seq_len(nrow(nullt))) R <- c(R, sprintf("  %-5s lambda %.4f pooled; mean over replicates %.4f (SD %.4f, SE %.4f); observed / expected at P < 0.05, 0.01, 0.001, 0.0001: %s",
  nullt$test[j], nullt$lambda_pooled[j], nullt$lambda_mean_over_replicates[j], nullt$lambda_sd_over_replicates[j], nullt$lambda_se_of_mean[j],
  paste(sprintf("%.3f", unlist(nullt[j, grep("^ratio_", names(nullt))])), collapse = ", ")))
writeLines(R, file.path(OUT, "pdac_demo_08_tests_summary.txt")); cat(R, sep = "\n")
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
