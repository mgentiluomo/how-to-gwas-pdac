#!/usr/bin/env Rscript

################################################################################
# Section 8, step 10: the smallest detectable OR by simulation, at one frequency
#
# PURPOSE
#   Steps 03 and 05 show that the Section 3B formula overstates power at the
#   simulated effects whichever test is used. The detectable effects the
#   formula gives (its curve of the smallest OR detectable at 80% power) are
#   therefore too small. This step measures by how much at one frequency, in
#   the European analysis set, by simulation.
#
# WHAT IT DOES
#   One variant, the 13q22.1 causal variant of the replicates (the frequency
#   of its minor allele in the European analysis set is the frequency of the
#   experiment), on the real genotypes of the European analysis set, in two
#   orientations: the risk allele is the minor allele, or the major allele.
#   The formula is symmetric in the allele frequency and cannot tell them
#   apart; in the replicates the risk allele of ABO and of 13q22.1 is the
#   major one (frequency 0.62 and 0.64 in the European set). For each
#   orientation and each OR of a grid, N_SIM phenotypes are drawn from a
#   logistic model with that variant as the only causal locus, male sex OR 1.30 and age OR 1.50 per 10
#   years as in the generator, age drawn as in the generator, and the
#   intercept solved so that the expected case fraction is that of the
#   canonical analysis set (so the expected N_eff is the one the formula is
#   given). Each phenotype is tested with the model of the scan (SEX, AGE,
#   PC1..PCk) by the Wald, score and likelihood-ratio tests. Power at
#   P < 5e-8 per OR and test; the OR at 80% power by linear interpolation on
#   the logit of power between the two grid values that bracket it.
#   The formula value is computed for the mean realised N_eff and the same
#   frequency, with the expression of Section 3B (non-centrality
#   N_eff / 4 * 2 p (1 - p) * log(OR)^2).
#   Limits: one variant and one frequency; one causal locus (the replicates
#   have three); the European scan only.
#
# NOT PART OF THE CANONICAL RUN: reads its outputs in results/.
# SETTINGS: N_SIM 2,000 per OR, seed 20261124; GWAS_08_PARALLEL cores (default 6).
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_simulated_power_grid.tsv     orientation x OR x test: power, 95% CI
#   pdac_demo_08_simulated_detectable_or.tsv  orientation x test: OR at 80% power; the formula value
#   pdac_demo_08_simulated_detectable_or.txt  plain-text report
################################################################################

suppressPackageStartupMessages(library(parallel))
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(normalizePath(file.path(script_dir, "..", "..", "..")))
OUT <- "sections/08_replicates/results"
VARIANT <- "13:73362027:C:T"; LOCUS <- "13q22.1"
OR_GRID <- seq(1.6, 3.6, by = 0.1); ORIENT <- c("minor", "major"); N_SIM <- 2000L; SEED <- 20261124L
OR_MALE <- 1.30; OR_AGE <- 1.50; AGE_MEAN <- 66.7; AGE_SD <- 10.0; AGE_MIN <- 33; AGE_MAX <- 90   # as step 01
ALPHA <- 5e-8; TARGET <- 0.80; crit <- qchisq(1 - ALPHA, 1)
N_PAR <- as.integer(Sys.getenv("GWAS_08_PARALLEL", "6"))
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
k <- as.integer(read.table("results/assoc/pdac_demo_04A_npc_chosen.txt")[1, 2])
keep <- read.table("results/pca/pdac_demo_02_eur_keep.txt", stringsAsFactors = FALSE)[[2]]
cv <- rd("results/assoc/pdac_demo_04A_covar.txt"); cv <- cv[match(keep, cv$IID), ]
ph <- read.delim("demo_data/phenotype.txt"); st <- ph$PHENO[match(keep, ph$IID)]
case_fraction <- mean(st == 2); neff_canonical <- 4 / (1 / sum(st == 2) + 1 / sum(st == 1))
tmp <- tempfile("sd"); dir.create(tmp); writeLines(VARIANT, file.path(tmp, "v"))
if (system2("plink2", c("--bfile", "results/pca/pdac_demo_02_hwe_filt", "--extract", file.path(tmp, "v"), "--export", "A",
                        "--out", file.path(tmp, "g")), stdout = FALSE, stderr = FALSE) != 0) stop("plink2 export failed")
raw <- read.table(file.path(tmp, "g.raw"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
g <- raw[match(keep, raw$IID), 7]; ok <- !is.na(g); g <- g[ok]
if (mean(g) / 2 > 0.5) g <- 2 - g          # g counts the minor allele
maf <- mean(g) / 2; n <- length(g)
male <- as.integer(cv$SEX[ok] == 1); PC <- as.matrix(cv[ok, paste0("PC", seq_len(k)), drop = FALSE])

fit_ml <- function(X, y) { m <- suppressWarnings(glm.fit(X, y, family = binomial())); mu <- m$fitted.values
  list(b = m$coefficients, mu = mu, w = mu * (1 - mu), dev = m$deviance) }
one <- function(or, seed, risk) {
  set.seed(seed)
  age <- round(pmin(pmax(rnorm(n, AGE_MEAN, AGE_SD), AGE_MIN), AGE_MAX)); age_z <- (age - AGE_MEAN) / 10
  eta <- log(or) * (if (risk == "minor") g else 2 - g) + log(OR_MALE) * male + log(OR_AGE) * age_z
  a <- uniroot(function(a) mean(plogis(a + eta)) - case_fraction, interval = c(-25, 25), tol = 1e-10)$root
  y <- rbinom(n, 1, plogis(a + eta))
  X0 <- cbind(1, male, age, PC); X <- cbind(g, X0); f1 <- fit_ml(X, y); f0 <- fit_ml(X0, y)
  se <- sqrt(solve(crossprod(X * sqrt(f1$w)))[1, 1])
  U <- sum(g * (y - f0$mu)); A <- crossprod(X0 * f0$w, g)
  V <- sum(f0$w * g^2) - drop(crossprod(A, solve(crossprod(X0 * sqrt(f0$w)), A)))
  c(wald = unname((f1$b[1] / se)^2), score = U^2 / V, lrt = f0$dev - f1$dev, neff = 4 / (1 / sum(y) + 1 / (n - sum(y))))
}
wilson <- function(x, n) { p <- x / n; z <- 1.959964
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
t0 <- Sys.time()
jobs <- expand.grid(j = seq_along(OR_GRID), o = seq_along(ORIENT))
res_all <- mclapply(seq_len(nrow(jobs)), function(q) { j <- jobs$j[q]; o <- jobs$o[q]
  t(sapply(seq_len(N_SIM), function(i) one(OR_GRID[j], SEED + 10000000L * o + 100000L * j + i, ORIENT[o]))) },
  mc.cores = N_PAR, mc.preschedule = FALSE)
if (any(!sapply(res_all, is.matrix))) stop("simulation failed")
neff_mean <- mean(sapply(res_all, function(m) mean(m[, "neff"])))
p_formula <- function(or, N) pchisq(crit, 1, ncp = N / 4 * 2 * maf * (1 - maf) * log(or)^2, lower.tail = FALSE)
TESTS <- c("wald", "score", "lrt")
grid <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(q) do.call(rbind, lapply(TESTS, function(t) { j <- jobs$j[q]
  s <- sum(res_all[[q]][, t] > crit); ci <- wilson(s, N_SIM)
  data.frame(risk_allele = ORIENT[jobs$o[q]], or = OR_GRID[j], test = t, simulations = N_SIM, power = s / N_SIM, l95 = signif(ci[1], 4), u95 = signif(ci[2], 4),
             power_formula = signif(p_formula(OR_GRID[j], neff_mean), 4)) }))))
or_at <- function(pw, target = TARGET) { lo <- max(which(pw < target)); hi <- lo + 1
  if (!is.finite(lo) || hi > length(pw) || pw[hi] < target) stop("the OR grid does not bracket ", target)
  q <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  OR_GRID[lo] + (q(target) - q(pw[lo])) / (q(pw[hi]) - q(pw[lo])) * (OR_GRID[hi] - OR_GRID[lo]) }
or_formula <- uniroot(function(o) p_formula(o, neff_mean) - TARGET, c(1.01, 20), tol = 1e-8)$root
det <- do.call(rbind, lapply(ORIENT, function(ra) do.call(rbind, lapply(TESTS, function(t) { x <- grid[grid$test == t & grid$risk_allele == ra, ]
  o <- or_at(x$power); ol <- or_at(x$u95); ou <- or_at(x$l95)
  data.frame(risk_allele = ra, test = t, locus = LOCUS, variant = VARIANT, maf = signif(maf, 4), n = n, case_fraction = signif(case_fraction, 4),
             n_eff_mean = signif(neff_mean, 5), n_eff_canonical = signif(neff_canonical, 5), simulations_per_or = N_SIM,
             or_at_80_power = signif(o, 4), or_at_80_power_l95 = signif(ol, 4), or_at_80_power_u95 = signif(ou, 4),
             or_at_80_power_formula = signif(or_formula, 4), ratio_log_or_simulated_to_formula = signif(log(o) / log(or_formula), 4),
             power_formula_at_simulated_or = signif(p_formula(o, neff_mean), 4)) }))))
wr <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
wr(grid, "pdac_demo_08_simulated_power_grid.tsv"); wr(det, "pdac_demo_08_simulated_detectable_or.tsv")
R <- c("SECTION 8: SMALLEST OR DETECTABLE AT 80% POWER, BY SIMULATION, AT ONE FREQUENCY",
       sprintf("variant %s (%s), minor allele frequency %.3f; European analysis set, %d individuals with a genotype; case fraction %.3f, mean N_eff %.1f (canonical set %.1f)",
               VARIANT, LOCUS, maf, n, case_fraction, neff_mean, neff_canonical),
       sprintf("model fitted: SEX, AGE, PC1-PC%d; %d simulations per OR, OR %.1f to %.1f; P < %g", k, N_SIM, min(OR_GRID), max(OR_GRID), ALPHA), "",
       sprintf("  Section 3B formula: OR %.3f", or_formula))
for (j in seq_len(nrow(det))) with(det[j, ], R <<- c(R, sprintf("  risk allele %s, %-5s simulated: OR %.3f (%.3f to %.3f); log OR %.0f%% above the formula; the formula gives %.1f%% power at that OR",
  risk_allele, test, or_at_80_power, or_at_80_power_l95, or_at_80_power_u95, 100 * (ratio_log_or_simulated_to_formula - 1), 100 * power_formula_at_simulated_or)))
R <- c(R, "", "  power by OR (formula; risk allele minor: wald, score, lrt; risk allele major: wald, score, lrt):")
for (o in OR_GRID) { x <- grid[abs(grid$or - o) < 1e-9, ]; pv <- function(ra, t) 100 * x$power[x$test == t & x$risk_allele == ra]
  R <- c(R, sprintf("  OR %.1f  %5.1f%%;  %5.1f%%  %5.1f%%  %5.1f%%;  %5.1f%%  %5.1f%%  %5.1f%%", o, 100 * x$power_formula[1],
                    pv("minor", "wald"), pv("minor", "score"), pv("minor", "lrt"), pv("major", "wald"), pv("major", "score"), pv("major", "lrt"))) }
writeLines(R, file.path(OUT, "pdac_demo_08_simulated_detectable_or.txt")); cat(R, sep = "\n")
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
