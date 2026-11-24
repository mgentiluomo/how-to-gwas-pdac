#!/usr/bin/env Rscript

################################################################################
# Section 8, step 04: why realised power and effect sizes differ from the
# prediction, and why TERT/CLPTM1L falls short in the European set
#
# 1. STANDARD ERROR AND POWER, European scan, per causal locus.
#    (a) Section 3B formula: information N_eff / 4 * 2 p (1 - p), the
#        information at the null with every individual at the average risk.
#    (b) Expected information of the JOINT generative model (the three loci,
#        sex, age): I = X' W X, W = diag(pi_i (1 - pi_i)), pi_i generative.
#    (c) Expected information of the model the scan actually FITS: one locus,
#        SEX, AGE and PC1..PCk, the other causal loci omitted. Its coefficient
#        beta* is the population limit of that misspecified model, found by
#        fitting it to the generative probabilities of the same individuals
#        (quasi-binomial, fractional response); its standard error is the
#        model-based one PLINK reports, from W* = diag(pi*_i (1 - pi*_i)).
#        beta* differs from the generative log OR through non-collapsibility
#        and through chance correlation with the omitted loci.
#    Power at P < 5e-8 from each, and the realised power of step 02.
# 2. GENOTYPE CORRELATION between the risk-allele counts of the three causal
#    variants within each analysis stratum (non-missing genotypes).
# 3. TERT/CLPTM1L ADJUSTED FOR ABO: in every replicate and in the released
#    phenotype, the European TERT estimate with SEX, AGE and PC1..PCk, then
#    also adjusted for the ABO causal variant, then for ABO and 13q22.1.
#
# Reuses the generator of step 01 (same seeds, same code); reads step 02.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_expected_se.tsv     per locus: SE and power, (a) (b) (c), realised
#   pdac_demo_08_genotype_r.tsv      per stratum and pair: n and r
#   pdac_demo_08_tert_adjusted.tsv   TERT estimate / generative OR by model
#   pdac_demo_08_expected_se.txt     plain-text report of the three
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
here <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()

# load the generator of step 01 (everything before its release check)
src <- readLines(file.path(here, "01_simulate_replicates.R"))
src <- sub("^script_dir <- .*", sprintf("script_dir <- \"%s\"", here), src)
stop_at <- grep("^## check against the release", src)[1] - 1
eval(parse(text = src[1:stop_at]))   # sets the working directory to the project root

OUT <- "sections/08_replicates/results"
ALPHA <- 5e-8
loci <- c("G1", "G2", "G3"); LABEL <- c(G1 = "ABO", G2 = "TERT/CLPTM1L", G3 = "13q22.1")
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
keep <- read.table("results/pca/pdac_demo_02_eur_keep.txt", stringsAsFactors = FALSE)[[2]]
k <- as.integer(read.table("results/assoc/pdac_demo_04A_npc_chosen.txt")[1, 2])
cv <- rd("results/assoc/pdac_demo_04A_covar.txt")
idx <- IID %in% keep
PC <- as.matrix(cv[match(IID[idx], cv$IID), paste0("PC", seq_len(k)), drop = FALSE])
if (anyNA(PC)) stop("components missing for European individuals")
per <- read.delim(file.path(OUT, "pdac_demo_08_per_replicate.tsv"), stringsAsFactors = FALSE)
reps <- sort(unique(per$replicate))
Ge <- G[idx, loci]; male_e <- male[idx]

## ---- 1. standard errors ------------------------------------------------------------------
se_of <- function(X, p) sqrt(diag(solve(crossprod(X * sqrt(p * (1 - p))))))
res <- lapply(reps, function(i) {
  s <- simulate(SEED_BASE + i, loci)
  eta <- log(OR_MALE) * male + log(OR_AGE) * s$age_z
  for (l in loci) eta <- eta + log(OR[[l]]) * G[, l]
  e <- anc == "eur"
  a <- uniroot(function(a) mean(plogis(a + eta[e])) - TARGET_CASE_FRACTION[["eur"]],
               interval = c(-25, 25), tol = 1e-10)$root
  p <- plogis(a + eta[idx]); age_e <- s$age_z[idx]
  joint <- se_of(cbind(1, Ge, male_e, age_e), p)[2:4]
  fit <- t(sapply(loci, function(l) {
    X <- cbind(1, g = Ge[, l], male_e, age_e, PC)
    m <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
    c(beta = unname(m$coefficients[2]), se = unname(se_of(X, m$fitted.values)[2]))
  }))
  y <- s$status[idx]; neff <- 4 / (1 / sum(y == 1) + 1 / sum(y == 0))
  q <- colMeans(Ge) / 2
  list(joint = joint, fit_beta = fit[, "beta"], fit_se = fit[, "se"], formula = 1 / sqrt(neff / 4 * 2 * q * (1 - q)),
       y = y, age = s$age[idx])
})
pick <- function(nm, l) mean(sapply(res, function(r) r[[nm]][[l]]))
pow <- function(b, se) pchisq(qchisq(1 - ALPHA, 1), 1, ncp = (b / se)^2, lower.tail = FALSE)
tab <- do.call(rbind, lapply(loci, function(l) {
  x <- per[per$locus == l & per$analysis == "eur_4a", ]
  b0 <- log(OR[[l]]); se_f <- pick("formula", l); se_j <- pick("joint", l)
  b_c <- pick("fit_beta", l); se_c <- pick("fit_se", l)
  z <- x$log_or_hat / x$se
  data.frame(locus = l, label = LABEL[[l]], or_generative = OR[[l]], n_replicates = nrow(x), n_pcs = k,
             se_formula = se_f, se_expected_joint = se_j, se_expected_fitted = se_c,
             or_expected_fitted = exp(b_c), se_realised_mean = mean(x$se), sd_log_or_realised = sd(x$log_or_hat),
             mean_log_or_realised = mean(x$log_or_hat),
             power_formula = pow(b0, se_f), power_expected_joint = pow(b0, se_j), power_expected_fitted = pow(b_c, se_c),
             power_realised = mean(x$significant), z_realised_mean = mean(z), z_realised_sd = sd(z),
             cor_est_se = cor(x$log_or_hat, x$se),
             row.names = NULL)
}))
write.table(tab, file.path(OUT, "pdac_demo_08_expected_se.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 2. genotype correlation within strata --------------------------------------------
Graw <- sapply(CAUSAL, function(id) raw[[col_of(id)]])   # counted allele = risk allele, NA kept
strata <- list(eur = "results/pca/pdac_demo_02_eur_keep.txt", afr = "results/meta/afr/keep.txt", eas = "results/meta/eas/keep.txt")
pairs <- list(c("G1", "G2"), c("G1", "G3"), c("G2", "G3"))
rtab <- do.call(rbind, lapply(names(strata), function(g) {
  ids <- read.table(strata[[g]], stringsAsFactors = FALSE)[[2]]; X <- Graw[IID %in% ids, ]
  do.call(rbind, lapply(pairs, function(pp) { ok <- complete.cases(X[, pp])
    data.frame(stratum = g, pair = paste0(LABEL[[pp[1]]], " x ", LABEL[[pp[2]]]), n_individuals = sum(IID %in% ids),
               n_complete = sum(ok), r = cor(X[ok, pp[1]], X[ok, pp[2]])) }))
}))
write.table(rtab, file.path(OUT, "pdac_demo_08_genotype_r.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 3. TERT/CLPTM1L adjusted for ABO ---------------------------------------------------
models <- list(unadjusted = character(), adjusted_abo = "G1", adjusted_abo_13q = c("G1", "G3"))
tert_fit <- function(y, age, extra) {
  d <- data.frame(y = y, g = Ge[, "G2"], male = male_e, age = age, PC, Ge[, extra, drop = FALSE])
  ok <- complete.cases(d)
  unname(coef(glm(y ~ ., binomial, d[ok, ]))["g"])
}
rel_ph <- read.delim("demo_data/phenotype.txt"); rel_cv <- read.delim("demo_data/covariates.txt")
y_rel <- as.integer(rel_ph$PHENO[match(IID[idx], rel_ph$IID)] == 2); age_rel <- rel_cv$AGE[match(IID[idx], rel_cv$IID)]
ttab <- do.call(rbind, lapply(names(models), function(m) {
  b <- sapply(res, function(r) tert_fit(r$y, r$age, models[[m]]))
  b_rel <- tert_fit(y_rel, age_rel, models[[m]])
  rbind(data.frame(data = "replicates", model = m, n = length(b), mean_or_ratio = mean(exp(b) / OR[["G2"]]),
                   mean_log_or = mean(b), log_or_generative = log(OR[["G2"]])),
        data.frame(data = "released_phenotype", model = m, n = 1, mean_or_ratio = exp(b_rel) / OR[["G2"]],
                   mean_log_or = b_rel, log_or_generative = log(OR[["G2"]])))
}))
write.table(ttab, file.path(OUT, "pdac_demo_08_tert_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- report --------------------------------------------------------------------------
R <- c("SECTION 8: STANDARD ERROR, POWER AND BIAS OF THE CAUSAL ESTIMATES, EUROPEAN SCAN",
       sprintf("%d replicates; European set %d individuals; scan model SEX, AGE, PC1-PC%d; P < %g", length(reps), sum(idx), k, ALPHA), "",
       "1. STANDARD ERROR OF THE LOG OR AND POWER",
       "   (a) Section 3B formula  (b) expected information, joint generative model",
       "   (c) expected information, model fitted by the scan (one locus, other loci omitted)", "")
for (j in seq_len(nrow(tab))) with(tab[j, ], R <<- c(R,
  sprintf("  %-13s OR %.2f  SE (a) %.4f  (b) %.4f  (c) %.4f  realised %.4f (SD of log OR %.4f)",
          label, or_generative, se_formula, se_expected_joint, se_expected_fitted, se_realised_mean, sd_log_or_realised),
  sprintf("  %-13s power (a) %5.1f%%  (b) %5.1f%%  (c) %5.1f%%  realised %5.1f%%;  OR of the fitted model (c) %.3f, realised mean log OR %.3f (generative %.3f); realised z mean %.2f, SD %.2f",
          "", 100 * power_formula, 100 * power_expected_joint, 100 * power_expected_fitted, 100 * power_realised,
          or_expected_fitted, mean_log_or_realised, log(or_generative), z_realised_mean, z_realised_sd)))
R <- c(R, "", "2. CORRELATION OF RISK-ALLELE COUNTS WITHIN STRATA (non-missing genotypes)")
for (j in seq_len(nrow(rtab))) with(rtab[j, ], R <<- c(R, sprintf("  %-4s %-28s n %4d  r %+.3f", stratum, pair, n_complete, r)))
R <- c(R, "", "3. TERT/CLPTM1L, EUROPEAN SET: estimated / generative OR (1.50)")
for (j in seq_len(nrow(ttab))) with(ttab[j, ], R <<- c(R, sprintf("  %-19s %-17s OR ratio %.3f  log OR %.3f (generative %.3f)",
  data, model, mean_or_ratio, mean_log_or, log_or_generative)))
R <- c(R, "", "HOW TO READ THIS",
  "  (a) is the information at the null with every individual at the average",
  "  risk. Under the generator risk varies with the loci, sex and age, and",
  "  individuals far from risk one half carry less information: (b) and (c) are",
  "  larger than (a). (c) is the standard error of the model the scan fits; its",
  "  coefficient is the limit of that model, which omits the other causal loci",
  "  and so differs from the generative log OR (non-collapsibility, and chance",
  "  correlation with the omitted loci, section 2). The realised z has SD below",
  "  1 mainly because the Wald standard error of a logistic coefficient grows",
  "  with the estimate (correlation of estimate and SE below), so z varies less",
  "  than the estimate; power computed as if z ~ N(mean, 1) is then too low",
  "  when the mean is above the threshold and too high when it is below.",
  "  Section 3: a negative correlation between the TERT and ABO risk alleles",
  "  in the European set biases the single-locus TERT estimate toward the null.",
  "  Over the replicates, adjusting for ABO and 13q22.1 removes the shortfall;",
  "  in the released phenotype it removes part of it, the rest being the",
  "  sampling error of one draw.")
for (j in seq_len(nrow(tab))) with(tab[j, ], R <<- c(R, sprintf("  %-13s correlation of estimate and SE over replicates %.3f", label, cor_est_se)))
writeLines(R, file.path(OUT, "pdac_demo_08_expected_se.txt"))
cat(R, sep = "\n")
