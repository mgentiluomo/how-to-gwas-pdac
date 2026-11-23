#!/usr/bin/env Rscript

################################################################################
# Section 8, step 04: why realised power falls short of predicted power
#
# The predicted power of step 02 uses the formula of Section 3B: information
# per allele N_eff / 4 * 2 p (1 - p), that is N * phi * (1 - phi) * 2 p (1 - p)
# with phi the case fraction. That is the information at the null. Under the
# generator the risk of each individual varies with the causal loci, sex and
# age, and the information about one locus is
#     I = X' W X,  W = diag(pi_i (1 - pi_i)),
# with pi_i the generative probability and X the design (intercept, the three
# loci, sex, age). This step computes, for every replicate, the standard
# error of each causal log OR from that expected information, the standard
# error of the Section 3B formula, and the power each implies at P < 5e-8 with
# the generative OR. Both use the European analysis set of the canonical run.
#
# It reuses the generator of step 01 (same seeds, same code) and reads the
# realised standard errors of step 02 for comparison.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_expected_se.tsv       per locus: mean SE (formula, expected
#                                      information, realised) and power
#   pdac_demo_08_expected_se.txt       plain-text report
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
keep <- read.table("results/pca/pdac_demo_02_eur_keep.txt", stringsAsFactors = FALSE)[[2]]
idx <- IID %in% keep
loci <- c("G1", "G2", "G3")
per <- read.delim(file.path(OUT, "pdac_demo_08_per_replicate.tsv"), stringsAsFactors = FALSE)
reps <- sort(unique(per$replicate))

res <- t(sapply(reps, function(i) {
  s <- simulate(SEED_BASE + i, loci)
  eta <- log(OR_MALE) * male + log(OR_AGE) * s$age_z
  for (l in loci) eta <- eta + log(OR[[l]]) * G[, l]
  e <- anc == "eur"
  a <- uniroot(function(a) mean(plogis(a + eta[e])) - TARGET_CASE_FRACTION[["eur"]],
               interval = c(-25, 25), tol = 1e-10)$root
  p <- plogis(a + eta[idx])
  X <- cbind(1, G[idx, loci], male[idx], s$age_z[idx])
  se_info <- sqrt(diag(solve(crossprod(X * sqrt(p * (1 - p))))))[2:4]
  y <- s$status[idx]; neff <- 4 / (1 / sum(y == 1) + 1 / sum(y == 0))
  q <- colMeans(G[idx, loci]) / 2
  se_formula <- 1 / sqrt(neff / 4 * 2 * q * (1 - q))
  c(se_info, se_formula)
}))
colnames(res) <- c(paste0("info_", loci), paste0("formula_", loci))

pow <- function(se, or) pchisq(qchisq(1 - ALPHA, 1), 1, ncp = (log(or) / se)^2, lower.tail = FALSE)
tab <- do.call(rbind, lapply(loci, function(l) {
  x <- per[per$locus == l & per$analysis == "eur_4a", ]
  se_f <- mean(res[, paste0("formula_", l)]); se_i <- mean(res[, paste0("info_", l)])
  data.frame(locus = l, or_generative = OR[[l]], n_replicates = nrow(x),
             se_formula = se_f, se_expected_information = se_i,
             se_realised_mean = mean(x$se), sd_log_or_realised = sd(x$log_or_hat),
             power_formula = pow(se_f, OR[[l]]), power_expected_information = pow(se_i, OR[[l]]),
             z_realised_mean = mean(x$log_or_hat / x$se), z_realised_sd = sd(x$log_or_hat / x$se),
             power_realised = mean(x$significant))
}))
write.table(tab, file.path(OUT, "pdac_demo_08_expected_se.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

R <- c("SECTION 8: STANDARD ERROR OF THE CAUSAL LOG OR, EUROPEAN SCAN",
       sprintf("%d replicates; power at P < %g with the generative OR", length(reps), ALPHA), "",
       sprintf("  %-4s %5s  %-28s %-28s %-14s %s", "", "OR", "Section 3B formula", "expected information", "realised SE", "realised power"))
for (k in seq_len(nrow(tab))) with(tab[k, ], R <<- c(R, sprintf(
  "  %-4s %5.2f  SE %.4f  power %5.1f%%      SE %.4f  power %5.1f%%      %.4f (sd %.4f)  %5.1f%%",
  locus, or_generative, se_formula, 100 * power_formula, se_expected_information,
  100 * power_expected_information, se_realised_mean, sd_log_or_realised, 100 * power_realised)))
R <- c(R, "", "HOW TO READ THIS",
  "  The Section 3B formula gives the information at the null with every",
  "  individual at the average risk. Under the generator, risk varies with the",
  "  three loci, sex and age, and individuals far from a probability of one",
  "  half carry less information. The expected-information standard error",
  "  matches the realised one; the formula's is smaller, so its power is",
  "  optimistic. The expected-information power assumes z ~ N(log OR / SE, 1).",
  "  Across replicates the realised z has the mean and SD below; an SD under 1",
  "  puts less mass above the threshold when the mean is below it.")
for (k in seq_len(nrow(tab))) with(tab[k, ], R <<- c(R, sprintf("  %-4s realised z mean %.2f, SD %.2f", locus, z_realised_mean, z_realised_sd)))
writeLines(R, file.path(OUT, "pdac_demo_08_expected_se.txt"))
cat(R, sep = "\n")
