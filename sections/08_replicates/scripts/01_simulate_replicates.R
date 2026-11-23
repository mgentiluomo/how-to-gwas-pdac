#!/usr/bin/env Rscript

################################################################################
# Section 8, step 01: simulate replicate phenotypes
#
# PURPOSE
#   The released dataset is one draw. A single draw can show that a result is
#   possible; it cannot show how often it happens. This step draws N_REP new
#   case-control phenotypes on the same genotypes, with the same generator as
#   the release (demo_dataset/scripts/04_simulate_phenotype.R) and a third
#   causal locus, so that power, winner's curse, credible-set coverage and the
#   behaviour of I-squared and lambda can be measured over replicates.
#   The released dataset does not change.
#
# GENERATIVE MODEL (logistic, per copy of the counted allele, as the release)
#   G1  9:133273682:A:T   ABO 9q34.2              OR 2.60   (as the release)
#   G2  5:1286401:C:A     TERT/CLPTM1L 5p15.33    OR 1.50   (as the release)
#   G3  13:73362027:C:T   13q22.1                 OR 1.97   (new)
#   male sex OR 1.30; age OR 1.50 per 10 years; intercept solved per ancestry
#   group for the case fractions of the release (eur 0.333, afr 0.456,
#   eas 0.414); age drawn as in the release.
#
#   G3 was chosen by a rule fixed before it was simulated: the array variant
#   closest to rs9543325 (13:73,342,491, GRCh38; the 13q22.1 PDAC locus, a
#   stand-in chosen for familiarity like ABO and TERT) that is biallelic, not
#   palindromic, with MAF 0.25 to 0.5 in the European analysis set and at
#   least 0.10 in the African and East Asian strata. Its MAF there: 0.364,
#   0.479, 0.328. OR 1.97 gives a predicted power of 50.7% in the European
#   analysis set (N_eff 562, P < 5e-8), with the power formula of Section 3B.
#
# DIFFERENCES FROM THE RELEASE GENERATOR, AND WHY
#   - Three causal loci instead of two.
#   - Seeds: replicate i uses SEED_BASE + i (SEED_BASE 20261124).
#   - The release generator refuses to write data whose re-estimated
#     parameters fall outside their confidence intervals. Applied to
#     replicates that check would discard about one draw in five at random
#     and condition the experiment on the realised effects being close to the
#     generative ones, which is exactly what the experiment measures. It is
#     therefore not applied; the realised odds ratios on all 1,461 individuals
#     are recorded for every replicate instead.
#   - Survival times are not drawn.
#   CHECK: with seed 2026 and the two release loci, this code reproduces the
#   released phenotype.txt and AGE in covariates.txt exactly; the step stops
#   otherwise, because the replicates are only comparable to the release if
#   they come from the same generator.
#
# OUTPUT (sections/08_replicates/results/sim/, not committed: regenerated)
#   rep_NNN.pheno   FID IID PHENO
#   rep_NNN.age     IID AGE
#   and sections/08_replicates/results/pdac_demo_08_realised.tsv (committed)
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "..", "..", ".."))
setwd(root)

N_REP     <- as.integer(Sys.getenv("GWAS_08_N_REP", "200"))
SEED_BASE <- 20261124L
OUT       <- "sections/08_replicates/results"
SIM       <- file.path(OUT, "sim")
dir.create(SIM, showWarnings = FALSE, recursive = TRUE)

CAUSAL <- c(G1 = "9:133273682:A:T", G2 = "5:1286401:C:A", G3 = "13:73362027:C:T")
OR     <- c(G1 = 2.60, G2 = 1.50, G3 = 1.97)
OR_MALE <- 1.30; OR_AGE <- 1.50
TARGET_CASE_FRACTION <- c(eur = 0.333, afr = 0.456, eas = 0.414)
AGE_MEAN <- 66.7; AGE_SD <- 10.0; AGE_MIN <- 33; AGE_MAX <- 90

## genotypes of the causal variants, extracted exactly as the release generator does
tmp <- tempfile("sim"); dir.create(tmp)
writeLines(unname(CAUSAL), file.path(tmp, "snps"))
if (system2("plink2", c("--bfile", "demo_data/pdac_demo", "--extract", file.path(tmp, "snps"),
                        "--recode", "A", "--out", file.path(tmp, "causal")), stdout = FALSE, stderr = FALSE) != 0)
  stop("plink2 extraction failed")
raw <- read.table(file.path(tmp, "causal.raw"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
col_of <- function(id) { h <- grep(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", id), "_"), names(raw), value = TRUE)
  if (length(h) != 1) stop("no unique column for ", id); h }
G <- sapply(CAUSAL, function(id) { x <- raw[[col_of(id)]]; x[is.na(x)] <- mean(x, na.rm = TRUE); x })
IID <- raw$IID; n <- length(IID)
anc <- read.table("demo_data/sample_ancestry.tsv", col.names = c("IID", "group"), stringsAsFactors = FALSE)
anc <- anc$group[match(IID, anc$IID)]
fam <- read.table("demo_data/pdac_demo.fam", stringsAsFactors = FALSE)
male <- as.integer(fam$V5[match(IID, fam$V2)] == 1)

simulate <- function(seed, loci) {
  set.seed(seed)
  age <- round(pmin(pmax(rnorm(n, AGE_MEAN, AGE_SD), AGE_MIN), AGE_MAX))
  age_z <- (age - AGE_MEAN) / 10
  # summed term by term, in the order of the release generator, so that the
  # release check is exact to the last bit
  eta <- log(OR[[loci[1]]]) * G[, loci[1]]
  for (l in loci[-1]) eta <- eta + log(OR[[l]]) * G[, l]
  eta <- eta + log(OR_MALE) * male + log(OR_AGE) * age_z
  p <- numeric(n)
  for (g in names(TARGET_CASE_FRACTION)) {
    idx <- which(anc == g)
    a <- uniroot(function(a) mean(plogis(a + eta[idx])) - TARGET_CASE_FRACTION[[g]],
                 interval = c(-25, 25), tol = 1e-10)$root
    p[idx] <- plogis(a + eta[idx])
  }
  status <- rbinom(n, 1, p)
  list(age = age, age_z = age_z, status = status)
}

## check against the release
chk <- simulate(2026L, c("G1", "G2"))
rel_ph <- read.delim("demo_data/phenotype.txt"); rel_cv <- read.delim("demo_data/covariates.txt")
same_ph <- identical(as.integer(ifelse(chk$status == 1, 2L, 1L)), as.integer(rel_ph$PHENO[match(IID, rel_ph$IID)]))
same_age <- identical(as.numeric(chk$age), as.numeric(rel_cv$AGE[match(IID, rel_cv$IID)]))
cat(sprintf("release check, seed 2026, two loci: phenotype identical %s, age identical %s\n", same_ph, same_age))
if (!same_ph || !same_age) stop("this code does not reproduce the release generator; replicates would not be comparable")

## replicates
realised <- vector("list", N_REP)
for (i in seq_len(N_REP)) {
  s <- simulate(SEED_BASE + i, names(CAUSAL))
  tag <- sprintf("rep_%03d", i)
  write.table(data.frame(FID = IID, IID = IID, PHENO = ifelse(s$status == 1, 2L, 1L)),
              file.path(SIM, paste0(tag, ".pheno")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(IID = IID, AGE = s$age), file.path(SIM, paste0(tag, ".age")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  fit <- glm(s$status ~ G[, "G1"] + G[, "G2"] + G[, "G3"] + male + s$age_z, family = binomial())
  b <- coef(fit)
  realised[[i]] <- data.frame(replicate = i, seed = SEED_BASE + i, cases = sum(s$status), controls = n - sum(s$status),
                              realised_or_G1 = exp(b[2]), realised_or_G2 = exp(b[3]), realised_or_G3 = exp(b[4]))
}
realised <- do.call(rbind, realised)
write.table(realised, file.path(OUT, "pdac_demo_08_realised.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d replicates to %s (seeds %d to %d)\n", N_REP, SIM, SEED_BASE + 1L, SEED_BASE + N_REP))
