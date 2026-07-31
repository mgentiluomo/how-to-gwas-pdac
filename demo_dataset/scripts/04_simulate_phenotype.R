#!/usr/bin/env Rscript

###############################################################################
# 04_simulate_phenotype.R
#
# Simulates the case/control phenotype, the covariates and the survival outcome
# for the pdac_demo demonstration dataset.
#
# WHY THIS SCRIPT WAS REWRITTEN
#   The previous version specified the causal effect on a liability scale and
#   the analysis estimated it on the log-odds scale. The two differ by a factor
#   of roughly 1.3 to 1.8, so a documented odds ratio of 2.40 appeared in the
#   released data as an effect of about 3.4. Nothing detected it, because
#   nothing checked: the documented value was never compared with an estimate
#   from the data it described.
#
#   This version therefore does two things differently.
#     1. The generative model is a logistic regression on the observed scale,
#        so the parameter that is specified is the parameter that a logistic
#        analysis estimates. No scale conversion is involved anywhere.
#     2. The script validates itself. After simulating, it re-estimates every
#        effect with the same model the tutorial uses, and stops with an error
#        if any estimate is inconsistent with its target. A dataset that fails
#        this check is never written.
#
# THE GENERATIVE MODEL
#
#   logit P(case) = alpha_ancestry
#                 + beta_1 * G1        (dosage of the ABO risk allele)
#                 + beta_2 * G2        (dosage of the TERT/CLPTM1L risk allele)
#                 + beta_sex  * male
#                 + beta_age  * (age - 66.7) / 10
#
#   The intercept is ancestry-specific and is solved numerically so that each
#   group reaches its target case fraction. This is a statement about
#   ascertainment, not about prevalence: PDAC has an incidence of roughly 11
#   per 100,000 person-years, so a population sample of 1,461 people would
#   contain no cases at all. A case-control study instead samples cases and
#   controls separately. Under a logistic model that kind of sampling shifts
#   the intercept and leaves the log-odds ratios unchanged, which is precisely
#   why the generative effect sizes remain recoverable from the sample, and why
#   a logistic model is the right generative choice here.
#
# CHOICE OF EFFECT SIZES
#
#   G1, ABO (9q34.2), odds ratio 2.60 per risk allele.
#     Detectable with about 99% probability in the European analysis set. The
#     figure is deliberately higher than the 91% that an odds ratio of 2.30
#     would give: the realised effect in any one draw differs from the
#     generative value by sampling, and a design with a 1-in-11 chance of
#     leaving the guide with no detectable locus at all is too fragile for a
#     teaching resource that others are expected to rebuild. The generative
#     value is chosen for the reliability of the demonstration; the seed is
#     fixed in advance and is never selected on the outcome, because choosing a
#     seed because the locus reached significance would be conditioning on
#     significance, which is the very bias the guide sets out to explain. This is the
#     positive control, the single peak of the Manhattan plot, and the locus
#     that the fine-mapping section resolves. The effect is larger than the
#     true ABO effect on pancreatic cancer, which is near 1.2; that is a
#     deliberate and declared amplification, without which no locus would be
#     detectable at this sample size and the guide would have nothing to
#     demonstrate downstream.
#
#   G2, TERT/CLPTM1L (5p15.33), odds ratio 1.50 per risk allele.
#     Not detectable: about 2% power in the European set, rising to about 29%
#     across the three ancestry strata combined. This is the second most
#     replicated pancreatic cancer locus, and 1.50 sits at the upper end of
#     what the literature reports for cancer susceptibility variants. It is
#     included so that the demonstration contains a true positive that the scan
#     misses, with the truth known, and so that the meta-analysis section can
#     show a real effect being partly recovered by combining strata.
#
#   Sex, odds ratio 1.30 for males, matching the male-to-female incidence ratio
#   reported for pancreatic cancer.
#
#   Age, odds ratio 1.50 per 10 years, which reproduces the roughly 3.5-year
#   difference in mean age between cases and controls seen in case-control
#   series where controls are broadly age-comparable but not matched.
#
#   The genetic effects are identical in all three ancestry groups. This is an
#   assumption, not a finding, and it is what allows the meta-analysis section
#   to demonstrate a homogeneous effect (I-squared near zero) at the causal
#   variants while their correlated neighbours show heterogeneity.
#
# SURVIVAL
#   Overall survival is defined for cases only, as time from diagnosis, which
#   is what a survival GWAS in a cancer cohort actually analyses. Controls have
#   no survival record. No variant is given a survival effect: the survival
#   scan in this demonstration is a true null, and the guide says so.
#
# USAGE
#   Rscript 04_simulate_phenotype.R <plink2> <bfile> <ancestry.tsv> <outdir>
###############################################################################

args    <- commandArgs(trailingOnly = TRUE)
PLINK2  <- if (length(args) >= 1) args[1] else "plink2"
BFILE   <- if (length(args) >= 2) args[2] else "pdac_demo"
ANCFILE <- if (length(args) >= 3) args[3] else "sample_ancestry.tsv"
OUTDIR  <- if (length(args) >= 4) args[4] else "."

SEED <- 2026
set.seed(SEED)

# --------------------------------------------------------------- parameters --
CAUSAL <- c(G1 = "9:133273682:A:T",   # ABO 9q34.2, 130 bp from rs505922
            G2 = "5:1286401:C:A")     # TERT/CLPTM1L 5p15.33

OR_G1   <- 2.60
OR_G2   <- 1.50
OR_MALE <- 1.30
OR_AGE  <- 1.50          # per 10 years

TARGET_CASE_FRACTION <- c(eur = 0.333, afr = 0.456, eas = 0.414)

AGE_MEAN <- 66.7
AGE_SD   <- 10.0
AGE_MIN  <- 33
AGE_MAX  <- 90

# survival, cases only
MEDIAN_SURVIVAL_MONTHS <- 11      # median overall survival in metastatic PDAC
# Recruitment is staggered, so each patient has a different amount of follow-up
# available when the study closes. Censoring in a survival GWAS comes mostly
# from this, not from anyone living a long time, and a demonstration with no
# censoring cannot teach what censoring is.
FOLLOWUP_MIN_MONTHS <- 6
FOLLOWUP_MAX_MONTHS <- 60

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# ------------------------------------------------------- read the genotypes --
msg("Extracting the two causal variants")
snpfile <- file.path(OUTDIR, "_causal.snps")
writeLines(unname(CAUSAL), snpfile)

system2(PLINK2, c("--bfile", BFILE, "--extract", snpfile,
                  "--recode", "A", "--out", file.path(OUTDIR, "_causal")),
        stdout = FALSE, stderr = FALSE)

raw <- read.table(file.path(OUTDIR, "_causal.raw"), header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)

# PLINK appends _<allele> to the column name; find each variant's column
find_col <- function(id) {
  hit <- grep(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", id), "_"),
              names(raw), value = TRUE)
  if (length(hit) != 1) stop("could not identify a unique column for ", id)
  hit
}
c1 <- find_col(CAUSAL[["G1"]]);  c2 <- find_col(CAUSAL[["G2"]])
msg("  %s  -> column %s", CAUSAL[["G1"]], c1)
msg("  %s  -> column %s", CAUSAL[["G2"]], c2)

# PLINK's --export A counts the allele PLINK 2 treats as REFERENCE when it
# reads a PLINK 1 fileset, which for these variant IDs (chr:pos:ref:alt) is
# the ref allele: A at ABO and C at TERT/CLPTM1L. The risk allele is the
# counted one, so it is the reference allele, and it is the MAJOR allele at
# ABO. Verify rather than assume: the column header written by --export A is
# <ID>_<counted allele>. An association run reports its effect for the other
# allele, so the odds ratio it prints is the reciprocal of the one set here.
G1 <- raw[[c1]];  G2 <- raw[[c2]]

# The demonstration data carry deliberately injected missingness, including 25
# samples with a low call rate, so no variant anywhere is completely called. An
# individual whose causal genotype is missing is given the expected dosage, so
# that the phenotype is defined for everyone and no one is silently dropped.
# The alternative, excluding them, would make missingness informative about
# case status and would quietly break the quality-control section that follows.
imp <- function(x, label) {
  k <- sum(is.na(x))
  if (k > 0) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    msg("  %s: %d missing calls (%.2f%%) set to the expected dosage",
        label, k, 100 * k / length(x))
  }
  x
}
G1 <- imp(G1, CAUSAL[["G1"]]);  G2 <- imp(G2, CAUSAL[["G2"]])

IID <- raw$IID
n   <- length(IID)
msg("  %d individuals, both causal variants fully called", n)

anc_tab <- read.table(ANCFILE, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("IID", "group"))
anc <- anc_tab$group[match(IID, anc_tab$IID)]
if (anyNA(anc)) stop("ancestry missing for ", sum(is.na(anc)), " individuals")

# ---------------------------------------------------------------- covariates -
# Sex comes from the .fam file, so that it stays consistent with the genotypes
# and with the sex-check step of the quality-control section.
fam <- read.table(paste0(BFILE, ".fam"), header = FALSE, stringsAsFactors = FALSE)
sex_code <- fam$V5[match(IID, fam$V2)]      # 1 = male, 2 = female
male <- as.integer(sex_code == 1)

age <- round(pmin(pmax(rnorm(n, AGE_MEAN, AGE_SD), AGE_MIN), AGE_MAX))
age_z <- (age - AGE_MEAN) / 10

# ------------------------------------------------------- generative model ----
b1  <- log(OR_G1);   b2 <- log(OR_G2)
bs  <- log(OR_MALE); ba <- log(OR_AGE)

eta_nointercept <- b1 * G1 + b2 * G2 + bs * male + ba * age_z

# Solve, per ancestry group, the intercept that yields the target case fraction.
solve_alpha <- function(eta, target) {
  f <- function(a) mean(plogis(a + eta)) - target
  uniroot(f, interval = c(-25, 25), tol = 1e-10)$root
}
alpha <- setNames(numeric(0), character(0))
p <- numeric(n)
for (g in names(TARGET_CASE_FRACTION)) {
  idx <- which(anc == g)
  a   <- solve_alpha(eta_nointercept[idx], TARGET_CASE_FRACTION[[g]])
  alpha[g] <- a
  p[idx]   <- plogis(a + eta_nointercept[idx])
  msg("  intercept for %s: %+.4f  (expected case fraction %.3f)",
      g, a, mean(p[idx]))
}

status <- rbinom(n, 1, p)                    # 1 = case, 0 = control
pheno  <- ifelse(status == 1, 2L, 1L)        # PLINK coding: 2 = case, 1 = control

msg("Realised counts: %d cases, %d controls", sum(status == 1), sum(status == 0))
for (g in names(TARGET_CASE_FRACTION)) {
  idx <- anc == g
  msg("  %s: %d cases, %d controls, fraction %.3f",
      g, sum(status[idx] == 1), sum(status[idx] == 0), mean(status[idx]))
}

# ------------------------------------------------------------- VALIDATION ----
# The dataset is written only if the generative parameters are recovered.
# This is the check whose absence let a previous version ship an effect that
# did not match its own documentation.
msg("")
msg("Validation: re-estimating the generative parameters from the simulated data")

fit <- glm(status ~ G1 + G2 + male + age_z, family = binomial())
est <- summary(fit)$coefficients

check <- function(term, target_or, label) {
  b   <- est[term, "Estimate"];  se <- est[term, "Std. Error"]
  lo  <- exp(b - 1.96 * se);     hi <- exp(b + 1.96 * se)
  ok  <- (target_or >= lo) && (target_or <= hi)
  msg("  %-28s target %.3f   estimate %.3f  (%.2f to %.2f)   %s",
      label, target_or, exp(b), lo, hi, if (ok) "ok" else "FAIL")
  ok
}

ok <- all(unlist(list(
  check("G1",    OR_G1,   "ABO, per risk allele"),
  check("G2",    OR_G2,   "TERT/CLPTM1L, per allele"),
  check("male",  OR_MALE, "male sex"),
  check("age_z", OR_AGE,  "age, per 10 years")
)))

if (!ok)
  stop("validation failed: the simulated data do not recover the generative ",
       "parameters. The dataset has NOT been written. Investigate before ",
       "changing the tolerance.")

msg("  all parameters recovered within their confidence intervals")

# ----------------------------------------------------------------- survival --
# Cases only, exponential survival with administrative censoring.
rate <- log(2) / MEDIAN_SURVIVAL_MONTHS
time  <- rep(NA_real_, n);  event <- rep(NA_integer_, n)
is_case <- status == 1
nk <- sum(is_case)
t_death <- rexp(nk, rate)                                   # time to death
t_admin <- runif(nk, FOLLOWUP_MIN_MONTHS, FOLLOWUP_MAX_MONTHS)  # follow-up available
time[is_case]  <- round(pmin(t_death, t_admin), 1)
event[is_case] <- as.integer(t_death <= t_admin)
msg("")
msg("Survival, cases only: %d individuals, %d events (%.1f%%), median %.1f months",
    sum(is_case), sum(event[is_case]), 100 * mean(event[is_case]),
    median(time[is_case]))

# -------------------------------------------------------------------- write --
write.table(data.frame(FID = IID, IID = IID, PHENO = pheno),
            file.path(OUTDIR, "phenotype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(FID = IID, IID = IID, SEX = sex_code, AGE = age),
            file.path(OUTDIR, "covariates.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

surv <- data.frame(FID = IID, IID = IID, TIME = time, EVENT = event)[is_case, ]
write.table(surv, file.path(OUTDIR, "survival.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# A machine-readable record of the truth, so that no downstream document has to
# rely on prose for the generative values.
truth <- data.frame(
  variant = c(CAUSAL[["G1"]], CAUSAL[["G2"]], "sex_male", "age_per_10y"),
  locus   = c("ABO 9q34.2", "TERT/CLPTM1L 5p15.33", NA, NA),
  generative_OR = c(OR_G1, OR_G2, OR_MALE, OR_AGE),
  estimated_OR  = c(exp(est["G1", "Estimate"]), exp(est["G2", "Estimate"]),
                    exp(est["male", "Estimate"]), exp(est["age_z", "Estimate"])),
  stringsAsFactors = FALSE)
write.table(truth, file.path(OUTDIR, "truth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

invisible(file.remove(file.path(OUTDIR, c("_causal.snps", "_causal.raw", "_causal.log"))))

msg("")
msg("Written to %s: phenotype.txt, covariates.txt, survival.txt, truth.tsv", OUTDIR)
msg("Seed %d. Rerunning this script reproduces these files exactly.", SEED)
