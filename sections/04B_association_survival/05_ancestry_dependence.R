#!/usr/bin/env Rscript

################################################################################
# Section 4B: does the simulated survival depend on ancestry?
#
# PURPOSE
#   Test directly whether the survival outcome carries structure that the
#   principal components could remove. If it does not, inflation of the Cox
#   scan cannot be confounding by ancestry, and adding components can only
#   cost calibration (see 03_covariate_cost.R and 04_permutation_null.R).
#
# WHY
#   The simulation gives survival no dependence on genotype, ancestry, SEX or
#   AGE (demo_dataset/scripts/04_simulate_phenotype.R). This step checks that
#   claim on the released data rather than assuming it, the same rule the
#   generator applies to the case-control phenotype: a simulated truth is
#   used only once the data are shown to reproduce it.
#
# WHAT IT DOES
#   1. All cases with survival data (demo_data, before QC): events and median
#      time by ancestry group; Cox model on ancestry group, likelihood-ratio
#      test on 2 df, unadjusted and given SEX and AGE; log-rank test.
#   2. The European survival set of 01_cox_survival.R (individuals in the
#      within-set components): likelihood-ratio test of within-set PC1-PC10
#      jointly, given SEX and AGE (10 df); each component alone, given SEX and
#      AGE (Wald); the Section 4A components PC1-PC10 jointly, as a check that
#      the answer does not depend on which components are used.
#   3. Censoring pattern of the survival set: distinct times, tied event
#      times, zero times.
#
# INPUT
#   demo_data/survival.txt, sample_ancestry.tsv, covariates.txt, phenotype.txt
#   results/assoc/pdac_demo_04B_pca_within.eigenvec   (01_cox_survival.R)
#   results/assoc/pdac_demo_04A_covar.txt             (Section 4A)
#
# OUTPUT (results/assoc/, moved to sections/04B_association_survival/results/)
#   pdac_demo_04B_ancestry_dependence.tsv          one row per test
#   pdac_demo_04B_ancestry_dependence_groups.tsv   events and times by group
#   pdac_demo_04B_ancestry_dependence.txt          plain-text report
#
# No random numbers are drawn; the output is deterministic.
#
# Derived from the diagnostic of 22 November 2026 (ancestry_check.R, outside
# the repository). Added here: the likelihood-ratio test of ancestry group
# given SEX and AGE, in place of per-coefficient Wald P values only.
################################################################################

suppressPackageStartupMessages(library(survival))

## ---- project root, as in the other 4B scripts --------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
root <- script_dir
for (up in list(c("..", ".."), c("..", "..", ".."))) {
  cand <- normalizePath(do.call(file.path, c(list(script_dir), as.list(up))), mustWork = FALSE)
  if (dir.exists(file.path(cand, "scripts", "dev"))) { root <- cand; break }
}
if (!dir.exists(file.path(root, "scripts", "dev"))) root <- getwd()
setwd(root)

EIGENVEC <- "results/assoc/pdac_demo_04B_pca_within.eigenvec"
COVAR4A  <- "results/assoc/pdac_demo_04A_covar.txt"
OUT_DIR  <- "results/assoc"
N_PC     <- 10L
for (f in c(EIGENVEC, COVAR4A, "demo_data/survival.txt", "demo_data/sample_ancestry.tsv",
            "demo_data/covariates.txt", "demo_data/phenotype.txt"))
  if (!file.exists(f)) stop("missing input: ", f, " (run Section 4A and 01_cox_survival.R first)")

report <- character()
w <- function(...) { line <- paste0(...); report <<- c(report, line); cat(line, "\n", sep = "") }
tests <- list()
add <- function(set, d, test, df, stat, p)
  tests[[length(tests) + 1L]] <<- data.frame(set = set, n = nrow(d), events = sum(d$EVENT),
    test = test, df = df, statistic = round(stat, 4), P = signif(p, 4))
lrt <- function(f1, f0) 2 * (tail(f1$loglik, 1) - tail(f0$loglik, 1))   # a null model has one loglik

## ---- data -------------------------------------------------------------------------
surv <- read.delim("demo_data/survival.txt", stringsAsFactors = FALSE)
anc  <- read.delim("demo_data/sample_ancestry.tsv", header = FALSE, col.names = c("IID", "anc"),
                   stringsAsFactors = FALSE)
cov  <- read.delim("demo_data/covariates.txt", stringsAsFactors = FALSE)
ph   <- read.delim("demo_data/phenotype.txt", stringsAsFactors = FALSE)
d <- merge(merge(merge(surv[, c("IID", "TIME", "EVENT")], anc, by = "IID"),
                 cov[, c("IID", "SEX", "AGE")], by = "IID"), ph[, c("IID", "PHENO")], by = "IID")
d$anc <- factor(d$anc)

w("SECTION 4B: DOES THE SIMULATED SURVIVAL DEPEND ON ANCESTRY?")
w("")
w("1. ALL CASES WITH SURVIVAL DATA (demo_data, before QC)")
w(sprintf("   n %d, events %d; case status of these individuals: %s (2 = case)",
          nrow(d), sum(d$EVENT), paste(names(table(d$PHENO)), table(d$PHENO), collapse = ", ")))
groups <- do.call(rbind, lapply(levels(d$anc), function(g) {
  x <- d[d$anc == g, ]
  km <- survfit(Surv(TIME, EVENT) ~ 1, x)
  data.frame(group = g, n = nrow(x), events = sum(x$EVENT),
             pct_events = round(100 * mean(x$EVENT), 1),
             median_time = round(unname(summary(km)$table["median"]), 1),
             mean_time = round(mean(x$TIME), 1))
}))
for (i in seq_len(nrow(groups))) with(groups[i, ], w(sprintf(
  "   %-4s n %3d  events %3d (%.1f%%)  median time %.1f  mean time %.1f",
  group, n, events, pct_events, median_time, mean_time)))

f0  <- coxph(Surv(TIME, EVENT) ~ 1, d)
fa  <- coxph(Surv(TIME, EVENT) ~ anc, d)
fsa <- coxph(Surv(TIME, EVENT) ~ SEX + AGE, d)
fas <- coxph(Surv(TIME, EVENT) ~ anc + SEX + AGE, d)
s <- lrt(fa, f0);   add("all_cases", d, "ancestry group, Cox LR", 2L, s, pchisq(s, 2, lower.tail = FALSE))
w(sprintf("   Cox on ancestry group: LR chi2 %.3f, 2 df, P %.3g", s, pchisq(s, 2, lower.tail = FALSE)))
hr <- summary(fa)$conf.int
for (r in rownames(hr)) w(sprintf("     HR %-8s %.3f (95%% CI %.3f to %.3f), reference %s",
                                  sub("^anc", "", r), hr[r, 1], hr[r, 3], hr[r, 4], levels(d$anc)[1]))
s <- lrt(fas, fsa); add("all_cases", d, "ancestry group given SEX + AGE, Cox LR", 2L, s, pchisq(s, 2, lower.tail = FALSE))
w(sprintf("   ancestry group given SEX and AGE: LR chi2 %.3f, 2 df, P %.3g", s, pchisq(s, 2, lower.tail = FALSE)))
cf <- summary(fas)$coefficients
for (r in c("SEX", "AGE")) {
  add("all_cases", d, paste(r, "given ancestry group and the other, Cox Wald"), 1L, cf[r, "z"]^2, cf[r, "Pr(>|z|)"])
  w(sprintf("     %s: HR %.4f, P %.3g", r, cf[r, "exp(coef)"], cf[r, "Pr(>|z|)"]))
}
lt <- survdiff(Surv(TIME, EVENT) ~ anc, d)
add("all_cases", d, "ancestry group, log-rank", 2L, lt$chisq, pchisq(lt$chisq, 2, lower.tail = FALSE))
w(sprintf("   log-rank across groups: chi2 %.3f, 2 df, P %.3g", lt$chisq, pchisq(lt$chisq, 2, lower.tail = FALSE)))

## ---- the European survival set of 01_cox_survival.R -----------------------------------
w("")
w("2. THE EUROPEAN SURVIVAL SET OF 01_cox_survival.R")
pcs <- paste0("PC", seq_len(N_PC))
ev <- read.table(EIGENVEC, header = TRUE, comment.char = "", check.names = FALSE)
names(ev)[1:2] <- c("FID", "IID")
e <- merge(d, ev[, c("IID", pcs)], by = "IID")
w(sprintf("   n %d, events %d", nrow(e), sum(e$EVENT)))
b   <- coxph(Surv(TIME, EVENT) ~ SEX + AGE, e)
fpc <- coxph(as.formula(paste("Surv(TIME, EVENT) ~ SEX + AGE +", paste(pcs, collapse = " + "))), e)
s <- lrt(fpc, b)
add("eur_survival_set", e, "within-set PC1-PC10 jointly given SEX + AGE, Cox LR", N_PC, s,
    pchisq(s, N_PC, lower.tail = FALSE))
w(sprintf("   within-set PC1-PC10 jointly, given SEX and AGE: LR chi2 %.3f, %d df, P %.3g",
          s, N_PC, pchisq(s, N_PC, lower.tail = FALSE)))
w("   each component alone, given SEX and AGE (Wald):")
for (p in pcs) {
  cp <- summary(coxph(as.formula(paste("Surv(TIME, EVENT) ~ SEX + AGE +", p)), e))$coefficients
  add("eur_survival_set", e, paste("within-set", p, "given SEX + AGE, Cox Wald"), 1L, cp[p, "z"]^2, cp[p, "Pr(>|z|)"])
  w(sprintf("     %-4s P %.3g", p, cp[p, "Pr(>|z|)"]))
}
cb <- summary(b)$coefficients
for (r in c("SEX", "AGE")) {
  add("eur_survival_set", e, paste(r, "given the other, Cox Wald"), 1L, cb[r, "z"]^2, cb[r, "Pr(>|z|)"])
  w(sprintf("   %s: HR %.4f, P %.3g", r, cb[r, "exp(coef)"], cb[r, "Pr(>|z|)"]))
}
c4 <- read.delim(COVAR4A, check.names = FALSE)
names(c4)[1] <- sub("^#", "", names(c4)[1])
e4 <- merge(d, c4[, c("IID", pcs)], by = "IID")
b4 <- coxph(Surv(TIME, EVENT) ~ SEX + AGE, e4)
f4 <- coxph(as.formula(paste("Surv(TIME, EVENT) ~ SEX + AGE +", paste(pcs, collapse = " + "))), e4)
s <- lrt(f4, b4)
add("eur_survival_set_4A_components", e4, "Section 4A PC1-PC10 jointly given SEX + AGE, Cox LR", N_PC, s,
    pchisq(s, N_PC, lower.tail = FALSE))
w(sprintf("   Section 4A PC1-PC10 jointly, given SEX and AGE (n %d): LR chi2 %.3f, %d df, P %.3g",
          nrow(e4), s, N_PC, pchisq(s, N_PC, lower.tail = FALSE)))

## ---- censoring pattern --------------------------------------------------------------
w("")
w("3. CENSORING PATTERN OF THE EUROPEAN SURVIVAL SET")
te <- e$TIME[e$EVENT == 1]
w(sprintf("   distinct times among %d individuals: %d", nrow(e), length(unique(e$TIME))))
w(sprintf("   events sharing their time with another event: %d", sum(duplicated(te) | duplicated(te, fromLast = TRUE))))
w(sprintf("   times equal to zero: %d", sum(e$TIME == 0)))

w("")
w("HOW TO READ THIS")
w("  Large P values in sections 1 and 2 mean the outcome carries no detectable")
w("  ancestry structure, at the group level or along the components. Inflation")
w("  of the Cox scan then cannot be confounding by ancestry, and components")
w("  added to the model remove nothing; compare 04_permutation_null.R.")

write.table(do.call(rbind, tests), file.path(OUT_DIR, "pdac_demo_04B_ancestry_dependence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(groups, file.path(OUT_DIR, "pdac_demo_04B_ancestry_dependence_groups.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(report, file.path(OUT_DIR, "pdac_demo_04B_ancestry_dependence.txt"))
