#!/usr/bin/env Rscript

################################################################################
# Section 8, step 06: the number of principal components over the replicates
#
# PURPOSE
#   The released phenotype gives one value of lambda per number of components
#   (04A step 01b), with no measure of its variation, and the replicates of
#   step 02 use only the number chosen there. This step repeats the European
#   scan of every replicate with 0, 2, 5, 10 and 20 within-European components
#   and applies the selection rule of 04A step 01b to every replicate.
#
# WHAT IT DOES, per replicate phenotype of step 01
#   1. Genome-wide European scan (autosomes, logistic regression with
#      firth-fallback, SEX, AGE and PC1..PCk) for k = 0, 2, 5, 10, 20: the
#      command of 04A step 01b, with the replicate's phenotype and age.
#   2. lambda on the LD-pruned European set (the set step 02 uses), on all
#      variants (the set the rule of 01b uses), and on the 5% of the pruned
#      variants most associated with PC1 and PC2 (see below); P, estimate and
#      standard error at the three causal variants; lead variant and the ten
#      smallest P.
#   3. The rule of 04A step 01b, with its constants and its candidates
#      (2, 5, 10, 20): scree floor e (the same for every replicate: it depends
#      on the genotypes only), association floor a (each of the 20 components
#      tested against the replicate's case status given SEX and AGE, at
#      P < 0.05 / 20), one-sided stability of lambda (tolerance 0.01), of the
#      lead variant and of the ten smallest P. Zero is not a candidate of the
#      rule; whether the scan without components would pass the stability
#      criterion is reported separately.
#   CHECK: at the number of components of the canonical run the scan must
#   reproduce step 02 (lambda on the pruned set, P at the causal variants);
#   the step stops otherwise. The agreement is not exact to the last digit:
#   step 02 scans an extracted subset, this step all variants, and PLINK 2
#   returns a P that differs from the fourth significant digit at about one
#   variant in twenty (4,186 of 84,816 in replicate 1; those inspected had
#   missing genotypes) when the variants loaded beside it change. We
#   attribute it to the iterative fit being started from the previous
#   variant in single precision; this was not verified in the PLINK source.
#   Tolerances: 0.001 on lambda, 0.01 on log10 P.
#
# VARIANTS THAT FOLLOW PC1 AND PC2
#   Each pruned variant is regressed on PC1 and on PC2 in the European set
#   (plink2 --glm, the component as a quantitative phenotype, no covariates);
#   R2 = t2 / (t2 + n - 2), summed over the two components, which are
#   orthogonal. The 5% with the largest sum are the variants a scan without
#   components would get wrong first if case status followed those axes.
#
# NOT PART OF THE CANONICAL RUN: reads its outputs in results/, as step 02.
#
# SETTINGS
#   GWAS_08_PARALLEL (default 4) replicates at a time, GWAS_08_THREADS (4) each.
#   About 1,000 genome-wide scans: two hours on 16 cores.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_npc_per_replicate.tsv   one row per replicate and number of components
#   pdac_demo_08_npc_rule.tsv            one row per replicate: e, a, stability, choice
#   pdac_demo_08_npc_summary.tsv         per number of components: lambda and power
#   pdac_demo_08_npc_rule_summary.tsv    distribution of the number chosen
#   pdac_demo_08_npc_summary.txt         plain-text report
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
EIGENVAL <- "results/pca/pdac_demo_02_pca_eur.eigenval"
PRUNE <- "results/pca/pdac_demo_02_prune_eur.prune.in"
CAUSAL <- c(G1 = "9:133273682:A:T", G2 = "5:1286401:C:A", G3 = "13:73362027:C:T")
LABEL <- c(G1 = "ABO", G2 = "TERT/CLPTM1L", G3 = "13q22.1")
KS <- c(0L, 2L, 5L, 10L, 20L)
## constants of 04A step 01b
CANDIDATES <- c(2L, 5L, 10L, 20L); N_TESTED <- 20L; ALPHA_ASSOC <- 0.05 / N_TESTED
LAMBDA_TOL <- 0.01; TOP_N <- 10L; TOP_MIN_SHARED <- 8L
ALPHA <- 5e-8; TOP_FRAC <- 0.05
N_PAR <- as.integer(Sys.getenv("GWAS_08_PARALLEL", "4"))
THREADS <- as.integer(Sys.getenv("GWAS_08_THREADS", "4"))
for (f in c(paste0(IN, ".bed"), EUR_KEEP, EUR_COV, EIGENVAL, PRUNE, file.path(OUT, "pdac_demo_08_per_replicate.tsv"),
            file.path(OUT, "pdac_demo_08_lambda.tsv")))
  if (!file.exists(f)) stop("missing input: ", f)
reps <- sort(Sys.glob(file.path(SIM, "rep_*.pheno")))
if (!length(reps)) stop("no replicates in ", SIM, ": run step 01")
rd <- function(f, ...) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, ...)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
lam <- function(p) { p <- p[is.finite(p) & p > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
tmp <- file.path(tempdir(), "npc"); dir.create(tmp, showWarnings = FALSE)

cov_e <- rd(EUR_COV)
if (!all(paste0("PC", 1:N_TESTED) %in% names(cov_e))) stop(EUR_COV, " lacks PC1-PC", N_TESTED)
keep_ids <- read.table(EUR_KEEP, stringsAsFactors = FALSE)[[2]]
pruned <- readLines(PRUNE)
ev <- scan(EIGENVAL, quiet = TRUE); e_floor <- which.max(ev[1:10] / ev[2:11])

## counted allele of the generator, to orient the estimates (as step 02)
writeLines(unname(CAUSAL), file.path(tmp, "c.snps"))
system2("plink2", c("--bfile", "demo_data/pdac_demo", "--extract", file.path(tmp, "c.snps"), "--recode", "A",
                    "--out", file.path(tmp, "c")), stdout = FALSE, stderr = FALSE)
hdr <- names(read.table(file.path(tmp, "c.raw"), header = TRUE, check.names = FALSE, nrows = 1))
counted <- sapply(CAUSAL, function(id) sub(".*_", "", grep(paste0("^", gsub("([.:])", "\\\\\\1", id), "_"), hdr, value = TRUE)))

## ---- pruned variants that follow PC1 and PC2 -------------------------------------------
pcf <- file.path(tmp, "pc.pheno")
x <- cov_e[cov_e$IID %in% keep_ids, c("IID", "PC1", "PC2")]
write.table(data.frame(`#FID` = x$IID, x, check.names = FALSE), pcf, sep = "\t", quote = FALSE, row.names = FALSE)
if (system2("plink2", c("--bfile", IN, "--keep", EUR_KEEP, "--extract", PRUNE, "--chr", "1-22", "--pheno", pcf,
                        "--glm", "allow-no-covars", "--threads", N_PAR * THREADS, "--out", file.path(tmp, "pc")),
            stdout = FALSE, stderr = FALSE) != 0) stop("plink2 failed: components as phenotypes")
r2 <- lapply(c("PC1", "PC2"), function(pc) { d <- rd(file.path(tmp, paste0("pc.", pc, ".glm.linear")))
  if ("TEST" %in% names(d)) d <- d[d$TEST == "ADD", ]
  setNames(d$T_STAT^2 / (d$T_STAT^2 + d$OBS_CT - 2), d$ID) })
ids12 <- intersect(names(r2[[1]]), names(r2[[2]]))
r2sum <- r2[[1]][ids12] + r2[[2]][ids12]; r2sum <- r2sum[is.finite(r2sum)]
pc12 <- names(r2sum)[r2sum >= quantile(r2sum, 1 - TOP_FRAC)]
cat(sprintf("pruned variants %d; following PC1 and PC2 (top %g%% of R2): %d, R2 at least %.4f, median R2 of all %.5f\n",
            length(r2sum), 100 * TOP_FRAC, length(pc12), min(r2sum[pc12]), median(r2sum)))

## ---- one replicate ------------------------------------------------------------------------
one <- function(f) {
  i <- as.integer(sub(".*rep_([0-9]+)\\.pheno$", "\\1", f)); tag <- sprintf("r%03d", i)
  age <- read.delim(sub("\\.pheno$", ".age", f)); ph <- read.delim(f)
  cv <- cov_e; cv$AGE <- age$AGE[match(cv$IID, age$IID)]
  cf <- file.path(tmp, paste0(tag, ".cov")); cvw <- cv; names(cvw)[1] <- paste0("#", names(cvw)[1])
  write.table(cvw, cf, sep = "\t", quote = FALSE, row.names = FALSE)
  ## association floor a, as 01b (b)
  d <- merge(cv, ph[, c("IID", "PHENO")], by = "IID"); d$status <- as.integer(d$PHENO == 2)
  pa <- sapply(paste0("PC", 1:N_TESTED), function(pc)
    summary(glm(as.formula(paste("status ~ SEX + AGE +", pc)), data = d, family = binomial()))$coefficients[pc, 4])
  a_idx <- which(pa < ALPHA_ASSOC); a <- if (length(a_idx)) max(a_idx) else 0L
  sc <- lapply(KS, function(k) {
    out <- file.path(tmp, paste0(tag, "_k", k))
    st <- system2("plink2", c("--bfile", IN, "--chr", "1-22", "--keep", EUR_KEEP, "--pheno", f, "--covar", cf,
                              "--covar-name", paste(c("SEX", "AGE", if (k > 0) paste0("PC", seq_len(k))), collapse = ","),
                              "--covar-variance-standardize", "--threads", THREADS,
                              "--glm", "firth-fallback", "hide-covar", "cols=+a1freq,+beta,+orbeta,+nobs,+err",
                              "--out", out), stdout = FALSE, stderr = FALSE)
    if (st != 0) stop("plink2 failed for ", tag, " k ", k)
    gf <- Sys.glob(paste0(out, ".*.glm.logistic.hybrid"))[1]
    h <- sub("^#", "", strsplit(readLines(gf, n = 1), "\t")[[1]])
    want <- c("ID", "A1", "TEST", "BETA", "SE", "P")
    g <- rd(gf, sep = "\t", colClasses = ifelse(h %in% want, NA, "NULL"))
    unlink(Sys.glob(paste0(out, "*")))
    if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
    ca <- g[match(CAUSAL, g$ID), ]
    g <- g[is.finite(g$P) & g$P > 0, ]; o <- order(g$P)
    sg <- ifelse(ca$A1 == counted, 1, -1)
    list(k = k, top = g$ID[o][seq_len(TOP_N)],
         row = data.frame(replicate = i, n_pcs = k, variants = nrow(g), lambda_pruned = lam(g$P[g$ID %in% pruned]),
                          lambda_all = lam(g$P), lambda_pc12 = lam(g$P[g$ID %in% pc12]),
                          n_p_lt_5e8 = sum(g$P < ALPHA), lead_id = g$ID[o][1], lead_p = g$P[o][1],
                          setNames(as.list(ca$P), paste0("p_", names(CAUSAL))),
                          setNames(as.list(sg * ca$BETA), paste0("log_or_", names(CAUSAL))),
                          setNames(as.list(ca$SE), paste0("se_", names(CAUSAL))), stringsAsFactors = FALSE))
  })
  names(sc) <- KS
  tab <- do.call(rbind, lapply(sc, `[[`, "row"))
  ## stability, as 01b (c); lambda rounded to four decimals as in its table
  l4 <- setNames(round(tab$lambda_all, 4), KS)
  stable <- sapply(KS, function(k) { larger <- KS[KS > k]
    if (!length(larger)) return(TRUE)
    all(sapply(larger, function(k2) l4[[as.character(k2)]] >= l4[[as.character(k)]] - LAMBDA_TOL &&
          tab$lead_id[tab$n_pcs == k] == tab$lead_id[tab$n_pcs == k2] &&
          length(intersect(sc[[as.character(k)]]$top, sc[[as.character(k2)]]$top)) >= TOP_MIN_SHARED)) })
  names(stable) <- KS
  ## why a candidate is not stable: lambda, lead or top ten, against any larger candidate
  why <- function(k) { larger <- KS[KS > k]; if (!length(larger)) return(c(FALSE, FALSE, FALSE))
    c(any(sapply(larger, function(k2) l4[[as.character(k2)]] < l4[[as.character(k)]] - LAMBDA_TOL)),
      any(sapply(larger, function(k2) tab$lead_id[tab$n_pcs == k] != tab$lead_id[tab$n_pcs == k2])),
      any(sapply(larger, function(k2) length(intersect(sc[[as.character(k)]]$top, sc[[as.character(k2)]]$top)) < TOP_MIN_SHARED))) }
  w2 <- why(2L)
  floor_k <- max(e_floor, a)
  ok <- CANDIDATES[CANDIDATES >= floor_k & stable[as.character(CANDIDATES)]]
  chosen <- if (length(ok)) min(ok) else max(CANDIDATES)
  tab$stable <- stable[as.character(tab$n_pcs)]
  list(tab = tab, rule = data.frame(replicate = i, scree_floor_e = e_floor, association_floor_a = a,
    min_assoc_p = min(pa), n_components_nominal_p_lt_0.05 = sum(pa < 0.05),
    stable_0 = stable[["0"]], stable_2 = stable[["2"]], stable_5 = stable[["5"]], stable_10 = stable[["10"]],
    k2_fails_lambda = w2[1], k2_fails_lead = w2[2], k2_fails_top10 = w2[3], chosen_n_pcs = chosen))
}

t0 <- Sys.time()
cat(sprintf("%d replicates, %d in parallel, %d threads each; components %s\n", length(reps), N_PAR, THREADS, paste(KS, collapse = ", ")))
res <- mclapply(reps, function(f) tryCatch(one(f), error = function(e) e), mc.cores = N_PAR, mc.preschedule = FALSE)
bad <- which(sapply(res, function(x) inherits(x, "error") || inherits(x, "try-error")))
if (length(bad)) stop("replicates failed: ", paste(basename(reps[bad]), collapse = ", "), ": ", as.character(res[[bad[1]]]))
tab <- do.call(rbind, lapply(res, `[[`, "tab")); rule <- do.call(rbind, lapply(res, `[[`, "rule"))
tab <- tab[order(tab$replicate, tab$n_pcs), ]; rule <- rule[order(rule$replicate), ]

## ---- check against step 02 at k = 2 ---------------------------------------------------------
k_run <- as.integer(read.table("results/assoc/pdac_demo_04A_npc_chosen.txt")[1, 2])
if (k_run %in% KS) {
  per <- read.delim(file.path(OUT, "pdac_demo_08_per_replicate.tsv"), stringsAsFactors = FALSE)
  l2 <- read.delim(file.path(OUT, "pdac_demo_08_lambda.tsv"))
  t2 <- tab[tab$n_pcs == k_run, ]
  dl <- max(abs(t2$lambda_pruned - l2$lambda_eur_4a[match(t2$replicate, l2$replicate)]))
  dp <- max(sapply(names(CAUSAL), function(L) { x <- per[per$locus == L & per$analysis == "eur_4a", ]
    max(abs(log10(t2[[paste0("p_", L)]]) - log10(x$p[match(t2$replicate, x$replicate)]))) }))
  chk_line <- sprintf("check against step 02 at %d components: max |lambda difference| %.2g, max |log10 P difference| %.2g", k_run, dl, dp); cat(chk_line, "\n")
  if (dl > 1e-3 || dp > 1e-2) stop("the scan at ", k_run, " components does not reproduce step 02")
}

## ---- summaries ----------------------------------------------------------------------------
wilson <- function(x, n) { p <- x / n; z <- 1.959964
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
summ <- do.call(rbind, lapply(KS, function(k) { x <- tab[tab$n_pcs == k, ]; n <- nrow(x)
  o <- data.frame(n_pcs = k, replicates = n)
  for (cn in c("lambda_pruned", "lambda_all", "lambda_pc12")) {
    o[[paste0(cn, "_mean")]] <- mean(x[[cn]]); o[[paste0(cn, "_sd")]] <- sd(x[[cn]])
    o[[paste0(cn, "_se_of_mean")]] <- sd(x[[cn]]) / sqrt(n) }
  o$lambda_pruned_q025 <- q(x$lambda_pruned, 0.025); o$lambda_pruned_q975 <- q(x$lambda_pruned, 0.975)
  o$frac_lambda_pruned_gt_1.05 <- mean(x$lambda_pruned > 1.05)
  for (L in names(CAUSAL)) { s <- sum(x[[paste0("p_", L)]] < ALPHA); ci <- wilson(s, n)
    o[[paste0("power_", L)]] <- s / n; o[[paste0("power_", L, "_l95")]] <- ci[1]; o[[paste0("power_", L, "_u95")]] <- ci[2]
    o[[paste0("mean_se_", L)]] <- mean(x[[paste0("se_", L)]]); o[[paste0("mean_log_or_", L)]] <- mean(x[[paste0("log_or_", L)]]) }
  o$mean_n_p_lt_5e8 <- mean(x$n_p_lt_5e8); o$frac_stable <- mean(x$stable); o }))
## paired differences in lambda against the scan without components
base0 <- tab[tab$n_pcs == 0, ]
for (cn in c("lambda_pruned", "lambda_pc12")) {
  dd <- sapply(KS, function(k) { x <- tab[tab$n_pcs == k, ]; d <- x[[cn]] - base0[[cn]][match(x$replicate, base0$replicate)]
    c(mean(d), sd(d) / sqrt(length(d))) })
  summ[[paste0(cn, "_minus_0pc_mean")]] <- dd[1, ]; summ[[paste0(cn, "_minus_0pc_se")]] <- dd[2, ] }
num <- sapply(summ, is.numeric); summ[num] <- lapply(summ[num], function(v) signif(v, 4))
N <- nrow(rule)
rs <- do.call(rbind, lapply(CANDIDATES, function(k) { s <- sum(rule$chosen_n_pcs == k); ci <- wilson(s, N)
  data.frame(chosen_n_pcs = k, replicates = s, fraction = s / N, l95 = signif(ci[1], 4), u95 = signif(ci[2], 4)) }))
extra <- data.frame(quantity = c("replicates", "scree_floor_e", "replicates_with_association_floor_above_0",
    "replicates_with_a_component_at_nominal_p_lt_0.05", "mean_components_at_nominal_p_lt_0.05",
    "scan_without_components_passes_stability", "k2_not_stable", "k2_fails_lambda", "k2_fails_lead", "k2_fails_top10",
    "pc12_variants", "pc12_min_r2", "pruned_variants"),
  value = c(N, e_floor, sum(rule$association_floor_a > 0), sum(rule$n_components_nominal_p_lt_0.05 > 0),
            signif(mean(rule$n_components_nominal_p_lt_0.05), 4), sum(rule$stable_0), sum(!rule$stable_2),
            sum(rule$k2_fails_lambda), sum(rule$k2_fails_lead), sum(rule$k2_fails_top10),
            length(pc12), signif(min(r2sum[pc12]), 4), length(r2sum)))
wr <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
tabw <- tab; numt <- sapply(tabw, is.numeric); tabw[numt] <- lapply(tabw[numt], function(v) signif(v, 6))
wr(tabw, "pdac_demo_08_npc_per_replicate.tsv"); wr(rule, "pdac_demo_08_npc_rule.tsv")
wr(summ, "pdac_demo_08_npc_summary.tsv"); wr(rs, "pdac_demo_08_npc_rule_summary.tsv"); wr(extra, "pdac_demo_08_npc_rule_counts.tsv")

R <- c("SECTION 8: NUMBER OF PRINCIPAL COMPONENTS OVER THE REPLICATES, EUROPEAN SCAN",
       sprintf("%d replicates; SEX, AGE and PC1..PCk; P < %g", N, ALPHA), "",
       "LAMBDA (mean, SD over replicates; SE of the mean)")
for (j in seq_len(nrow(summ))) with(summ[j, ], R <<- c(R, sprintf(
  "  %2d PCs  pruned set %.4f, %.4f; %.4f   all variants %.4f, %.4f   variants following PC1-PC2 %.4f, %.4f; %.4f   lambda > 1.05 in %.3f",
  n_pcs, lambda_pruned_mean, lambda_pruned_sd, lambda_pruned_se_of_mean, lambda_all_mean, lambda_all_sd,
  lambda_pc12_mean, lambda_pc12_sd, lambda_pc12_se_of_mean, frac_lambda_pruned_gt_1.05)))
R <- c(R, "", "REALISED POWER (95% CI)")
for (j in seq_len(nrow(summ))) R <- c(R, sprintf("  %2d PCs  %s", summ$n_pcs[j], paste(sapply(names(CAUSAL), function(L)
  sprintf("%s %.1f%% (%.1f to %.1f)", LABEL[[L]], 100 * summ[[paste0("power_", L)]][j], 100 * summ[[paste0("power_", L, "_l95")]][j],
          100 * summ[[paste0("power_", L, "_u95")]][j])), collapse = "   ")))
R <- c(R, "", sprintf("RULE OF 04A STEP 01b (candidates %s; scree floor e = %d in every replicate)", paste(CANDIDATES, collapse = ", "), e_floor))
for (j in seq_len(nrow(rs))) with(rs[j, ], R <<- c(R, sprintf("  chooses %2d in %3d replicates (%.3f, 95%% CI %.3f to %.3f)", chosen_n_pcs, replicates, fraction, l95, u95)))
R <- c(R, sprintf("  a component associated with case status at P < %.4f in %d replicates; at nominal P < 0.05 in %d (mean %.2f of 20 components)",
                  ALPHA_ASSOC, sum(rule$association_floor_a > 0), sum(rule$n_components_nominal_p_lt_0.05 > 0), mean(rule$n_components_nominal_p_lt_0.05)),
       sprintf("  two components not stable in %d replicates: lambda %d, lead variant %d, ten smallest P %d (a replicate can fail more than one)",
               sum(!rule$stable_2), sum(rule$k2_fails_lambda), sum(rule$k2_fails_lead), sum(rule$k2_fails_top10)),
       sprintf("  the scan without components would pass the stability criterion in %d replicates (zero is not a candidate of the rule)", sum(rule$stable_0)),
       "", sprintf("variants following PC1-PC2: %d of %d pruned variants, R2 on the two components at least %.4f", length(pc12), length(r2sum), min(r2sum[pc12])),
       if (exists("chk_line")) chk_line else "no check against step 02",
       sprintf("elapsed %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(R[-length(R)], file.path(OUT, "pdac_demo_08_npc_summary.txt")); cat(R, sep = "\n")
