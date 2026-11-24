#!/usr/bin/env Rscript

################################################################################
# write_facts_tsv.R
#
# Writes one flat table of the figures of a canonical run: key, value, source.
# Every key is unique (the script stops otherwise), every value is read from a
# results file written by the run, and the source column names that file, so a
# manuscript or a page can quote a number by key and anyone can trace it back.
#
# USAGE (from the repository root, after run_all.sh and before or after the
# results are committed):
#   Rscript scripts/dev/write_facts_tsv.R [OUT]
#   OUT defaults to ../facts_<date>.tsv, <date> being the "finished" date in
#   RUN_MANIFEST.txt (YYYYMMDD).
#
# KEYS
#   <section>.<table>.<row>.<column> for summary tables read whole, lower case,
#   non-alphanumerics replaced by "_". A few keys are computed explicitly
#   (causal variants, credible sets, survival test calibration, truth table,
#   run provenance); their names start with the section they belong to.
#   Values are written as they appear in the source; nothing is rounded here.
################################################################################

suppressPackageStartupMessages(library(utils))
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
here <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
root <- Sys.getenv("FACTS_ROOT", normalizePath(file.path(here, "..", "..")))
setwd(root)
# FACTS_ALLOW_MISSING=1 is for comparisons with an earlier run checked out from
# git, which lacks the intermediate files in results/: keys that need them are
# skipped with a message. A canonical facts table is never written with it.
ALLOW_MISSING <- Sys.getenv("FACTS_ALLOW_MISSING") == "1"
block <- function(expr) if (ALLOW_MISSING) tryCatch(expr, error = function(e) message("skipped: ", conditionMessage(e))) else expr

if (!file.exists("RUN_MANIFEST.txt")) stop("RUN_MANIFEST.txt not found: run scripts/dev/run_all.sh first")
man <- readLines("RUN_MANIFEST.txt")
if (any(grepl("NOT A CANONICAL RUN", man))) stop("RUN_MANIFEST.txt describes a non-canonical run")
mget1 <- function(pat) { x <- grep(pat, man, value = TRUE); if (!length(x)) NA else trimws(sub(pat, "", x[1])) }
finished <- mget1("^finished +")
out <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out)) out <- file.path("..", sprintf("facts_%s.tsv", gsub("-", "", substr(finished, 1, 10))))

rows <- list()
add <- function(key, value, source) {
  rows[[length(rows) + 1L]] <<- data.frame(key = key, value = as.character(value), source = source,
                                           stringsAsFactors = FALSE)
}
san <- function(x) gsub("^_+|_+$", "", gsub("[^a-z0-9.]+", "_", tolower(x)))
rd <- function(f, ...) { d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", ...)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
S <- function(sec, f) file.path("sections", sec, "results", f)
need <- function(f) if (!file.exists(f)) stop("missing results file: ", f)

## ---- provenance --------------------------------------------------------------------
add("run.commit", mget1("^commit +"), "RUN_MANIFEST.txt")
add("run.started", mget1("^started +"), "RUN_MANIFEST.txt")
add("run.finished", finished, "RUN_MANIFEST.txt")
add("run.data_release", mget1("^data release +"), "RUN_MANIFEST.txt")
add("run.steps", mget1("^steps +"), "RUN_MANIFEST.txt")
pl <- mget1("^4B permutation null +")
add("run.4b_n_perm", sub(".*n_perm +([0-9]+).*", "\\1", pl), "RUN_MANIFEST.txt")
add("run.4b_perm_seed", sub(".*seed +([0-9]+).*", "\\1", pl), "RUN_MANIFEST.txt")
add("run.susier", mget1("^susieR "), "RUN_MANIFEST.txt")
add("run.4a_components_chosen", sub(" .*", "", mget1("^4A components chosen +")), "RUN_MANIFEST.txt")
add("run.lrld_bed_md5", sub(".*md5 ([0-9a-f]{32}).*", "\\1", mget1("^long-range LD regions ")), "RUN_MANIFEST.txt")
add("run.plink2_version", grep("^PLINK v2", man, value = TRUE)[1], "RUN_MANIFEST.txt")
add("run.plink19_version", grep("^PLINK v1", man, value = TRUE)[1], "RUN_MANIFEST.txt")
add("run.regenie_version", sub("^REGENIE ", "", grep("^REGENIE", man, value = TRUE)[1]), "RUN_MANIFEST.txt")
add("run.r_version", grep("^R version", man, value = TRUE)[1], "RUN_MANIFEST.txt")
for (l in grep("^  [0-9a-f]{32}  ", man, value = TRUE)) {
  p <- strsplit(trimws(l), " +")[[1]]; add(paste0("run.md5.", san(p[2])), p[1], "RUN_MANIFEST.txt")
}

## ---- summary tables read whole ---------------------------------------------------------
# kv: two columns, name and value.  rows: first column (or `id`) identifies the row.
tables <- list(
  list("qc",   S("01B_genotyping_qc", "pdac_demo_09_qc_counts.tsv"), "rows", "step"),
  list("qc",   S("01B_genotyping_qc", "pdac_demo_03_sexcheck_summary.tsv"), "kv"),
  list("qc",   S("01B_genotyping_qc", "pdac_demo_04_het_summary.tsv"), "rows", "group"),
  list("qc",   S("01B_genotyping_qc", "pdac_demo_06_hwe_threshold_sweep.tsv"), "rows", "threshold"),
  list("qc",   S("01B_genotyping_qc", "pdac_demo_07_relatedness_summary.tsv"), "rows", "relationship_category"),
  list("qc",   S("01B_genotyping_qc", "pdac_demo_07_removal_summary.tsv"), "kv"),
  list("qc",   S("01B_genotyping_qc", "pdac_demo_07_relatedness_pruning_diagnostics.tsv"), "kv"),
  list("pca",  S("02_population_stratification", "pdac_demo_02_hwe_summary.tsv"), "rows", "group"),
  list("pca",  S("02_population_stratification", "pdac_demo_02_analysis_set.tsv"), "kv"),
  list("pca",  S("02_population_stratification", "pdac_demo_02_pca_variance.tsv"), "rows", NA),
  list("power", S("03B_statistical_power", "pdac_demo_07_power_table.tsv"), "rows", "MAF"),
  list("assoc", S("04A_association_binary", "pdac_demo_04A_lambda.tsv"), "kv"),
  list("assoc", S("04A_association_binary", "pdac_demo_04A_npc_scans.tsv"), "rows", "n_pcs"),
  list("assoc", S("04A_association_binary", "pdac_demo_04A_npc_criteria.tsv"), "kv"),
  list("assoc", S("04A_association_binary", "pdac_demo_04A_pc_model_comparison.tsv"), "rows", NA),
  list("assoc", S("04A_association_binary", "pdac_demo_04A_sens_relatedness.tsv"), "rows", "rule"),
  list("assoc", S("04A_association_binary", "pdac_demo_04A_sens_prune_r2.tsv"), "rows", "r2"),
  list("assoc", S("04A_association_binary", "pdac_demo_04A_sens_validation.tsv"), "rows", "quantity"),
  list("surv", S("04B_association_survival", "pdac_demo_04B_lambda.tsv"), "rows", NA),
  list("surv", S("04B_association_survival", "pdac_demo_04B_covariate_cost.tsv"), "rows", "covariates"),
  list("surv", S("04B_association_survival", "pdac_demo_04B_pc_comparison.tsv"), "rows", NA),
  list("surv", S("04B_association_survival", "pdac_demo_04B_perm_summary.tsv"), "rows", NA),
  list("surv", S("04B_association_survival", "pdac_demo_04B_perm_validation.tsv"), "rows", NA),
  list("surv", S("04B_association_survival", "pdac_demo_04B_ancestry_dependence.tsv"), "rows", "test"),
  list("surv", S("04B_association_survival", "pdac_demo_04B_ancestry_dependence_groups.tsv"), "rows", "group"),
  list("regenie", S("04C_association_other", "pdac_demo_regenie_comparison.tsv"), "rows", "analysis"),
  list("meta", S("05_meta_analysis", "pdac_demo_05_strata_counts.tsv"), "rows", "stratum"),
  list("finemap", S("06_finemapping_annotation", "pdac_demo_06_conditional.tsv"), "rows", "conditioned_on"),
  list("finemap", S("06_finemapping_annotation", "pdac_demo_06_susie_summary.tsv"), "kv", optional = TRUE),
  list("finemap", S("06_finemapping_annotation", "pdac_demo_06_susie_cs.tsv"), "rows", "cs", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_power.tsv"), "rows", c("locus", "analysis"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_winners_curse.tsv"), "rows", c("locus", "analysis"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_finemap.tsv"), "rows", c("locus", "analysis"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_i2.tsv"), "rows", "locus", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_lambda_summary.tsv"), "rows", "lambda", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_expected_se.tsv"), "rows", "locus", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_genotype_r.tsv"), "rows", c("stratum", "pair"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_tert_adjusted.tsv"), "rows", c("data", "model"), optional = TRUE),
  # Section 8 steps 05 to 09 and 4B step 06: outside the canonical run, on its outputs
  list("rep", S("08_replicates", "pdac_demo_08_tests_power.tsv"), "rows", c("locus", "analysis", "test"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_tests_null.tsv"), "rows", "test", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_npc_summary.tsv"), "rows", "n_pcs", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_randcov_summary.tsv"), "rows", "k", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_npc_rule_summary.tsv"), "rows", "chosen_n_pcs", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_npc_rule_counts.tsv"), "kv", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_causal_structure.tsv"), "rows", c("locus", "n_pcs"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_causal_partial_r.tsv"), "rows", c("pair", "given_n_pcs"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_expected_gain.tsv"), "rows", "comparison", optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_set_composition.tsv"), "rows", c("set", "population"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_simulated_detectable_or.tsv"), "rows", c("risk_allele", "test"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_set_composition_by_project.tsv"), "rows", c("set", "project"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_label_rule.tsv"), "rows", c("label", "gnomad_inferred_group", "genetic_region"), optional = TRUE),
  list("rep", S("08_replicates", "pdac_demo_08_plink2_variant_set.tsv"), "rows", "fit", optional = TRUE),
  list("surv", S("04B_association_survival", "pdac_demo_04B_randcov_summary.tsv"), "rows", "model", optional = TRUE),
  list("surv", S("04B_association_survival", "pdac_demo_04B_randcov_paired.tsv"), "rows", "contrast", optional = TRUE),
  list("surv", S("04B_association_survival", "pdac_demo_04B_ties.tsv"), "kv", optional = TRUE)
)
for (t in tables) {
  sec <- t[[1]]; f <- t[[2]]
  if (!file.exists(f) && (isTRUE(t$optional) || ALLOW_MISSING)) { message("table absent, no keys: ", f); next }
  need(f)
  d <- rd(f); stem <- san(sub("^pdac_demo_", "", sub("\\.tsv$", "", basename(f))))
  if (t[[3]] == "kv") {
    for (i in seq_len(nrow(d))) add(paste(sec, stem, san(d[i, 1]), sep = "."), d[i, 2], f)
  } else {
    idc <- if (is.na(t[[4]][1])) names(d)[1] else t[[4]]
    for (i in seq_len(nrow(d))) for (cn in setdiff(names(d), idc))
      add(paste(sec, stem, paste(san(unlist(d[i, idc])), collapse = "."), san(cn), sep = "."), d[i, cn], f)
  }
}

## ---- explicit keys ---------------------------------------------------------------------------
truth <- rd("demo_data/truth.tsv")
for (i in seq_len(nrow(truth))) {
  k <- san(if (grepl(":", truth$variant[i])) sub(" .*", "", truth$locus[i]) else truth$variant[i])
  add(paste0("truth.", k, ".generative_or"), truth$generative_OR[i], "demo_data/truth.tsv")
  add(paste0("truth.", k, ".realised_or"), truth$estimated_OR[i], "demo_data/truth.tsv")
}
causal <- c(abo = "9:133273682:A:T", tert = "5:1286401:C:A")

block({
# 04A primary scan at the two causal variants (risk allele = reference at both)
gf <- Sys.glob("results/assoc/pdac_demo_04A_gwas.*.glm.logistic.hybrid")[1]; need(gf)
g <- rd(gf); if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
for (nm in names(causal)) {
  r <- g[g$ID == causal[[nm]], ]; if (nrow(r) != 1) stop("causal variant not in 04A scan: ", causal[[nm]])
  ref <- strsplit(causal[[nm]], ":")[[1]][3]
  s <- if (r$A1 == ref) 1 else -1
  add(paste0("assoc.causal.", nm, ".a1"), r$A1, gf)
  add(paste0("assoc.causal.", nm, ".or_risk_allele"), exp(s * r$BETA), gf)
  add(paste0("assoc.causal.", nm, ".or_risk_l95"), exp(s * r$BETA - 1.959964 * r$SE), gf)
  add(paste0("assoc.causal.", nm, ".or_risk_u95"), exp(s * r$BETA + 1.959964 * r$SE), gf)
  add(paste0("assoc.causal.", nm, ".p"), r$P, gf)
}
th <- rd(S("04A_association_binary", "pdac_demo_04A_top_hits.tsv"))
th <- th[order(th$P), ]
add("assoc.lead.id", th$ID[1], S("04A_association_binary", "pdac_demo_04A_top_hits.tsv"))
add("assoc.lead.p", th$P[1], S("04A_association_binary", "pdac_demo_04A_top_hits.tsv"))
add("assoc.n_p_lt_5e8", sum(g$P < 5e-8, na.rm = TRUE), gf)
add("assoc.n_p_lt_1e5", sum(g$P < 1e-5, na.rm = TRUE), gf)
})

block({
# meta-analysis at the causal variants and genome-wide heterogeneity
mf <- "results/meta/pdac_demo_05_meta.tsv.gz"; need(mf)
m <- read.delim(gzfile(mf), check.names = FALSE, stringsAsFactors = FALSE)
for (nm in names(causal)) {
  r <- m[m$ID == causal[[nm]], ]; if (nrow(r) != 1) stop("causal variant not in meta: ", causal[[nm]])
  for (cn in intersect(names(r), c("A1", "n_strata", "BETA", "SE", "P", "Q", "P_het", "I2", "BETA_re", "SE_re", "P_re",
                                   paste0(rep(c("BETA_", "SE_", "P_", "FREQ_"), each = 3), c("eur", "afr", "eas")))))
    add(paste0("meta.causal.", nm, ".", san(cn)), r[[cn]], mf)
  sg <- if (r$A1 == strsplit(causal[[nm]], ":")[[1]][3]) 1 else -1   # orient to the risk (reference) allele
  orr <- function(b, se) c(exp(sg * b), exp(sg * b - 1.959964 * se), exp(sg * b + 1.959964 * se))
  for (v in list(c("fixed", "BETA", "SE"), c("random", "BETA_re", "SE_re"), c("eur", "BETA_eur", "SE_eur"),
                 c("afr", "BETA_afr", "SE_afr"), c("eas", "BETA_eas", "SE_eas"))) {
    o <- orr(r[[v[2]]], r[[v[3]]])
    add(paste0("meta.causal.", nm, ".or_risk_", v[1]), o[1], mf)
    add(paste0("meta.causal.", nm, ".or_risk_", v[1], "_l95"), o[2], mf)
    add(paste0("meta.causal.", nm, ".or_risk_", v[1], "_u95"), o[3], mf)
  }
}
add("meta.variants", nrow(m), mf)
add("meta.n_p_lt_5e8", sum(m$P < 5e-8, na.rm = TRUE), mf)
add("meta.n_i2_gt_50", sum(m$I2 > 50, na.rm = TRUE), mf)
add("meta.n_phet_lt_0.05", sum(m$P_het < 0.05, na.rm = TRUE), mf)
add("meta.lambda", median(qchisq(m$P[is.finite(m$P) & m$P > 0], 1, lower.tail = FALSE)) / qchisq(0.5, 1), mf)
})
# lambda per stratum and of the meta-analysis as the meta-analysis step printed them
local({
  sf <- S("05_meta_analysis", "pdac_demo_05_meta_summary.txt")
  if (!file.exists(sf) && ALLOW_MISSING) return(invisible())
  need(sf)
  l <- regmatches(readLines(sf), regexec("^\\s+(eur|afr|eas|meta)\\s+lambda = ([0-9.]+)", readLines(sf)))
  for (x in l[lengths(l) == 3]) add(paste0("meta.printed_lambda.", x[2]), x[3], sf)
})

block({
# within-European PCA: variance share of each extracted component
ef <- "results/pca/pdac_demo_02_pca_eur.eigenval"; need(ef)
ev <- scan(ef, quiet = TRUE)
for (i in seq_along(ev)) add(sprintf("pca.eur.pc%d.variance_explained_pct_of_extracted", i), 100 * ev[i] / sum(ev), ef)
})

# fine-mapping: credible sets and the rank of the causal variant
for (set in c("eur", "meta")) {
  f <- S("06_finemapping_annotation", sprintf("pdac_demo_06_pp_%s.tsv", set)); need(f)
  d <- rd(f); idc <- intersect(c("ID", "SNP", "variant"), names(d))[1]
  ppc <- intersect(c("PP", "pp", "posterior"), names(d))[1]
  d <- d[order(-d[[ppc]]), ]; cs <- cumsum(d[[ppc]])
  n_cs <- which(cs >= 0.95)[1]
  add(paste0("finemap.", set, ".variants"), nrow(d), f)
  add(paste0("finemap.", set, ".credible_set_size"), n_cs, f)
  add(paste0("finemap.", set, ".top_id"), d[[idc]][1], f)
  add(paste0("finemap.", set, ".top_pp"), d[[ppc]][1], f)
  k <- which(d[[idc]] == causal[["abo"]])
  add(paste0("finemap.", set, ".causal_rank"), k, f)
  add(paste0("finemap.", set, ".causal_pp"), d[[ppc]][k], f)
  add(paste0("finemap.", set, ".causal_in_set"), k <= n_cs, f)
}

# survival test calibration on the 5,000-variant subsample (lambda per statistic)
cf <- S("04B_association_survival", "pdac_demo_04B_calibration.tsv"); need(cf)
cal <- rd(cf)
for (cn in grep("^P_", names(cal), value = TRUE)) {
  p <- cal[[cn]]; p <- p[is.finite(p) & p > 0]
  add(paste0("surv.calibration.lambda_", san(sub("^P_", "", cn))), median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1), cf)
}

## ---- write ----------------------------------------------------------------------------------------
out_d <- do.call(rbind, rows)
dup <- out_d$key[duplicated(out_d$key)]
if (length(dup)) stop("duplicate keys: ", paste(unique(dup), collapse = ", "))
if (ALLOW_MISSING) message("FACTS_ALLOW_MISSING=1: comparison table only, not a canonical facts table")
if (any(is.na(out_d$value) | out_d$value == "NA"))
  warning(sum(is.na(out_d$value) | out_d$value == "NA"), " keys have NA values; see the table")
write.table(out_d, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s: %d keys, all unique\n", out, nrow(out_d)))
