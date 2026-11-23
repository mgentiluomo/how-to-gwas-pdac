#!/usr/bin/env Rscript

################################################################################
# Section 6, step 04: conditional analysis at ABO on individual genotypes
#
# PURPOSE
#   Ask how many independent signals the ABO region contains. With individual
#   genotypes this is an exact regression: the lead variant of the primary
#   scan is added as a covariate (plink2 --condition) and the region is tested
#   again. A residual signal at genome-wide significance would mean a second
#   independent association; none means the region is explained by one.
#
# MODEL
#   The primary scan of Section 4A in the same European analysis set: logistic
#   regression with firth-fallback, SEX, AGE and the number of components
#   chosen by 04A step 01b, plus the lead variant. Region: the fine-mapping
#   region of step 01, chr9:132,900,000-133,600,000 (GRCh38).
#   For teaching, the same is repeated conditioning on the simulated causal
#   variant (demo_data/truth.tsv), which the analyst would not know.
#
# INPUT
#   results/pca/pdac_demo_02_hwe_filt.*, results/pca/pdac_demo_02_eur_keep.txt
#   results/assoc/pdac_demo_04A_covar.txt, pdac_demo_04A_npc_chosen.txt,
#   results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid
#   demo_data/phenotype.txt, demo_data/truth.tsv
#
# OUTPUT (results/finemap/, collected into sections/06_finemapping_annotation/results/)
#   pdac_demo_06_conditional.tsv        one row per analysis (none / lead / causal)
#   pdac_demo_06_conditional_top.tsv    ten smallest P in the region, per analysis
#   pdac_demo_06_conditional.txt        plain-text report
################################################################################

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

IN <- "results/pca/pdac_demo_02_hwe_filt"; KEEP <- "results/pca/pdac_demo_02_eur_keep.txt"
COVAR <- "results/assoc/pdac_demo_04A_covar.txt"; NPCF <- "results/assoc/pdac_demo_04A_npc_chosen.txt"
GWAS <- "results/assoc/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid"
PHENO <- "demo_data/phenotype.txt"; OUT_DIR <- "results/finemap"
CHR <- 9; FROM <- 132900000; TO <- 133600000
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (f in c(paste0(IN, ".bed"), KEEP, COVAR, NPCF, GWAS, PHENO, "demo_data/truth.tsv"))
  if (!file.exists(f)) stop("missing input: ", f)
TMP <- file.path(tempdir(), "cond"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) { d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1]); d }
k <- as.integer(read.table(NPCF)[1, 2])

g <- rd(GWAS); if ("TEST" %in% names(g)) g <- g[g$TEST == "ADD", ]
reg <- g[g$CHROM == CHR & g$POS >= FROM & g$POS <= TO & is.finite(g$P), ]
lead <- reg$ID[which.min(reg$P)]
truth <- read.delim("demo_data/truth.tsv", stringsAsFactors = FALSE)
causal <- truth$variant[grepl("^9:", truth$variant)][1]

run <- function(tag, cond) {
  out <- file.path(TMP, tag)
  cmd <- paste("plink2 --bfile", IN, "--chr", CHR, "--from-bp", FROM, "--to-bp", TO, "--keep", KEEP,
               "--pheno", PHENO, "--covar", COVAR,
               "--covar-name", paste0("SEX,AGE,", paste0("PC", 1:k, collapse = ",")),
               "--covar-variance-standardize",
               if (!is.na(cond)) paste("--condition", cond) else "",
               "--glm firth-fallback hide-covar cols=+a1freq,+beta,+orbeta,+nobs,+err --out", out, "> /dev/null 2>&1")
  if (system(cmd) != 0) stop("plink2 failed: ", cmd)
  r <- rd(Sys.glob(paste0(out, ".*.glm.logistic.hybrid"))[1])
  if ("TEST" %in% names(r)) r <- r[r$TEST == "ADD", ]
  r <- r[is.na(cond) | r$ID != cond, ]
  r[order(r$P, na.last = TRUE), ]
}
# P of one variant, or its PLINK error code when the fit returned no P
# (for example a variant almost collinear with the conditioning variant)
p_of <- function(r, id) { x <- r[r$ID == id, ]
  if (!nrow(x)) return("absent")
  if (is.finite(x$P[1])) format(signif(x$P[1], 3)) else paste0("NA (", x$ERRCODE[1], ")") }
res <- list(none = run("none", NA), lead = run("lead", lead), causal = run("causal", causal))
summ <- do.call(rbind, lapply(names(res), function(t) {
  r <- res[[t]]; ok <- is.finite(r$P) & r$P > 0
  data.frame(conditioned_on = if (t == "none") "none" else if (t == "lead") paste("lead", lead) else paste("causal", causal),
             n_pcs = k, variants_tested = sum(ok), variants_no_p = sum(!ok),
             min_p = signif(r$P[1], 3), min_p_id = r$ID[1],
             n_p_lt_5e8 = sum(r$P[ok] < 5e-8), n_p_lt_1e5 = sum(r$P[ok] < 1e-5), n_p_lt_0.001 = sum(r$P[ok] < 1e-3),
             causal_p = if (t == "causal") NA else p_of(r, causal),
             lead_p = if (t == "lead") NA else p_of(r, lead))
}))
top <- do.call(rbind, lapply(names(res), function(t) data.frame(analysis = t, head(res[[t]][is.finite(res[[t]]$P), c("ID", "POS", "A1", "BETA", "SE", "P")], 10))))
write.table(summ, file.path(OUT_DIR, "pdac_demo_06_conditional.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top, file.path(OUT_DIR, "pdac_demo_06_conditional_top.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rep <- c("SECTION 6: CONDITIONAL ANALYSIS AT ABO (individual genotypes, plink2 --condition)",
         sprintf("region chr%d:%d-%d, European analysis set, SEX + AGE + %d PCs", CHR, FROM, TO, k), "",
         sprintf("%-40s %8s %10s %-18s %6s %6s %8s", "conditioned on", "variants", "min P", "at", "P<5e-8", "P<1e-5", "P<0.001"),
         sapply(seq_len(nrow(summ)), function(i) with(summ[i, ], sprintf("%-40s %8d %10.3g %-18s %6d %6d %8d",
           conditioned_on, variants_tested, min_p, min_p_id, n_p_lt_5e8, n_p_lt_1e5, n_p_lt_0.001))), "",
         sprintf("P of the causal variant after conditioning on the lead: %s; P of the lead after conditioning on the causal: %s",
                 summ$causal_p[summ$conditioned_on != "none" & grepl("^lead", summ$conditioned_on)],
                 summ$lead_p[grepl("^causal", summ$conditioned_on)]), "",
         "HOW TO READ THIS",
         "  No variant at genome-wide significance after conditioning on the lead means",
         "  the region holds one detectable signal in this sample. The simulation has one",
         "  causal variant here, so a residual signal would be a false second signal.")
writeLines(rep, file.path(OUT_DIR, "pdac_demo_06_conditional.txt")); cat(rep, sep = "\n")
