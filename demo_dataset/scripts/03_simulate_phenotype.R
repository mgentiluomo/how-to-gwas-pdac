#!/usr/bin/env Rscript
# =============================================================================
# 03_simulate_phenotype.R
# -----------------------------------------------------------------------------
# Simulates a binary PDAC-like phenotype (1,000 cases + 5,000 controls) on top
# of the real LD/allele-frequency structure of the 6,000 HAPNEST individuals.
#
# HOW IT WORKS (liability-threshold model)
#   1. For each of six known PDAC regions we pick, from the genotype data, the
#      COMMON variant nearest the target position and call it "causal".
#   2. Each causal variant is given a realistic per-allele odds ratio
#      (OR 1.15-1.25). We work on the log scale: beta = log(OR).
#   3. We build a genetic liability  g = sum_j beta_j * (genotype_j)  and scale
#      it so that the genetic variance equals the target heritability h2.
#   4. We add Gaussian "environmental" noise for the remaining (1 - h2) variance
#      to obtain each individual's total liability L.
#   5. The 1,000 individuals with the HIGHEST liability become CASES; the
#      remaining 5,000 become CONTROLS. (This mimics ascertained case-control
#      sampling: cases are enriched relative to population prevalence.)
#
# Because the environmental noise dominates (1 - h2 = 0.85), each individual
# variant explains only a little of the outcome — exactly as in real GWAS — so
# most loci will NOT be genome-wide significant at n = 6,000. That lesson is the
# whole point (see README).
#
# Output:
#   data/subset/pheno.txt          FID IID PHENO   (1 = control, 2 = case)
#   data/subset/causal_truth.tsv   the variants we made causal + their OR
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---- command-line options ----------------------------------------------------
option_list <- list(
  make_option("--bfile", type = "character",
              default = "./data/subset/hapnest_6k",
              help = "PLINK1 prefix of the merged 6k dataset [default %default]"),
  make_option("--out_dir", type = "character", default = "./data/subset",
              help = "output directory [default %default]"),
  make_option("--n_cases", type = "integer", default = 1000L,
              help = "number of cases [default %default]"),
  make_option("--h2", type = "double", default = 0.15,
              help = "liability-scale heritability [default %default]"),
  make_option("--or_scale", type = "double", default = 1.0,
              help = "multiply all log-ORs by this factor (>1 = stronger, easier signal) [default %default]"),
  make_option("--window_kb", type = "integer", default = 100L,
              help = "half-window (kb) to search for a causal variant [default %default]"),
  make_option("--seed", type = "integer", default = 2026L,
              help = "random seed [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- the six target PDAC regions (GRCh38) ------------------------------------
# Per-allele OR are illustrative values within the realistic common-variant
# range for PDAC. Adjust if you update the underlying references.
targets <- data.table(
  chr   = c(1L, 5L, 9L, 13L, 16L, 17L),
  bp    = c(200000000, 1320000, 136149000, 73330000, 75180000, 73261000),
  gene  = c("NR5A2", "TERT/CLPTM1L", "ABO", "KLF5", "BCAR1", "LINC00673"),
  OR    = c(1.20, 1.25, 1.22, 1.15, 1.16, 1.18)
)

# ---- 1. read the .bim and choose one causal variant per region ---------------
bim <- fread(paste0(opt$bfile, ".bim"),
             col.names = c("chr", "snp", "cm", "bp", "a1", "a2"))

pick_causal <- function(t) {
  win <- opt$window_kb * 1000L
  cand <- bim[chr == t$chr & bp >= (t$bp - win) & bp <= (t$bp + win)]
  if (nrow(cand) == 0L)
    stop(sprintf("No variant within %d kb of %s (chr%d:%d). ",
                 opt$window_kb, t$gene, t$chr, t$bp),
         "Widen --window_kb or check the build.")
  # nearest to the target position
  cand[which.min(abs(bp - t$bp))]
}
causal <- rbindlist(lapply(seq_len(nrow(targets)),
                           function(i) pick_causal(targets[i])))
causal[, `:=`(gene = targets$gene, OR = targets$OR,
              beta = log(targets$OR) * opt$or_scale)]

cat("Chosen causal variants:\n")
print(causal[, .(chr, snp, bp, gene, OR)])

# ---- 2. extract genotype dosages for the causal variants via PLINK2 ----------
# We export an additive (0/1/2) coding for just these six variants.
extract_file <- file.path(opt$out_dir, "causal_snps.txt")
fwrite(causal[, .(snp)], extract_file, col.names = FALSE)

raw_prefix <- file.path(opt$out_dir, "causal_dosages")
cmd <- sprintf(
  "plink2 --bfile %s --extract %s --export A --out %s",
  shQuote(opt$bfile), shQuote(extract_file), shQuote(raw_prefix))
cat("\nRunning:\n  ", cmd, "\n")
status <- system(cmd)
if (status != 0L) stop("PLINK2 export failed (is plink2 on your PATH?).")

# PLINK2 --export A writes a .raw file: columns FID IID PAT MAT SEX PHENOTYPE
# followed by one column per variant named "<snp>_<allele>".
raw <- fread(paste0(raw_prefix, ".raw"))
geno_cols <- setdiff(names(raw),
                     c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
G <- as.matrix(raw[, ..geno_cols])
# mean-impute the (rare) missing genotype calls
for (j in seq_len(ncol(G))) {
  m <- is.na(G[, j]); if (any(m)) G[m, j] <- mean(G[, j], na.rm = TRUE)
}

# ---- 3. build the genetic liability ------------------------------------------
# Align the beta vector to the column order of G (columns are "<snp>_<allele>").
snp_of_col <- sub("_[ACGT]+$", "", geno_cols)
beta <- causal$beta[match(snp_of_col, causal$snp)]

g <- as.numeric(G %*% beta)               # raw genetic liability
g <- scale(g)[, 1]                        # standardise to mean 0, var 1

# ---- 4. add environmental noise to hit the target heritability ---------------
# L = sqrt(h2) * g_std + sqrt(1 - h2) * e,   e ~ N(0,1)
e <- rnorm(length(g))
L <- sqrt(opt$h2) * g + sqrt(1 - opt$h2) * e

# ---- 5. top-liability individuals become cases -------------------------------
n_cases <- opt$n_cases
ord  <- order(L, decreasing = TRUE)
case <- rep(0L, length(L)); case[ord[seq_len(n_cases)]] <- 1L

pheno <- data.table(
  FID   = raw$FID,
  IID   = raw$IID,
  PHENO = ifelse(case == 1L, 2L, 1L)      # PLINK coding: 1 = control, 2 = case
)

cat(sprintf("\nPhenotype: %d cases / %d controls (n = %d).\n",
            sum(pheno$PHENO == 2L), sum(pheno$PHENO == 1L), nrow(pheno)))

# ---- write outputs -----------------------------------------------------------
fwrite(pheno, file.path(opt$out_dir, "pheno.txt"), sep = "\t")
fwrite(causal[, .(chr, snp, bp, gene, OR, beta)],
       file.path(opt$out_dir, "causal_truth.tsv"), sep = "\t")

cat("\nWrote pheno.txt and causal_truth.tsv to", opt$out_dir, "\n")
cat("Keep causal_truth.tsv private from students during the exercise —\n",
    "it is the answer key for which loci they should rediscover.\n", sep = "")
cat("\nNext step:  bash scripts/04_run_gwas.sh\n")
