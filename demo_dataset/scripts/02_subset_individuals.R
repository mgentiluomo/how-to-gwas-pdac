#!/usr/bin/env Rscript
# =============================================================================
# 02_subset_individuals.R
# -----------------------------------------------------------------------------
# Selects 6,000 synthetic individuals from HAPNEST with a defined genetic-
# ancestry composition, so that the demo cohort looks like a realistic, mostly
# European PDAC case-control study while still containing enough non-European
# samples to TEACH population stratification (Section 2).
#
#   target mix (n = 6,000):   EUR 80%  -> 4,800
#                             AFR 10%  ->   600
#                             EAS  5%  ->   300
#                             SAS  5%  ->   300
#
# At this stage we choose WHO is in the study. We do NOT yet know who is a case
# and who is a control — the phenotype is simulated later, in
# 03_simulate_phenotype.R. We also apply NO quality control here: QC is a
# teaching exercise in Section 1B (Murat Güler).
#
# Output:
#   data/subset/keep_ids.txt    individuals to keep (for PLINK --keep)
#   data/subset/sample_meta.tsv full metadata for the selected individuals
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---- command-line options ----------------------------------------------------
option_list <- list(
  make_option("--sample_file", type = "character",
              default = "./data/hapnest_raw/hapnest_samples.tsv",
              help = "HAPNEST sample description table [default %default]"),
  make_option("--out_dir", type = "character",
              default = "./data/subset",
              help = "output directory [default %default]"),
  make_option("--n_total", type = "integer", default = 6000L,
              help = "total individuals to select [default %default]"),
  make_option("--seed", type = "integer", default = 2026L,
              help = "random seed for reproducibility [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

# ---- desired ancestry composition --------------------------------------------
# Proportions must sum to 1. Edit here if you want a different teaching mix.
mix <- c(EUR = 0.80, AFR = 0.10, EAS = 0.05, SAS = 0.05)
stopifnot(abs(sum(mix) - 1) < 1e-9)

n_target <- round(mix * opt$n_total)
# fix any rounding drift so the totals add up exactly
n_target["EUR"] <- n_target["EUR"] + (opt$n_total - sum(n_target))

cat("Target ancestry composition:\n")
print(n_target)

# ---- read the HAPNEST sample table -------------------------------------------
if (!file.exists(opt$sample_file)) {
  stop("Sample file not found: ", opt$sample_file,
       "\nCheck the path / the file name used on EBI (see 01_download_hapnest.sh).")
}
meta <- fread(opt$sample_file)

# HAPNEST tables vary slightly between releases. We try to locate, robustly:
#   - an individual identifier column   (sample_id / IID / ID / id)
#   - a genetic-ancestry / population column (population / ancestry / pop / group)
id_col  <- intersect(c("sample_id", "IID", "ID", "id", "#IID"), names(meta))[1]
pop_col <- intersect(c("population", "ancestry", "pop", "group",
                       "genetic_ancestry", "superpop"), names(meta))[1]

if (is.na(id_col) || is.na(pop_col)) {
  cat("Columns found in the sample file:\n"); print(names(meta))
  stop("Could not auto-detect the ID and/or population columns. ",
       "Edit id_col / pop_col in this script to match your HAPNEST release.")
}
cat(sprintf("Using ID column '%s' and population column '%s'.\n",
            id_col, pop_col))

# standardise the two columns we need
meta <- meta[, .(sample_id = get(id_col), population = toupper(get(pop_col)))]

cat("\nAvailable individuals per population in HAPNEST:\n")
print(meta[, .N, by = population][order(-N)])

# ---- sample the requested number from each ancestry group --------------------
selected <- rbindlist(lapply(names(n_target), function(pop) {
  pool <- meta[population == pop]
  k    <- n_target[[pop]]
  if (nrow(pool) < k) {
    stop(sprintf("Only %d '%s' individuals available but %d requested.",
                 nrow(pool), pop, k))
  }
  pool[sample(.N, k)]
}))

# shuffle so the file is not ordered by ancestry
selected <- selected[sample(.N)]

cat(sprintf("\nSelected %d individuals:\n", nrow(selected)))
print(selected[, .N, by = population][order(-N)])

# ---- write the PLINK keep file -----------------------------------------------
# PLINK's --keep expects two whitespace-separated columns: FID and IID.
# HAPNEST uses the same identifier for both; if your .psam stores FID = 0,
# change the FID column below to 0 to match (see comment in 02b_apply_subset.sh).
keep_dt <- data.table(FID = selected$sample_id, IID = selected$sample_id)
fwrite(keep_dt, file.path(opt$out_dir, "keep_ids.txt"),
       sep = "\t", col.names = FALSE)

# keep the full metadata for the phenotype-simulation step
fwrite(selected, file.path(opt$out_dir, "sample_meta.tsv"), sep = "\t")

cat(sprintf("\nWrote %d individuals to %s\n",
            nrow(selected), file.path(opt$out_dir, "keep_ids.txt")))
cat("Next step:  bash scripts/02b_apply_subset.sh\n")
