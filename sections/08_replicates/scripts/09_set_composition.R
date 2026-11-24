#!/usr/bin/env Rscript

################################################################################
# Section 8, step 09: source populations of the analysis sets (counts only)
#
# PURPOSE
#   The analysis sets are defined by genetic similarity (Section 2), and the
#   released identifiers are anonymised. This step reports from which source
#   populations of the HGDP and 1000 Genomes panels the individuals of each
#   analysis set come, as counts per population. It is the description of the
#   set every replicate is analysed in.
#
# PRIVATE INPUT, NOT IN THE REPOSITORY AND NOT IN THE RELEASE
#   The mapping between the anonymised identifiers and the source identifiers
#   is withheld (demo_dataset/scripts/README.md). This step can therefore be
#   run only where that file exists. It reads two paths from the environment:
#     GWAS_PRIVATE_ID_MAP   tab-separated, columns original_IID and new_ID
#     GWAS_SOURCE_META      gnomAD v3.1 HGDP+1KG sample metadata (tab-separated;
#                           columns s, hgdp_tgp_meta.Project,
#                           hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Population)
#   If either is unset or missing the step says so and writes nothing.
#   It writes counts per population only: no identifier, released or source,
#   is written or printed, and no row describes an individual.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_set_composition.tsv   analysis set x source population: n
#   pdac_demo_08_set_composition.txt   plain-text report
################################################################################

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(normalizePath(file.path(script_dir, "..", "..", "..")))
OUT <- "sections/08_replicates/results"
fm <- Sys.getenv("GWAS_PRIVATE_ID_MAP"); fs <- Sys.getenv("GWAS_SOURCE_META")
if (!nzchar(fm) || !nzchar(fs) || !file.exists(fm) || !file.exists(fs)) {
  cat("GWAS_PRIVATE_ID_MAP and GWAS_SOURCE_META are not both set to existing files: the composition cannot be computed here; nothing written\n")
  quit(status = 0)
}
map <- read.delim(fm, stringsAsFactors = FALSE)
meta <- read.delim(fs, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
need <- c("s", "hgdp_tgp_meta.Project", "hgdp_tgp_meta.Genetic.region", "hgdp_tgp_meta.Population")
if (!all(c("original_IID", "new_ID") %in% names(map)) || !all(need %in% names(meta))) stop("unexpected columns in the private inputs")
sets <- list(eur = "results/pca/pdac_demo_02_eur_keep.txt", afr = "results/meta/afr/keep.txt", eas = "results/meta/eas/keep.txt")
label <- c(eur = "European analysis set", afr = "African stratum", eas = "East Asian stratum")
tab <- do.call(rbind, lapply(names(sets), function(g) {
  ids <- read.table(sets[[g]], stringsAsFactors = FALSE)[[2]]
  src <- map$original_IID[match(ids, map$new_ID)]; m <- meta[match(src, meta$s), need[-1]]
  if (anyNA(src) || anyNA(m[[1]])) stop("individuals of ", g, " not found in the private inputs")
  x <- as.data.frame(table(project = m[[1]], genetic_region = m[[2]], population = m[[3]]), stringsAsFactors = FALSE)
  x <- x[x$Freq > 0, ]; x <- x[order(-x$Freq, x$population), ]
  data.frame(set = g, set_label = label[[g]], set_n = length(ids), project = x$project, genetic_region = x$genetic_region,
             population = x$population, n = x$Freq, fraction = signif(x$Freq / length(ids), 4)) }))
write.table(tab, file.path(OUT, "pdac_demo_08_set_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
R <- c("SECTION 8: SOURCE POPULATIONS OF THE ANALYSIS SETS (counts; labels of the gnomAD v3.1 HGDP+1KG metadata)", "")
for (g in names(sets)) { x <- tab[tab$set == g, ]
  R <- c(R, sprintf("%s: %d individuals, %d populations; by project: %s; by genetic region: %s", label[[g]], x$set_n[1], nrow(x),
                    paste(sprintf("%s %d", names(tapply(x$n, x$project, sum)), tapply(x$n, x$project, sum)), collapse = ", "),
                    paste(sprintf("%s %d", names(tapply(x$n, x$genetic_region, sum)), tapply(x$n, x$genetic_region, sum)), collapse = ", ")),
         paste0("  ", paste(strwrap(paste(sprintf("%s %d", x$population, x$n), collapse = ", "), 110), collapse = "\n  ")), "") }
writeLines(R, file.path(OUT, "pdac_demo_08_set_composition.txt")); cat(R, sep = "\n")
