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
#   It also checks the rule by which demo_data/sample_ancestry.tsv was
#   written when the dataset was built (the build script is recorded in
#   scripts/data_release/build_gwas_tutorial_dataset.sh): the label is the gnomAD
#   v3.1 inferred genetic-ancestry group of the sample, column
#   population_inference.pop of the HGDP+1KG metadata, with nfe and fin
#   merged into eur, afr kept as afr and eas as eas; samples with
#   high_quality true only. It is the group gnomAD inferred from the
#   genotypes, not the population or region the sample was collected under,
#   so a few individuals carry a label that differs from the region of their
#   population. The analysis sets are these labels restricted to the
#   individuals that pass quality control (02 step 04, 05 step 01).
#
# PRIVATE INPUT, NOT IN THE REPOSITORY AND NOT IN THE RELEASE
#   The mapping between the anonymised identifiers and the source identifiers
#   is withheld (demo_dataset/scripts/README.md). This step can therefore be
#   run only where that file exists. It reads two paths from the environment:
#     GWAS_PRIVATE_ID_MAP   tab-separated, columns original_IID and new_ID
#     GWAS_SOURCE_META      gnomAD v3.1 HGDP+1KG sample metadata (tab-separated;
#                           columns s, hgdp_tgp_meta.Project,
#                           hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Population,
#                           population_inference.pop, high_quality)
#   If either is unset or missing the step says so and writes nothing.
#   It writes counts per population only: no identifier, released or source,
#   is written or printed, and no row describes an individual.
#
# OUTPUT (sections/08_replicates/results/)
#   pdac_demo_08_set_composition.tsv              analysis set x source population: n
#   pdac_demo_08_set_composition_by_project.tsv   analysis set x panel: n
#   pdac_demo_08_label_rule.tsv                   released label x gnomAD inferred group x genetic region: n;
#                                                 all released individuals, and those in an analysis set
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
need_rule <- c("population_inference.pop", "high_quality")
if (!all(c("original_IID", "new_ID") %in% names(map)) || !all(c(need, need_rule) %in% names(meta))) stop("unexpected columns in the private inputs")
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
byp <- aggregate(n ~ set + set_label + set_n + project, tab, sum); byp <- byp[order(match(byp$set, names(sets)), byp$project), ]
byp$fraction <- signif(byp$n / byp$set_n, 4)
write.table(byp, file.path(OUT, "pdac_demo_08_set_composition_by_project.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## the rule of the released labels, on all released individuals
anc <- read.table("demo_data/sample_ancestry.tsv", col.names = c("IID", "label"), stringsAsFactors = FALSE)
ma <- meta[match(map$original_IID[match(anc$IID, map$new_ID)], meta$s), ]
if (anyNA(ma$s)) stop("released individuals not found in the private inputs")
inferred <- ma[["population_inference.pop"]]; by_rule <- ifelse(inferred %in% c("nfe", "fin"), "eur", inferred)
in_set <- anc$IID %in% unlist(lapply(sets, function(f) read.table(f, stringsAsFactors = FALSE)[[2]]))
rule <- as.data.frame(table(label = anc$label, gnomad_inferred_group = inferred, genetic_region = ma[["hgdp_tgp_meta.Genetic.region"]]), stringsAsFactors = FALSE)
names(rule)[4] <- "n_released"
rs <- as.data.frame(table(label = anc$label[in_set], gnomad_inferred_group = inferred[in_set], genetic_region = ma[["hgdp_tgp_meta.Genetic.region"]][in_set]), stringsAsFactors = FALSE)
rule$n_in_analysis_sets <- rs$Freq[match(paste(rule$label, rule$gnomad_inferred_group, rule$genetic_region), paste(rs$label, rs$gnomad_inferred_group, rs$genetic_region))]
rule$n_in_analysis_sets[is.na(rule$n_in_analysis_sets)] <- 0L
rule <- rule[rule$n_released > 0, ]; rule$region_matches_label <- toupper(rule$label) == rule$genetic_region
rule_ok <- sum(by_rule == anc$label); hq <- sum(ma$high_quality %in% c("true", "TRUE", TRUE))
write.table(rule, file.path(OUT, "pdac_demo_08_label_rule.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
R <- c("SECTION 8: SOURCE POPULATIONS OF THE ANALYSIS SETS (counts; labels of the gnomAD v3.1 HGDP+1KG metadata)", "")
for (g in names(sets)) { x <- tab[tab$set == g, ]
  R <- c(R, sprintf("%s: %d individuals, %d populations; by project: %s; by genetic region: %s", label[[g]], x$set_n[1], nrow(x),
                    paste(sprintf("%s %d", names(tapply(x$n, x$project, sum)), tapply(x$n, x$project, sum)), collapse = ", "),
                    paste(sprintf("%s %d", names(tapply(x$n, x$genetic_region, sum)), tapply(x$n, x$genetic_region, sum)), collapse = ", ")),
         paste0("  ", paste(strwrap(paste(sprintf("%s %d", x$population, x$n), collapse = ", "), 110), collapse = "\n  ")), "") }
R <- c(R, "RULE OF THE RELEASED LABELS (demo_data/sample_ancestry.tsv)",
       sprintf("  label = gnomAD population_inference.pop with nfe and fin merged into eur: holds for %d of %d released individuals; high_quality true for %d",
               rule_ok, nrow(anc), hq),
       "  label, gnomAD inferred group, genetic region of the population: released (in an analysis set)")
for (j in seq_len(nrow(rule))) R <- c(R, sprintf("  %s  %s  %s: %d (%d)%s", rule$label[j], rule$gnomad_inferred_group[j], rule$genetic_region[j],
  rule$n_released[j], rule$n_in_analysis_sets[j], if (rule$region_matches_label[j]) "" else "   region differs from the label"))
writeLines(R, file.path(OUT, "pdac_demo_08_set_composition.txt")); cat(R, sep = "\n")
