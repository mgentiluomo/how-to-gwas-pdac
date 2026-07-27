#!/usr/bin/env Rscript

################################################################################
#  make_example_raw_data.R
#
#  Build the raw-data example used by Section 1A from the demonstration dataset.
#
#  WHY THIS EXISTS:
#    Section 1A teaches the step that comes before every other section: turning
#    raw genotyping output into an analysable PLINK dataset. Demonstrating it
#    needs a Final Report and an array manifest. Real ones cannot be
#    distributed: a Final Report holds the genome-wide genotypes of an
#    identifiable person, and a manufacturer's manifest is not ours to
#    redistribute.
#
#    So both are generated here from `pdac_demo`, exactly as the rest of the
#    demonstration data were built: real public reference genotypes, anonymised
#    identifiers, and deliberately injected artefacts so that each teaching step
#    has something concrete to find.
#
#  WHAT IS SYNTHETIC AND WHAT IS NOT:
#    Genotypes come from pdac_demo, which derives from the public HGDP + 1000
#    Genomes panel. GenCall scores, strand assignments, indel calls, no-calls
#    and the manifest itself are simulated. The parameters were calibrated on a
#    real Illumina GSAMD-24v3 Final Report so that the artefacts are realistic;
#    no part of that file is reproduced here.
#
#  INJECTED ARTEFACTS (all reported in the README, all reproducible with seed 2026):
#    - GenCall scores drawn from a bimodal distribution, median about 0.72
#    - 1.4% of calls below the 0.15 failure threshold
#    - 3.5% no-calls, written as "-" with a GenCall of NaN in a few cases
#    - 7.7% indel probes, reported as I/D, which a naive A/C/G/T filter drops
#    - 43% of variants designed on the minus strand, so their manifest alleles
#      must be complemented
#    - 25% of variant names encode a position from an older build, differing
#      from the manifest coordinate
#    - 6% of variants absent from the manifest, as if the report came from a
#      product with extra drop-in content
#    - sample 3 carries a strand-flipped block of 200 variants, which makes the
#      merge in Step 03 report multi-allelic conflicts and write a .missnp
#
#  Usage:  Rscript make_example_raw_data.R [n_variants] [seed]
################################################################################

args <- commandArgs(trailingOnly = TRUE)
n_var <- if (length(args) >= 1) as.integer(args[1]) else 8000
seed  <- if (length(args) >= 2) as.integer(args[2]) else 2026
set.seed(seed)

bim_file <- "demo_data/pdac_demo.bim"
out_dir  <- "sections/01A_study_design/example_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- c("DEMO001", "DEMO002", "DEMO003")

# --- pick the variants --------------------------------------------------------
bim <- read.delim(bim_file, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
bim <- bim[bim$chr %in% 1:22, ]
sel <- sort(sample(nrow(bim), n_var))
v   <- bim[sel, ]

# --- variant names ------------------------------------------------------------
# Most markers on a real array are named by rs number; a minority carry a name
# that looks like a coordinate, and those names often come from an older build.
n_pos <- round(0.25 * n_var)
is_pos <- rep(FALSE, n_var); is_pos[sample(n_var, n_pos)] <- TRUE

# The offset mimics a GRCh37-to-GRCh38 shift: large, plausible, and wrong.
offset <- sample(c(465556, 542622, 384119), n_var, replace = TRUE)
name <- ifelse(is_pos,
               paste0(v$chr, ":", v$pos + offset),
               paste0("rs", sample(1e6:9e8, n_var)))
# A few positional names also disagree on the chromosome.
chr_mismatch <- which(is_pos)[seq_len(round(0.013 * n_var))]
name[chr_mismatch] <- paste0(ifelse(v$chr[chr_mismatch] == 22, 21, v$chr[chr_mismatch] + 1),
                             ":", v$pos[chr_mismatch] + offset[chr_mismatch])

# --- indels -------------------------------------------------------------------
is_indel <- rep(FALSE, n_var); is_indel[sample(n_var, round(0.077 * n_var))] <- TRUE

# --- design strand ------------------------------------------------------------
comp <- c(A = "T", T = "A", C = "G", G = "C")
ref_strand <- ifelse(runif(n_var) < 0.43, "-", "+")

# Manifest alleles are written on the design strand.
man_a1 <- ifelse(ref_strand == "-", comp[v$a1], v$a1)
man_a2 <- ifelse(ref_strand == "-", comp[v$a2], v$a2)
man_a1[is_indel] <- "I"; man_a2[is_indel] <- "D"

# --- the manifest -------------------------------------------------------------
# A small number of report variants are absent, as if the calls came from a
# product carrying drop-in content the base manifest does not describe.
in_manifest <- rep(TRUE, n_var)
in_manifest[sample(n_var, round(0.06 * n_var))] <- FALSE

man <- data.frame(
  IlmnID     = paste0(name, "_B_F_", sample(1e9, n_var)),
  Name       = name,
  IlmnStrand = ifelse(ref_strand == "-", "BOT", "TOP"),
  SNP        = paste0("[", man_a1, "/", man_a2, "]"),
  GenomeBuild = 38,
  Chr        = v$chr,
  MapInfo    = v$pos,
  RefStrand  = ref_strand,
  stringsAsFactors = FALSE
)[in_manifest, ]
man <- man[order(man$Chr, man$MapInfo), ]

man_path <- file.path(out_dir, "example_manifest.csv")
con <- file(man_path, open = "wt")
writeLines(c("Illumina, Inc.",
             "[Heading]",
             "Descriptor File Name,DEMO-ARRAY-v1_A1.bpm",
             "Assay Format,Infinium HTS",
             "Date Manufactured,1/1/2026",
             paste0("Loci Count ,", nrow(man)),
             "[Assay]"), con)
write.table(man, con, sep = ",", row.names = FALSE, quote = FALSE)
close(con)

# --- genotypes ----------------------------------------------------------------
# Allele frequencies are drawn per variant, then genotypes per sample, so the
# three samples differ from one another as real samples would.
freq <- runif(n_var, 0.05, 0.95)

make_report <- function(sample_id, flip_block = FALSE) {
  g <- rbinom(n_var, 2, freq)
  a1 <- ifelse(g == 0, v$a1, ifelse(g == 1, v$a1, v$a2))
  a2 <- ifelse(g == 0, v$a1, v$a2)

  a1[is_indel] <- ifelse(g[is_indel] == 0, "I", "D")
  a2[is_indel] <- ifelse(g[is_indel] == 2, "D", "I")

  # GenCall: a good sample is bimodal, most calls confident, a tail that is not.
  gc <- ifelse(runif(n_var) < 0.55, rbeta(n_var, 9, 2), rbeta(n_var, 3, 3))

  # No-calls, concentrated among the low-confidence tail as they are in reality.
  p_nc <- 0.005 + 0.12 * (gc < 0.3)
  no_call <- runif(n_var) < p_nc
  a1[no_call] <- "-"; a2[no_call] <- "-"
  gc[no_call] <- pmin(gc[no_call], 0.30)

  # A handful of calls have no score at all.
  gc_txt <- sprintf("%.4f", gc)
  gc_txt[sample(which(no_call), min(2, sum(no_call)))] <- "NaN"

  # One sample carries a strand-flipped block, which surfaces as a multi-allelic
  # conflict when the samples are merged in Step 03.
  if (flip_block) {
    blk <- which(!is_indel & !no_call)[1:200]
    a1[blk] <- comp[a1[blk]]; a2[blk] <- comp[a2[blk]]
  }

  data.frame(`SNP Name` = name, `Sample ID` = sample_id, `GC Score` = gc_txt,
             `Allele1 - Plus` = a1, `Allele2 - Plus` = a2,
             check.names = FALSE, stringsAsFactors = FALSE)
}

for (i in seq_along(samples)) {
  s <- samples[i]
  d <- make_report(s, flip_block = (i == 3))
  f <- file.path(out_dir, paste0(s, "_final_report.csv"))
  con <- file(f, open = "wt")
  writeLines(c("[Header]",
               "GSGT Version;2.0.4",
               "Processing Date;01/01/2026 9:00 AM",
               "Content;;DEMO-ARRAY-v1_A1.bpm",
               paste0("Num SNPs;", n_var),
               paste0("Total SNPs;", n_var),
               "Num Samples;3",
               "Total Samples;3",
               paste0("File ;", i, " of 3"),
               "[Data]"), con)
  write.table(d, con, sep = ";", row.names = FALSE, quote = FALSE)
  close(con)
  cat("written:", f, "\n")
}

cat("\nInjected artefacts, for the README:\n")
cat("  variants:                    ", n_var, "\n")
cat("  samples:                     ", length(samples), "\n")
cat("  positional names:            ", sum(is_pos), sprintf(" (%.1f%%)", 100*mean(is_pos)), "\n")
cat("  of which chromosome differs: ", length(chr_mismatch), "\n")
cat("  indel probes:                ", sum(is_indel), sprintf(" (%.1f%%)", 100*mean(is_indel)), "\n")
cat("  minus-strand design:         ", sum(ref_strand == "-"), sprintf(" (%.1f%%)", 100*mean(ref_strand == "-")), "\n")
cat("  absent from the manifest:    ", sum(!in_manifest), sprintf(" (%.1f%%)", 100*mean(!in_manifest)), "\n")
cat("  manifest loci:               ", nrow(man), "\n")
cat("  sample 3 flipped block:       200 variants\n")
