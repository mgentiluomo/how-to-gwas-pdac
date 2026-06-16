#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop(
    paste(
      "Usage: Rscript 01_initial_qc_plots.R",
      "<bim_file> <afreq_file> <smiss_file> <vmiss_file> <het_file> <output_prefix>"
    ),
    call. = FALSE
  )
}

bim_file <- args[1]
afreq_file <- args[2]
smiss_file <- args[3]
vmiss_file <- args[4]
het_file <- args[5]
output_prefix <- args[6]

ensure_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

clean_names <- function(x) {
  names(x) <- sub("^#", "", names(x))
  x
}

read_plink_table <- function(path) {
  ensure_file(path)
  clean_names(read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  ))
}

first_matching_col <- function(data, candidates, label) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop(
      "Could not find ", label, " column. Available columns: ",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
  hit[1]
}

missing_rate <- function(data, label) {
  if ("F_MISS" %in% names(data)) {
    return(as.numeric(data[["F_MISS"]]))
  }

  if (all(c("MISSING_CT", "OBS_CT") %in% names(data))) {
    missing_ct <- as.numeric(data[["MISSING_CT"]])
    obs_ct <- as.numeric(data[["OBS_CT"]])
    return(missing_ct / obs_ct)
  }

  stop(
    "Could not calculate ", label,
    " missingness. Expected F_MISS or MISSING_CT/OBS_CT columns.",
    call. = FALSE
  )
}

finite_values <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

make_png <- function(path) {
  png(filename = path, width = 1200, height = 800, res = 150)
}

plot_empty <- function(title, message) {
  plot.new()
  title(main = title)
  text(0.5, 0.5, message)
}

chromosome_order <- function(chr) {
  chr_clean <- toupper(gsub("^CHR", "", as.character(chr)))
  chr_clean[chr_clean == "X"] <- "23"
  chr_clean[chr_clean == "Y"] <- "24"
  chr_clean[chr_clean %in% c("M", "MT")] <- "25"
  suppressWarnings(as.numeric(chr_clean))
}

message("Reading Step 01 QC outputs...")

ensure_file(bim_file)
bim <- read.table(
  bim_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("CHR", "ID", "CM", "POS", "A1", "A2")
)
bim$POS <- as.numeric(bim$POS)

afreq <- read_plink_table(afreq_file)
smiss <- read_plink_table(smiss_file)
vmiss <- read_plink_table(vmiss_file)
het <- read_plink_table(het_file)

# 1. Variant density per chromosome -----------------------------------------
variant_density_png <- paste0(output_prefix, "_variant_density_by_chromosome.png")
variant_density_tsv <- paste0(output_prefix, "_variant_density_by_chromosome.tsv")

chr_split <- split(bim, bim$CHR)
variant_density <- do.call(rbind, lapply(names(chr_split), function(chr) {
  x <- chr_split[[chr]]
  max_pos_mb <- max(x$POS, na.rm = TRUE) / 1e6
  max_pos_mb <- max(max_pos_mb, 1e-6)
  data.frame(
    CHR = chr,
    N_VARIANTS = nrow(x),
    OBSERVED_SPAN_MB = max_pos_mb,
    VARIANTS_PER_MB = nrow(x) / max_pos_mb,
    stringsAsFactors = FALSE
  )
}))

variant_density$CHR_ORDER <- chromosome_order(variant_density$CHR)
variant_density <- variant_density[order(
  is.na(variant_density$CHR_ORDER),
  variant_density$CHR_ORDER,
  variant_density$CHR
), ]

write.table(
  variant_density[, c("CHR", "N_VARIANTS", "OBSERVED_SPAN_MB", "VARIANTS_PER_MB")],
  file = variant_density_tsv,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

make_png(variant_density_png)
barplot(
  variant_density$VARIANTS_PER_MB,
  names.arg = variant_density$CHR,
  las = 2,
  col = "#4C78A8",
  border = NA,
  main = "Variant Density by Chromosome",
  xlab = "Chromosome",
  ylab = "Variants per observed Mb"
)
dev.off()

# 2. Sample missingness -------------------------------------------------------
sample_missing_png <- paste0(output_prefix, "_sample_missingness_histogram.png")
sample_miss <- finite_values(missing_rate(smiss, "sample-level"))

make_png(sample_missing_png)
if (length(sample_miss) == 0) {
  plot_empty("Sample Missingness", "No finite sample missingness values found")
} else {
  hist(
    sample_miss,
    breaks = 40,
    col = "#59A14F",
    border = "white",
    main = "Missingness Rate per Individual",
    xlab = "Sample missingness rate",
    ylab = "Number of samples"
  )
  abline(v = 0.02, col = "#E15759", lwd = 2, lty = 2)
  legend("topright", legend = "2% sample filter threshold", col = "#E15759", lty = 2, lwd = 2, bty = "n")
}
dev.off()

# 3. Variant missingness ------------------------------------------------------
variant_missing_png <- paste0(output_prefix, "_variant_missingness_histogram.png")
variant_miss <- finite_values(missing_rate(vmiss, "variant-level"))

make_png(variant_missing_png)
if (length(variant_miss) == 0) {
  plot_empty("Variant Missingness", "No finite variant missingness values found")
} else {
  hist(
    variant_miss,
    breaks = 50,
    col = "#F28E2B",
    border = "white",
    main = "Missingness Rate per Variant",
    xlab = "Variant missingness rate",
    ylab = "Number of variants"
  )
  abline(v = 0.05, col = "#E15759", lwd = 2, lty = 2)
  legend("topright", legend = "5% variant filter threshold", col = "#E15759", lty = 2, lwd = 2, bty = "n")
}
dev.off()

# 4. Allele frequency distribution ------------------------------------------
allele_freq_png <- paste0(output_prefix, "_allele_frequency_distribution.png")
freq_col <- first_matching_col(
  afreq,
  c("ALT_FREQS", "ALT_FREQ", "A1_FREQ", "MAF"),
  "allele frequency"
)
allele_freq <- as.numeric(sub(",.*", "", as.character(afreq[[freq_col]])))
allele_freq <- finite_values(allele_freq)
minor_allele_freq <- pmin(allele_freq, 1 - allele_freq)

make_png(allele_freq_png)
if (length(minor_allele_freq) == 0) {
  plot_empty("Allele Frequency Distribution", "No finite allele frequency values found")
} else {
  hist(
    minor_allele_freq,
    breaks = 50,
    col = "#B07AA1",
    border = "white",
    main = "Minor Allele Frequency Distribution",
    xlab = "Minor allele frequency",
    ylab = "Number of variants"
  )
  abline(v = 0.01, col = "#E15759", lwd = 2, lty = 2)
  legend("topright", legend = "1% MAF reference line", col = "#E15759", lty = 2, lwd = 2, bty = "n")
}
dev.off()

# 5. Heterozygote rate distribution -----------------------------------------
heterozygosity_png <- paste0(output_prefix, "_heterozygote_rate_distribution.png")
obs_hom_col <- first_matching_col(het, c("O(HOM)", "O_HOM"), "observed homozygote")
nonmissing_col <- first_matching_col(het, c("N(NM)", "N_NM"), "non-missing genotype count")

obs_hom <- as.numeric(het[[obs_hom_col]])
nonmissing <- as.numeric(het[[nonmissing_col]])
het_rate <- (nonmissing - obs_hom) / nonmissing
het_rate <- finite_values(het_rate)

make_png(heterozygosity_png)
if (length(het_rate) == 0) {
  plot_empty("Heterozygote Rate Distribution", "No finite heterozygote rate values found")
} else {
  hist(
    het_rate,
    breaks = 40,
    col = "#76B7B2",
    border = "white",
    main = "Heterozygote Rate Distribution",
    xlab = "Observed heterozygote rate per individual",
    ylab = "Number of samples"
  )
  abline(v = median(het_rate, na.rm = TRUE), col = "#4E79A7", lwd = 2)
  legend("topright", legend = "Median", col = "#4E79A7", lty = 1, lwd = 2, bty = "n")
}
dev.off()

message("Plots written:")
message("  ", variant_density_png)
message("  ", sample_missing_png)
message("  ", variant_missing_png)
message("  ", allele_freq_png)
message("  ", heterozygosity_png)
message("Summary table written:")
message("  ", variant_density_tsv)
