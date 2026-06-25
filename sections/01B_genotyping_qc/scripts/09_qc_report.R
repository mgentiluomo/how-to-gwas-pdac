#!/usr/bin/env Rscript

################################################################################
# Section 1B: Genotyping QC - Step 09 PDF report renderer
#
# This script is called by 09_qc_summary.sh. It reads the structured results from
# Steps 01-08 and renders one combined QC report PDF using base R only.
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(
    "Usage: 09_qc_report.R <dataset_name> <seed> <counts_file> <out_dir> ",
    "<summary_file> <report_pdf>",
    call. = FALSE
  )
}

dataset_name <- args[[1]]
seed <- as.integer(args[[2]])
counts_file <- args[[3]]
out_dir <- args[[4]]
summary_file <- args[[5]]
report_pdf <- args[[6]]

set.seed(seed)

theme <- list(
  navy = "#12355B",
  blue = "#1F77B4",
  green = "#2E7D32",
  orange = "#C77700",
  red = "#B23A48",
  gray = "#59636E",
  light_gray = "#EEF2F5",
  pale_blue = "#EAF3FB",
  ink = "#1F2933"
)

fmt_int <- function(x) {
  format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_pct <- function(numerator, denominator) {
  numerator <- as.numeric(numerator)
  denominator <- as.numeric(denominator)
  if (is.na(denominator) || denominator == 0) return("NA")
  sprintf("%.1f%%", 100 * numerator / denominator)
}

clean_names <- function(data) {
  names(data) <- sub("^#", "", names(data))
  data
}

read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- tryCatch(
    read.table(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  clean_names(out)
}

find_col <- function(data, candidates) {
  if (is.null(data)) return(NA_character_)
  nms <- names(data)
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

num_col <- function(data, candidates) {
  col <- find_col(data, candidates)
  if (is.na(col)) return(NULL)
  suppressWarnings(as.numeric(data[[col]]))
}

value_or_na <- function(data, metric) {
  if (is.null(data) || nrow(data) == 0) return("NA")
  if (!all(c("metric", "value") %in% names(data))) return("NA")
  hit <- data$value[data$metric == metric]
  if (length(hit) == 0) return("NA")
  as.character(hit[[1]])
}

start_page <- function(title, subtitle = NULL) {
  par(mfrow = c(1, 1), fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0), family = "sans", new = FALSE)
  plot.new()
  rect(0, 0, 1, 1, col = "#FBFCFE", border = NA)
  rect(0, 0.93, 1, 1, col = theme$navy, border = NA)
  text(0.04, 0.965, title, adj = c(0, 0.5), col = "white", font = 2, cex = 1.25)
  if (!is.null(subtitle)) {
    text(0.96, 0.965, subtitle, adj = c(1, 0.5), col = "#D7E7F5", cex = 0.72)
  }
}

draw_wrapped_text <- function(x, y, text_value, width = 85, cex = 0.72, line_height = 0.028,
                              col = theme$ink, font = 1) {
  lines <- unlist(strwrap(text_value, width = width))
  if (length(lines) == 0) return(y)
  for (line in lines) {
    text(x, y, line, adj = c(0, 1), cex = cex, col = col, font = font)
    y <- y - line_height
  }
  y
}

draw_table <- function(data, x = 0.05, y = 0.85, width = 0.90, row_height = 0.045,
                       max_rows = 14, cex = 0.58, title = NULL) {
  if (!is.null(title)) {
    text(x, y + 0.045, title, adj = c(0, 0.5), font = 2, cex = 0.82, col = theme$navy)
  }

  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    text(x, y, "Not available", adj = c(0, 1), cex = 0.72, col = theme$gray)
    return(invisible(y - row_height))
  }

  shown <- data
  truncated <- FALSE
  if (nrow(shown) > max_rows) {
    shown <- shown[seq_len(max_rows), , drop = FALSE]
    truncated <- TRUE
  }

  shown[] <- lapply(shown, function(col) {
    col <- as.character(col)
    col[is.na(col)] <- ""
    too_long <- nchar(col) > 26
    col[too_long] <- paste0(substr(col[too_long], 1, 23), "...")
    col
  })

  n_col <- ncol(shown)
  col_width <- width / n_col
  x_left <- x + (seq_len(n_col) - 1) * col_width
  x_text <- x_left + 0.008

  rect(x, y - row_height, x + width, y, col = theme$navy, border = "white")
  for (i in seq_len(n_col)) {
    text(x_text[i], y - row_height / 2, names(shown)[i],
         adj = c(0, 0.5), cex = cex, col = "white", font = 2)
  }

  current_y <- y - row_height
  for (r in seq_len(nrow(shown))) {
    fill <- if (r %% 2 == 1) "white" else theme$light_gray
    rect(x, current_y - row_height, x + width, current_y, col = fill, border = "white")
    for (c in seq_len(n_col)) {
      text(x_text[c], current_y - row_height / 2, shown[r, c],
           adj = c(0, 0.5), cex = cex, col = theme$ink)
    }
    current_y <- current_y - row_height
  }

  if (truncated) {
    text(x, current_y - 0.02, sprintf("Showing first %s rows.", max_rows),
         adj = c(0, 1), cex = 0.58, col = theme$gray)
  }

  invisible(current_y)
}

draw_mini_bar <- function(labels, values, x = 0.05, y = 0.10, width = 0.35, height = 0.20,
                          title = NULL, colors = NULL) {
  values <- suppressWarnings(as.numeric(values))
  keep <- is.finite(values)
  labels <- labels[keep]
  values <- values[keep]

  if (!is.null(title)) {
    text(x, y + height + 0.045, title, adj = c(0, 0.5), font = 2,
         cex = 0.82, col = theme$navy)
  }

  if (length(values) == 0) {
    rect(x, y, x + width, y + height, col = "white", border = "#CDD6DF")
    text(x + width / 2, y + height / 2, "Not available", cex = 0.72, col = theme$gray)
    return(invisible(FALSE))
  }

  if (is.null(colors)) {
    colors <- rep(theme$blue, length(values))
  }
  colors <- rep(colors, length.out = length(values))

  rect(x, y, x + width, y + height, col = "white", border = "#CDD6DF")
  max_value <- max(values, na.rm = TRUE)
  if (max_value <= 0) max_value <- 1

  n <- length(values)
  gap <- width * 0.08
  bar_width <- (width - gap * (n + 1)) / n

  for (i in seq_len(n)) {
    x0 <- x + gap + (i - 1) * (bar_width + gap)
    x1 <- x0 + bar_width
    bar_height <- (values[i] / max_value) * (height * 0.72)
    rect(x0, y + 0.045, x1, y + 0.045 + bar_height, col = colors[i], border = "white")
    text((x0 + x1) / 2, y + 0.035, labels[i], adj = c(0.5, 1), cex = 0.52, col = theme$ink)
    text((x0 + x1) / 2, y + 0.055 + bar_height, fmt_int(values[i]),
         adj = c(0.5, 0), cex = 0.58, col = theme$ink, font = 2)
  }

  invisible(TRUE)
}

plot_placeholder <- function(label) {
  plot.new()
  box(col = "#CDD6DF")
  text(0.5, 0.5, label, cex = 0.9, col = theme$gray)
}

plot_hist_if <- function(values, main, xlab, col = theme$blue, breaks = 40, transform = identity) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    plot_placeholder(paste(main, "\nnot available"))
    return(invisible(FALSE))
  }
  values <- transform(values)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    plot_placeholder(paste(main, "\nnot available"))
    return(invisible(FALSE))
  }
  hist(values, breaks = breaks, main = main, xlab = xlab, col = col, border = "white")
  grid(nx = NA, ny = NULL, col = "#D6DEE6")
  invisible(TRUE)
}

counts <- read_tsv(counts_file)
if (is.null(counts)) {
  stop("Could not read counts file: ", counts_file, call. = FALSE)
}

summary_lines <- if (file.exists(summary_file)) readLines(summary_file, warn = FALSE) else character()

smiss <- read_tsv(file.path(out_dir, paste0(dataset_name, "_01_qc.smiss")))
vmiss <- read_tsv(file.path(out_dir, paste0(dataset_name, "_01_qc.vmiss")))
afreq_initial <- read_tsv(file.path(out_dir, paste0(dataset_name, "_01_qc.afreq")))
het_initial <- read_tsv(file.path(out_dir, paste0(dataset_name, "_01_qc.het")))
het_step04 <- read_tsv(file.path(out_dir, paste0(dataset_name, "_04_het.het")))
hwe <- read_tsv(file.path(out_dir, paste0(dataset_name, "_06_hwe.hardy")))
afreq_final <- read_tsv(file.path(out_dir, paste0(dataset_name, "_08_afreq.afreq")))
density <- read_tsv(file.path(out_dir, paste0(dataset_name, "_01_qc_variant_density_by_chromosome.tsv")))
relatedness_summary <- read_tsv(file.path(out_dir, paste0(dataset_name, "_07_relatedness_summary.tsv")))
relatedness_diag <- read_tsv(file.path(out_dir, paste0(dataset_name, "_07_relatedness_pruning_diagnostics.tsv")))
king_pihat_summary <- read_tsv(file.path(out_dir, paste0(dataset_name, "_07_king_vs_pihat_summary.tsv")))
components <- read_tsv(file.path(out_dir, paste0(dataset_name, "_07_pruning_components.tsv")))
removal_counts <- read_tsv(file.path(out_dir, paste0(dataset_name, "_07_removal_phenotype_counts.tsv")))

final_row <- counts[nrow(counts), , drop = FALSE]
raw_row <- counts[1, , drop = FALSE]
sample_loss <- as.numeric(raw_row$samples) - as.numeric(final_row$samples)
variant_loss <- as.numeric(raw_row$variants) - as.numeric(final_row$variants)

pdf(report_pdf, width = 11, height = 8.5, onefile = TRUE)
on.exit(invisible(dev.off()), add = TRUE)

# Page 1: title and index
start_page("PDAC GWAS Genotyping QC Report", format(Sys.Date()))
text(0.05, 0.84, dataset_name, adj = c(0, 0.5), cex = 2.0, font = 2, col = theme$navy)
text(0.05, 0.79, "Section 1B: sample and variant-level genotyping quality control",
     adj = c(0, 0.5), cex = 0.92, col = theme$gray)

summary_cards <- data.frame(
  Metric = c("Raw samples", "Final samples", "Sample retention", "Raw variants",
             "Final variants", "Variant retention"),
  Value = c(
    fmt_int(raw_row$samples),
    fmt_int(final_row$samples),
    fmt_pct(final_row$samples, raw_row$samples),
    fmt_int(raw_row$variants),
    fmt_int(final_row$variants),
    fmt_pct(final_row$variants, raw_row$variants)
  ),
  stringsAsFactors = FALSE
)
draw_table(summary_cards, x = 0.05, y = 0.70, width = 0.46, row_height = 0.050,
           max_rows = 10, cex = 0.72, title = "Headline Counts")

index_text <- c(
  "1. Executive summary",
  "2. QC pipeline count table",
  "3. Sample and variant retention",
  "4. Initial QC distributions",
  "5. Heterozygosity, HWE, and MAF checks",
  "6. Relatedness pruning diagnostics",
  "7. Methods notes and interpretation"
)
text(0.62, 0.72, "Report Index", adj = c(0, 0.5), cex = 0.95, font = 2, col = theme$navy)
for (i in seq_along(index_text)) {
  text(0.62, 0.66 - (i - 1) * 0.055, index_text[i], adj = c(0, 0.5),
       cex = 0.78, col = theme$ink)
}

note <- paste(
  "This report is produced from the output files generated by Steps 01-08.",
  "It combines the count table, filtering rationale, and QC figures into one",
  "portable PDF that can be reviewed before moving to PCA and association testing."
)
invisible(draw_wrapped_text(0.05, 0.27, note, width = 120, cex = 0.78, line_height = 0.033))

# Page 2: count table
start_page("1. QC Pipeline Count Table", "Step 09")
qc_table <- data.frame(
  Step = sprintf("%02d. %s", counts$step_index, counts$step),
  Samples = fmt_int(counts$samples),
  Variants = fmt_int(counts$variants),
  `Samples lost` = fmt_int(counts$samples_lost),
  `Variants lost` = fmt_int(counts$variants_lost),
  check.names = FALSE
)
draw_table(qc_table, x = 0.04, y = 0.86, width = 0.92, row_height = 0.055,
           max_rows = 12, cex = 0.66, title = "Filtering trajectory")

text(0.05, 0.19, sprintf("Total sample loss: %s (%s)", fmt_int(sample_loss), fmt_pct(sample_loss, raw_row$samples)),
     adj = c(0, 0.5), cex = 0.84, col = theme$red, font = 2)
text(0.05, 0.14, sprintf("Total variant loss: %s (%s)", fmt_int(variant_loss), fmt_pct(variant_loss, raw_row$variants)),
     adj = c(0, 0.5), cex = 0.84, col = theme$orange, font = 2)
invisible(draw_wrapped_text(
  0.05, 0.09,
  "Interpretation: sample-level filters remove low-quality or non-independent individuals; variant-level filters remove unreliable markers before downstream population structure and association analyses.",
  width = 125,
  cex = 0.68,
  line_height = 0.025,
  col = theme$gray
))

# Page 3: retention plots
par(mfrow = c(1, 2), mar = c(8, 5, 4, 2), family = "sans")
sample_cols <- ifelse(counts$step == "Relatedness", theme$red, theme$blue)
barplot(
  counts$samples,
  names.arg = gsub(" ", "\n", counts$step),
  las = 2,
  col = sample_cols,
  border = "white",
  main = "Sample Retention Through QC",
  ylab = "Number of samples",
  cex.names = 0.70
)
grid(nx = NA, ny = NULL, col = "#D6DEE6")
variant_cols <- ifelse(counts$variants_lost > 0, theme$green, "#7FA9C9")
barplot(
  counts$variants / 1000,
  names.arg = gsub(" ", "\n", counts$step),
  las = 2,
  col = variant_cols,
  border = "white",
  main = "Variant Retention Through QC",
  ylab = "Number of variants (thousands)",
  cex.names = 0.70
)
grid(nx = NA, ny = NULL, col = "#D6DEE6")
mtext("2. Sample and Variant Retention", side = 3, outer = TRUE, line = -1.5,
      font = 2, cex = 1.2, col = theme$navy)

# Page 4: initial QC distributions
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), family = "sans")
if (!is.null(density) && all(c("CHR", "VARIANTS_PER_MB") %in% names(density))) {
  barplot(
    density$VARIANTS_PER_MB,
    names.arg = density$CHR,
    col = theme$blue,
    border = "white",
    main = "Variant Density by Chromosome",
    ylab = "Variants per Mb",
    xlab = "Chromosome",
    cex.names = 0.65
  )
  grid(nx = NA, ny = NULL, col = "#D6DEE6")
} else {
  plot_placeholder("Variant density\nnot available")
}
plot_hist_if(num_col(smiss, c("F_MISS")), "Sample Missingness", "Proportion missing genotypes",
             col = theme$blue, breaks = 30)
plot_hist_if(num_col(vmiss, c("F_MISS")), "Variant Missingness", "Proportion missing samples",
             col = theme$green, breaks = 50)
plot_hist_if(num_col(afreq_initial, c("ALT_FREQS", "ALT_FREQ", "MAF")), "Initial Allele Frequencies",
             "Alternate allele frequency", col = theme$orange, breaks = 50)
mtext("3. Initial QC Distributions", side = 3, outer = TRUE, line = -1.5,
      font = 2, cex = 1.2, col = theme$navy)

# Page 5: heterozygosity, HWE, and MAF checks
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), family = "sans")
if (!is.null(smiss) && !is.null(het_initial)) {
  fid_col_s <- find_col(smiss, c("FID"))
  iid_col_s <- find_col(smiss, c("IID", "ID"))
  fid_col_h <- find_col(het_initial, c("FID"))
  iid_col_h <- find_col(het_initial, c("IID", "ID"))
  miss_col <- find_col(smiss, c("F_MISS"))
  ohom_col <- find_col(het_initial, c("O(HOM)", "O.HOM."))
  obs_col <- find_col(het_initial, c("OBS_CT"))
  if (!any(is.na(c(fid_col_s, iid_col_s, fid_col_h, iid_col_h, miss_col, ohom_col, obs_col)))) {
    smiss$key <- paste(smiss[[fid_col_s]], smiss[[iid_col_s]], sep = "/")
    het_initial$key <- paste(het_initial[[fid_col_h]], het_initial[[iid_col_h]], sep = "/")
    merged <- merge(
      smiss[, c("key", miss_col)],
      het_initial[, c("key", ohom_col, obs_col)],
      by = "key"
    )
    names(merged) <- c("key", "missingness", "observed_hom", "observed_ct")
    het_rate <- 1 - suppressWarnings(as.numeric(merged$observed_hom)) /
      suppressWarnings(as.numeric(merged$observed_ct))
    missingness <- suppressWarnings(as.numeric(merged$missingness))
    keep <- is.finite(het_rate) & is.finite(missingness)
    if (sum(keep) > 0) {
      plot(
        missingness[keep],
        het_rate[keep],
        pch = 16,
        col = adjustcolor(theme$blue, alpha.f = 0.55),
        main = "Missingness vs Heterozygosity",
        xlab = "Proportion missing genotypes",
        ylab = "Heterozygosity rate"
      )
      grid(col = "#D6DEE6")
    } else {
      plot_placeholder("Missingness vs heterozygosity\nnot available")
    }
  } else {
    plot_placeholder("Missingness vs heterozygosity\nnot available")
  }
} else {
  plot_placeholder("Missingness vs heterozygosity\nnot available")
}
plot_hist_if(num_col(het_step04, c("F")), "Autosomal Inbreeding Coefficient",
             "F statistic", col = theme$blue, breaks = 35)
plot_hist_if(num_col(hwe, c("P")), "HWE P-value Distribution",
             "-log10(HWE p-value)", col = theme$green, breaks = 50,
             transform = function(x) -log10(pmax(x, .Machine$double.xmin)))
plot_hist_if(num_col(afreq_final, c("ALT_FREQS", "ALT_FREQ", "MAF")), "Final Allele Frequencies",
             "Alternate allele frequency after QC", col = theme$orange, breaks = 50)
mtext("4. Heterozygosity, HWE, and MAF Checks", side = 3, outer = TRUE, line = -1.5,
      font = 2, cex = 1.2, col = theme$navy)

# Page 6: relatedness diagnostics
start_page("5. Relatedness Pruning Diagnostics", "KING + PI_HAT")
diag_table <- NULL
if (!is.null(relatedness_diag) && all(c("metric", "value") %in% names(relatedness_diag))) {
  diag_table <- relatedness_diag
  names(diag_table) <- c("Metric", "Value")
}
draw_table(diag_table, x = 0.04, y = 0.84, width = 0.43, row_height = 0.041,
           max_rows = 13, cex = 0.57, title = "Headline diagnostics")

rel_table <- relatedness_summary
if (!is.null(rel_table) && all(c("relationship_category", "n_pairs", "n_unique_samples") %in% names(rel_table))) {
  names(rel_table) <- c("Relationship", "Pairs", "Unique samples")
  rel_table$Relationship <- gsub("_", " ", rel_table$Relationship)
}
draw_table(rel_table, x = 0.53, y = 0.84, width = 0.42, row_height = 0.041,
           max_rows = 7, cex = 0.55, title = "KING related-pair summary")

if (!is.null(removal_counts) && all(c("phenotype", "n_removed") %in% names(removal_counts))) {
  rc <- removal_counts[removal_counts$phenotype %in% c("control", "case", "unknown"), ]
  draw_mini_bar(
    labels = rc$phenotype,
    values = rc$n_removed,
    x = 0.04,
    y = 0.08,
    width = 0.43,
    height = 0.18,
    title = "Removed samples by phenotype",
    colors = c(theme$blue, theme$red, theme$gray)
  )
}

if (!is.null(king_pihat_summary) && all(c("metric", "value") %in% names(king_pihat_summary))) {
  kp <- king_pihat_summary
  names(kp) <- c("Metric", "Value")
  draw_table(kp, x = 0.53, y = 0.37, width = 0.42, row_height = 0.036,
             max_rows = 8, cex = 0.50, title = "KING vs PI_HAT agreement")
}

if (!is.null(components) && nrow(components) > 0) {
  comp <- components[seq_len(min(6, nrow(components))), , drop = FALSE]
  names(comp) <- gsub("_", " ", names(comp))
  draw_table(comp, x = 0.04, y = 0.38, width = 0.43, row_height = 0.035,
             max_rows = 6, cex = 0.43, title = "Largest relatedness components")
}

# Page 7: methods notes
start_page("6. Methods Notes and Interpretation", "Review before downstream analysis")
thresholds <- data.frame(
  Level = c("Sample", "Sample", "Sample", "Sample", "Variant", "Variant", "Variant"),
  Check = c("Call rate", "Sex check", "Heterozygosity", "Relatedness",
            "Call rate", "Hardy-Weinberg", "Minor allele frequency"),
  Threshold = c("--mind 0.02", "X chromosome F statistic", "F +/- 3 SD",
                "KING kinship > 0.1875", "--geno 0.05", "--hwe 1e-6 in controls",
                "--maf 0.01"),
  Rationale = c("Remove low-quality samples", "Detect swaps or metadata errors",
                "Detect contamination/outliers", "Preserve independence; keep cases where possible",
                "Remove poorly genotyped variants", "Remove likely genotyping errors",
                "Keep low-frequency signal in rare cancer"),
  stringsAsFactors = FALSE
)
draw_table(thresholds, x = 0.04, y = 0.84, width = 0.92, row_height = 0.050,
           max_rows = 9, cex = 0.50, title = "QC thresholds")

notes <- c(
  "Rare cancer setting: cases are precious, so phenotype-aware relatedness pruning is used.",
  "KING is used for final pruning because it is robust for relationship inference in GWAS QC; PI_HAT is reported for review and teaching.",
  "A large relatedness loss should trigger a review of recruitment design, family structure, and downstream relatedness-aware association options.",
  "The final dataset from Step 08 is the starting point for ancestry analysis, PCA, and association testing."
)
y <- 0.34
text(0.04, y, "Review notes", adj = c(0, 1), cex = 0.86, font = 2, col = theme$navy)
y <- y - 0.045
for (note in notes) {
  y <- draw_wrapped_text(0.06, y, paste0("- ", note), width = 130, cex = 0.68,
                         line_height = 0.030, col = theme$ink)
  y <- y - 0.012
}

# Page 8: plain-text appendix excerpt
start_page("7. Plain-Text Summary Appendix", basename(summary_file))
if (length(summary_lines) == 0) {
  text(0.05, 0.84, "Summary text file was not available.", adj = c(0, 1),
       cex = 0.75, col = theme$gray)
} else {
  appendix <- summary_lines[seq_len(min(length(summary_lines), 34))]
  y <- 0.86
  for (line in appendix) {
    text(0.04, y, line, adj = c(0, 1), cex = 0.49, family = "mono", col = theme$ink)
    y <- y - 0.023
  }
  if (length(summary_lines) > length(appendix)) {
    text(0.04, 0.06, "Appendix truncated in PDF; see the full text summary file for complete notes.",
         adj = c(0, 1), cex = 0.62, col = theme$gray)
  }
}

cat(sprintf("Combined QC report saved to: %s\n", report_pdf))
