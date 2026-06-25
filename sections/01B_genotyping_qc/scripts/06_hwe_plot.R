#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript 06_hwe_plot.R <hardy_file> <output_prefix>",
    call. = FALSE
  )
}

hardy_file <- args[1]
output_prefix <- args[2]
hwe_threshold <- 1e-6

ensure_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

clean_names <- function(x) {
  names(x) <- sub("^#", "", names(x))
  x
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

ensure_file(hardy_file)

message("Reading HWE results...")
hardy <- clean_names(read.table(
  hardy_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
))

if ("TEST" %in% names(hardy) && any(hardy$TEST == "ALL")) {
  hardy <- hardy[hardy$TEST == "ALL", ]
}

p_col <- first_matching_col(hardy, c("P", "P_HWE"), "HWE p-value")
p_values <- suppressWarnings(as.numeric(hardy[[p_col]]))
p_values <- p_values[is.finite(p_values) & p_values >= 0 & p_values <= 1]

hwe_plot_png <- paste0(output_prefix, "_pvalue_distribution.png")

png(filename = hwe_plot_png, width = 1200, height = 800, res = 150)

if (length(p_values) == 0) {
  plot.new()
  title(main = "Hardy-Weinberg P-value Distribution")
  text(0.5, 0.5, "No finite HWE p-values found")
} else {
  neg_log10_p <- -log10(pmax(p_values, .Machine$double.xmin))

  hist(
    neg_log10_p,
    breaks = 60,
    col = "#4C78A8",
    border = "white",
    main = "Hardy-Weinberg P-value Distribution",
    xlab = "-log10(HWE p-value)",
    ylab = "Number of variants"
  )
  abline(v = -log10(hwe_threshold), col = "#E15759", lwd = 2, lty = 2)
  legend(
    "topright",
    legend = "HWE filter threshold: p < 1e-6",
    col = "#E15759",
    lty = 2,
    lwd = 2,
    bty = "n"
  )
}

invisible(dev.off())

message("HWE p-value plot written:")
message("  ", hwe_plot_png)
