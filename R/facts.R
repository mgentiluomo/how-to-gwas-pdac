# Read the figures quoted in the pages from the committed results files.
#
# Every number a page states about the demo run should come from a file that
# scripts/dev/run_all.sh wrote, so that a rerun cannot leave the page behind.
# Pages load this file in a hidden setup chunk and quote values inline, e.g.
#   `r num(qc_counts("MAF filter", "variants"))`
# Base R only, so rendering adds no package dependency.

.facts_root <- local({
  root <- NULL
  function() {
    if (!is.null(root)) return(root)
    d <- Sys.getenv("QUARTO_PROJECT_DIR")
    if (!nzchar(d)) {
      d <- normalizePath(getwd())
      while (!file.exists(file.path(d, "_quarto.yml"))) {
        up <- dirname(d)
        if (up == d) stop("facts.R: project root (_quarto.yml) not found")
        d <- up
      }
    }
    root <<- d
    d
  }
})

# Path to a committed result, e.g. res("04A_association_binary", "x.tsv")
res <- function(section, file) {
  p <- file.path(.facts_root(), "sections", section, "results", file)
  if (!file.exists(p)) stop("facts.R: missing results file ", p)
  p
}

# Tab-separated table with a header
tab <- function(section, file)
  read.delim(res(section, file), stringsAsFactors = FALSE, check.names = FALSE)

# Two-column metric/value table as a named character vector
kv <- function(section, file) {
  d <- tab(section, file)
  setNames(as.character(d[[2]]), d[[1]])
}

# Value of one metric in a metric/value table; stops if absent
val <- function(section, file, key) {
  v <- kv(section, file)
  if (!key %in% names(v)) stop("facts.R: ", key, " not in ", file)
  v[[key]]
}

# Lines of a text file; count of non-empty lines
lines_of <- function(section, file) readLines(res(section, file), warn = FALSE)
n_lines  <- function(section, file) sum(nzchar(trimws(lines_of(section, file))))

# First capture group of a regular expression in a text report; stops if absent
grab <- function(section, file, pattern, n = 1L) {
  x <- lines_of(section, file)
  hit <- regmatches(x, regexec(pattern, x, perl = TRUE))
  hit <- Filter(length, hit)
  if (length(hit) < n) stop("facts.R: pattern not found in ", file, ": ", pattern)
  hit[[n]][2]
}

# ---- formatting --------------------------------------------------------------
num <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
dec <- function(x, d = 2) formatC(as.numeric(x), format = "f", digits = d, big.mark = ",")
pct <- function(x, d = 1) paste0(formatC(as.numeric(x), format = "f", digits = d), "%")
# P value in the style the pages use: 1.8e-9, not 1.75e-09
pv  <- function(x, d = 1) {
  s <- formatC(as.numeric(x), format = "e", digits = d)
  sub("e\\+?(-?)0*([0-9])", "e\\1\\2", s)
}
# P value as prose: 1.8 × 10⁻⁹
sci <- function(x, d = 1) {
  x <- as.numeric(x)
  e <- floor(log10(x)); m <- x / 10^e
  if (round(m, d) >= 10) { m <- m / 10; e <- e + 1 }
  sup <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴",
           "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹",
           "-" = "⁻")
  ex <- paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")
  paste0(formatC(m, format = "f", digits = d), " × 10", ex)
}

# ---- section helpers ---------------------------------------------------------
# Section 1B QC step table: qc_counts("MAF filter", "variants")
qc_counts <- function(step, col) {
  d <- tab("01B_genotyping_qc", "pdac_demo_09_qc_counts.tsv")
  r <- d[d$step == step, ]
  if (nrow(r) != 1L) stop("facts.R: QC step not found: ", step)
  r[[col]]
}
