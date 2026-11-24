## Shared settings for the manuscript figures (scripts/dev/figures/).
## Every plotted number is read from a results file of the canonical run or of
## Section 8; each figure script names its files in its header and in the
## manifest it writes next to the figure.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
fig_dir <- if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
ROOT <- normalizePath(file.path(fig_dir, "..", "..", ".."))
setwd(ROOT)
OUTDIR <- Sys.getenv("FIG_OUT", file.path(dirname(ROOT), "manuscript", "figures"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

WIDTH_MM <- 174
FONT <- "Liberation Sans"
PT <- 7                      # base point size
INK <- "grey15"; MID <- "grey45"; LIGHT <- "grey80"; ACCENT <- "#1F4E79"; ACCENT2 <- "#B5651D"

S <- function(sec, f) file.path("sections", sec, "results", f)
need <- function(f) { if (!file.exists(f)) stop("missing results file: ", f); f }
rd <- function(f) read.delim(need(f), stringsAsFactors = FALSE, check.names = FALSE)

## draw(fun) writes <name>.pdf (cairo, fonts embedded) and <name>.png (600 dpi)
save_fig <- function(name, height_mm, fun, sources) {
  w <- WIDTH_MM / 25.4; h <- height_mm / 25.4
  cairo_pdf(file.path(OUTDIR, paste0(name, ".pdf")), width = w, height = h, family = FONT, pointsize = PT)
  fun(); invisible(dev.off())
  png(file.path(OUTDIR, paste0(name, ".png")), width = w, height = h, units = "in", res = 600,
      type = "cairo", family = FONT, pointsize = PT)
  fun(); invisible(dev.off())
  src <- data.frame(figure = name, file = sources,
                    md5 = unname(tools::md5sum(sources)), stringsAsFactors = FALSE)
  write.table(src, file.path(OUTDIR, paste0(name, "_sources.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("wrote", file.path(OUTDIR, paste0(name, c(".pdf", ".png", "_sources.tsv"))), sep = "\n")
}

par_base <- function(...) par(las = 1, mgp = c(2.2, 0.5, 0), tcl = -0.25, col.axis = INK, col.lab = INK,
                              fg = INK, bty = "l", cex.axis = 1, cex.lab = 1, ...)
## panel letter at the top-left corner of the figure region of the current panel
panel_label <- function(lab) {
  x <- grconvertX(0.01, "nfc", "user"); y <- grconvertY(0.985, "nfc", "user")
  text(x, y, lab, font = 2, cex = 1.3, adj = c(0, 1), xpd = NA)
}
