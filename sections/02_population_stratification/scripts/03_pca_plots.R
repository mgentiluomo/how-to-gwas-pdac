################################################################################
# Section 2: Population stratification — Step 03: visualising the PCA
#
# PURPOSE:
#   Three plots, each answering a different question.
#
#     Scree plot        how many components carry real signal?
#     PC1 versus PC2    does the cohort separate into ancestry groups?
#     Pairwise panel    is there structure beyond the first two components?
#
#   The third is the one most often skipped and most often needed. In a cohort
#   that looks homogeneous on PC1 and PC2, genuine structure frequently appears
#   on higher-order components, and it confounds association testing just the
#   same.
#
# INPUT:
#   - results/pca/pdac_demo_02_pca_all.eigenvec / .eigenval
#   - demo_data/sample_ancestry.tsv   (labels, for colouring only)
#
# OUTPUT:
#   - results/pca/pdac_demo_02_pca_scree.png
#   - results/pca/pdac_demo_02_pca_pc1_pc2.png
#   - results/pca/pdac_demo_02_pca_pairwise.png
#   - results/pca/pdac_demo_02_pca_variance.tsv
#
# NOTE ON THE LABELS:
#   The demonstration dataset ships with true ancestry labels, which is a
#   teaching convenience. In a real study no such labels exist: the groups are
#   what the PCA reveals, or what projection onto a reference panel assigns.
#   Here the labels serve to confirm that the genetic data and the recorded
#   metadata agree, which is itself a quality check worth performing.
#
# Base R only, so that nothing has to be installed.
################################################################################

out_dir   <- "results/pca"
eigenvec  <- file.path(out_dir, "pdac_demo_02_pca_all.eigenvec")
eigenval  <- file.path(out_dir, "pdac_demo_02_pca_all.eigenval")
anc_file  <- "demo_data/sample_ancestry.tsv"

pcs <- read.table(eigenvec, header = TRUE, comment.char = "", check.names = FALSE)
names(pcs)[1] <- sub("^#", "", names(pcs)[1])
val <- scan(eigenval, quiet = TRUE)

anc <- read.table(anc_file, header = FALSE,
                  col.names = c("IID", "group"), stringsAsFactors = FALSE)
pcs <- merge(pcs, anc, by = "IID", all.x = TRUE)
pcs$group <- toupper(pcs$group)

groups <- sort(unique(pcs$group))
# Colour-blind safe, and distinguishable in greyscale print.
palette_map <- c(EUR = "#0072B2", AFR = "#D55E00", EAS = "#009E73")
cols <- palette_map[pcs$group]
cols[is.na(cols)] <- "grey50"

cat("Individuals plotted:", nrow(pcs), "\n")
print(table(pcs$group))

# --- 1. Scree plot -----------------------------------------------------------
pct <- 100 * val / sum(val)
png(file.path(out_dir, "pdac_demo_02_pca_scree.png"),
    width = 1600, height = 1200, res = 200)
plot(seq_along(pct), pct, type = "b", pch = 19, col = "#0072B2",
     xlab = "Principal component",
     ylab = "Variance explained (%)",
     main = "Scree plot: where the signal stops")
abline(h = 0, col = "grey80")
dev.off()

write.table(
  data.frame(PC = seq_along(val), eigenvalue = val,
             variance_explained_pct = round(pct, 3)),
  file.path(out_dir, "pdac_demo_02_pca_variance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nVariance explained by the first four components:\n")
cat(paste0("  PC", 1:4, ": ", sprintf("%.2f%%", pct[1:4]), collapse = "\n"), "\n")

# --- 2. PC1 versus PC2 -------------------------------------------------------
png(file.path(out_dir, "pdac_demo_02_pca_pc1_pc2.png"),
    width = 1600, height = 1400, res = 200)
plot(pcs$PC1, pcs$PC2, col = cols, pch = 19, cex = 0.7,
     xlab = sprintf("PC1 (%.1f%%)", pct[1]),
     ylab = sprintf("PC2 (%.1f%%)", pct[2]),
     main = "Principal components 1 and 2")
legend("topright", legend = groups, col = palette_map[groups],
       pch = 19, bty = "n")
dev.off()

# --- 3. Pairwise panel, PC1 to PC6 -------------------------------------------
png(file.path(out_dir, "pdac_demo_02_pca_pairwise.png"),
    width = 1800, height = 1300, res = 190)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
pairs_to_plot <- list(c(1, 2), c(3, 4), c(5, 6),
                      c(1, 3), c(2, 4), c(4, 6))
for (p in pairs_to_plot) {
  plot(pcs[[paste0("PC", p[1])]], pcs[[paste0("PC", p[2])]],
       col = cols, pch = 19, cex = 0.5,
       xlab = sprintf("PC%d (%.1f%%)", p[1], pct[p[1]]),
       ylab = sprintf("PC%d (%.1f%%)", p[2], pct[p[2]]))
}
par(op)
dev.off()

cat("\nWritten:\n")
cat("  ", file.path(out_dir, "pdac_demo_02_pca_scree.png"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_02_pca_pc1_pc2.png"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_02_pca_pairwise.png"), "\n")
cat("  ", file.path(out_dir, "pdac_demo_02_pca_variance.tsv"), "\n")
cat("\n=== NEXT STEP ===\n\n")
cat("  bash scripts/02_population_stratification/04_define_analysis_set.sh\n\n")
