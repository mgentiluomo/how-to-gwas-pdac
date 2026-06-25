#!/usr/bin/env Rscript

# Configuration
args <- commandArgs(trailingOnly = TRUE)
het_file <- args[1]
output_prefix <- args[2]
seed <- as.numeric(args[3])

set.seed(seed)

# Read heterozygosity data
het <- read.table(het_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
names(het) <- sub("^#", "", names(het))

if (!"F" %in% names(het)) {
  stop("The PLINK2 .het file does not contain an F column.")
}
het$F <- as.numeric(het$F)

# Identify outliers: mean(F) ± 3 * SD(F)
mean_f <- mean(het$F, na.rm = TRUE)
sd_f <- sd(het$F, na.rm = TRUE)
lower_bound <- mean_f - 3 * sd_f
upper_bound <- mean_f + 3 * sd_f

# Flag outliers
het$outlier <- het$F < lower_bound | het$F > upper_bound

# Summary statistics
n_total <- nrow(het)
n_outliers <- sum(het$outlier, na.rm = TRUE)
percent_outliers <- 100 * n_outliers / n_total

cat(sprintf("Heterozygosity Summary:\n"))
cat(sprintf("  Mean F: %.4f\n", mean_f))
cat(sprintf("  SD F:   %.4f\n", sd_f))
cat(sprintf("  Lower outlier bound (F < %.4f): %.4f\n", lower_bound, lower_bound))
cat(sprintf("  Upper outlier bound (F > %.4f): %.4f\n", upper_bound, upper_bound))
cat(sprintf("  Total samples: %d\n", n_total))
cat(sprintf("  Outliers (3-SD): %d (%.2f%%)\n\n", n_outliers, percent_outliers))

# Write outlier list (for removal)
outliers <- het[het$outlier, c("FID", "IID")]
write.table(
  outliers,
  file = paste0(output_prefix, "_outliers.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = " "
)
cat(sprintf("Outlier list written to: %s_outliers.txt\n\n", output_prefix))

# Create visualization: histogram + density + outliers highlighted
pdf(file = paste0(output_prefix, "_outliers.pdf"), width = 10, height = 6)

hist(het$F, breaks = 50, main = "Distribution of Heterozygosity (F statistic)",
     xlab = "Inbreeding Coefficient (F)", ylab = "Frequency", col = "lightblue", border = "white")
abline(v = mean_f, col = "blue", lwd = 2, lty = 1, label = "Mean")
abline(v = lower_bound, col = "red", lwd = 2, lty = 2, label = "Outlier threshold")
abline(v = upper_bound, col = "red", lwd = 2, lty = 2)
legend("topright", c("Mean", "Outlier bounds (±3 SD)"), 
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# Q-Q plot
qqnorm(het$F, main = "Q-Q Plot: Heterozygosity")
qqline(het$F, col = "red")

dev.off()
cat(sprintf("Plot saved to: %s_outliers.pdf\n", output_prefix))
