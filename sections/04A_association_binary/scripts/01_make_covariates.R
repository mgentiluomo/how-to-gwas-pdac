################################################################################
# Section 4A: Association testing — Step 01: assembling the covariates
#
# PURPOSE:
#   Combine the demographic covariates with the principal components from
#   Section 2 into the single file PLINK 2 expects.
#
#   Covariate adjustment is not a formality. It is the mechanism by which the
#   confounding identified in Section 2 is actually removed from the test. A
#   principal component that is computed and then not used has corrected
#   nothing.
#
# INPUT:
#   - demo_data/covariates.txt                    FID IID SEX AGE
#   - results/pca/pdac_demo_02_pca_eur.eigenvec   within-EUR components
#
# OUTPUT:
#   - results/assoc/pdac_demo_04A_covar.txt
#   - results/assoc/pdac_demo_04A_pc_screen.tsv   each PC tested against status
#
# ON THE NUMBER OF COMPONENTS:
#   Ten is a convention, not a result. This script also tests each component for
#   association with case status in a model containing sex and age, which is one
#   of the three checks the guide recommends. The others are the scree plot
#   (Section 2) and confirming that the genomic inflation factor and the leading
#   associations are stable when the number of components is varied.
#
#   A component that is strongly associated with case status is capturing
#   something that differs between cases and controls. In a well-designed study
#   that is ancestry; in a multi-centre study it may be recruitment, and no
#   number of components will fix a design in which cases and controls came from
#   different populations.
################################################################################

covar_in  <- "demo_data/covariates.txt"
eigenvec  <- "results/pca/pdac_demo_02_pca_eur.eigenvec"
pheno_in  <- "demo_data/phenotype.txt"
out_dir   <- "results/assoc"
covar_out <- file.path(out_dir, "pdac_demo_04A_covar.txt")
screen_out<- file.path(out_dir, "pdac_demo_04A_pc_screen.tsv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_tab <- function(f) {
  d <- read.table(f, header = TRUE, comment.char = "", check.names = FALSE)
  names(d)[1] <- sub("^#", "", names(d)[1])
  d
}

cov <- read_tab(covar_in)
pcs <- read_tab(eigenvec)
phe <- read_tab(pheno_in)

key <- intersect(c("FID", "IID"), intersect(names(cov), names(pcs)))
m <- merge(cov, pcs, by = key)

cat("covariate rows: ", nrow(cov), "\n")
cat("eigenvec rows:  ", nrow(pcs), "\n")
cat("merged rows:    ", nrow(m), "\n")
if (nrow(m) != nrow(pcs)) {
  warning("Merged rows differ from eigenvec rows: check that identifiers match.")
}

# --- screen each component against case status -------------------------------
pheno_col <- names(phe)[ncol(phe)]
d <- merge(m, phe[, c("IID", pheno_col)], by = "IID")
# PLINK coding is 1 = control, 2 = case; recode to 0/1 for glm.
d$status <- ifelse(d[[pheno_col]] == 2, 1, ifelse(d[[pheno_col]] == 1, 0, NA))
d <- d[!is.na(d$status), ]

pc_names <- grep("^PC[0-9]+$", names(d), value = TRUE)
res <- do.call(rbind, lapply(pc_names, function(pc) {
  f <- as.formula(paste("status ~ SEX + AGE +", pc))
  s <- summary(glm(f, data = d, family = binomial()))$coefficients
  data.frame(PC = pc,
             beta = unname(s[pc, 1]),
             se   = unname(s[pc, 2]),
             P    = unname(s[pc, 4]))
}))
res$nominally_associated <- ifelse(res$P < 0.05, "yes", "no")

write.table(res, screen_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nEach component tested against case status, adjusted for sex and age:\n")
print(transform(res, beta = round(beta, 4), se = round(se, 4),
                P = signif(P, 3)), row.names = FALSE)

n_assoc <- sum(res$P < 0.05)
cat("\nComponents nominally associated with status (P < 0.05):", n_assoc, "\n")
cat("Interpretation: in a single-ancestry analysis set, few or no components\n")
cat("should be associated with status. Several strongly associated components\n")
cat("would indicate that case and control ascertainment differ systematically,\n")
cat("which is a design problem rather than a covariate problem.\n")

# --- write the covariate file ------------------------------------------------
# PLINK 2 expects the first column to be named #FID (or #IID).
names(m)[which(names(m) == key[1])] <- paste0("#", key[1])
write.table(m, covar_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nWritten:\n")
cat("  ", covar_out, "\n")
cat("  ", screen_out, "\n")
cat("\n=== NEXT STEP ===\n\n")
cat("  bash scripts/04A_association_binary/02_association.sh\n\n")
