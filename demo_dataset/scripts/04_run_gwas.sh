#!/bin/bash
# =============================================================================
# 04_run_gwas.sh
# -----------------------------------------------------------------------------
# Reference end-to-end GWAS on the simulated dataset. This produces the
# "known-good" summary statistics and the Manhattan/QQ inputs that the
# manuscript reproduces. The detailed teaching of each step (QC in Section 1B,
# PCA in Section 2, association testing in Section 4) is done elsewhere — here
# we run a clean, minimal version so the example is fully reproducible.
#
# Steps:
#   1. Light QC          (MAF, genotype/sample call rate, HWE in controls)
#   2. LD pruning        (independent SNPs for principal-component analysis)
#   3. PCA               (10 principal components to control for ancestry)
#   4. Association test  (logistic regression, additive model, PCs as covars)
# =============================================================================

set -euo pipefail

# ---- configuration -----------------------------------------------------------
BFILE="${BFILE:-./data/subset/hapnest_6k}"
PHENO="${PHENO:-./data/subset/pheno.txt}"
OUT_DIR="${OUT_DIR:-./results}"
THREADS="${THREADS:-4}"
N_PC="${N_PC:-10}"

mkdir -p "${OUT_DIR}"
command -v plink2 >/dev/null 2>&1 || {
    echo "ERROR: plink2 not found on PATH." >&2; exit 1; }

# =============================================================================
# 1. Light quality control
#    --maf 0.01      drop very rare variants (unstable association estimates)
#    --geno 0.02     drop variants missing in >2% of samples
#    --mind 0.02     drop samples missing >2% of genotypes
#    --hwe 1e-6 ...   Hardy-Weinberg filter, evaluated in CONTROLS only
# =============================================================================
echo "=== Step 1/4: QC ==="
plink2 \
    --bfile "${BFILE}" \
    --pheno "${PHENO}" --pheno-name PHENO \
    --maf 0.01 --geno 0.02 --mind 0.02 \
    --hwe 1e-6 keep-fewhet \
    --make-bed \
    --threads "${THREADS}" \
    --out "${OUT_DIR}/qc"

# =============================================================================
# 2. LD pruning  — keep a set of approximately independent SNPs for PCA.
#    Window 200 kb, step 50 variants, r^2 threshold 0.2.
# =============================================================================
echo ""
echo "=== Step 2/4: LD pruning for PCA ==="
plink2 \
    --bfile "${OUT_DIR}/qc" \
    --indep-pairwise 200 50 0.2 \
    --threads "${THREADS}" \
    --out "${OUT_DIR}/prune"

# =============================================================================
# 3. Principal-component analysis on the pruned SNPs.
#    The first PCs capture genetic ancestry; we feed them into the GWAS as
#    covariates so that the EUR/AFR/EAS/SAS mix does not create false positives.
# =============================================================================
echo ""
echo "=== Step 3/4: PCA (${N_PC} PCs) ==="
plink2 \
    --bfile "${OUT_DIR}/qc" \
    --extract "${OUT_DIR}/prune.prune.in" \
    --pca "${N_PC}" \
    --threads "${THREADS}" \
    --out "${OUT_DIR}/pca"

# =============================================================================
# 4. Association testing: logistic regression, additive coding, PCs as covars.
#    --glm flags:
#      hide-covar   only report the SNP rows, not one row per PC
#      firth-fallback  use Firth correction when standard ML fails (rare alleles)
#      cols=...     write a tidy, complete summary-statistics table
#    --covar-variance-standardize stabilises the PC covariates.
# =============================================================================
echo ""
echo "=== Step 4/4: GWAS (logistic regression) ==="
plink2 \
    --bfile "${OUT_DIR}/qc" \
    --pheno "${PHENO}" --pheno-name PHENO \
    --covar "${OUT_DIR}/pca.eigenvec" \
    --covar-variance-standardize \
    --glm hide-covar firth-fallback \
        cols=chrom,pos,ref,alt,a1,ax,nobs,beta,orbeta,se,tz,p \
    --threads "${THREADS}" \
    --out "${OUT_DIR}/gwas"

# PLINK2 names the logistic output "<out>.PHENO.glm.logistic[.hybrid]".
RESULT=$(ls "${OUT_DIR}"/gwas.PHENO.glm.logistic* 2>/dev/null | head -n1 || true)
echo ""
echo "============================================================"
if [[ -n "${RESULT}" ]]; then
    echo " GWAS summary statistics: ${RESULT}"
    echo " Lines (incl. header):    $(wc -l < "${RESULT}")"
else
    echo " WARNING: could not locate the .glm.logistic output — check ${OUT_DIR}/gwas.log"
fi
echo "============================================================"
echo ""
echo "Next step:  Rscript scripts/05_plots.R"
