#!/bin/bash
# =============================================================================
# 02b_apply_subset.sh
# -----------------------------------------------------------------------------
# Uses PLINK2 to extract the 6,000 selected individuals (keep_ids.txt) from each
# per-chromosome HAPNEST file, then merges the chromosomes into a single
# PLINK1 binary fileset (.bed/.bim/.fam) that the rest of the tutorial uses.
#
# We deliberately keep the genotype calls intact here — NO quality control is
# applied. QC (call rate, HWE, MAF, sex check, relatedness ...) is the subject
# of Section 1B and is run later as a hands-on teaching exercise.
# =============================================================================

set -euo pipefail

# ---- configuration -----------------------------------------------------------
RAW_DIR="${RAW_DIR:-./data/hapnest_raw}"
SUB_DIR="${SUB_DIR:-./data/subset}"
KEEP_FILE="${KEEP_FILE:-${SUB_DIR}/keep_ids.txt}"
CHROMS="${CHROMS:-1 5 9 13 16 17}"
FORMAT="${FORMAT:-pgen}"          # 'pgen' (PLINK2 input) or 'bed' (PLINK1 input)
THREADS="${THREADS:-4}"
OUT_PREFIX="${OUT_PREFIX:-${SUB_DIR}/hapnest_6k}"

mkdir -p "${SUB_DIR}/by_chr"

# pick the right input flag for PLINK2
if [[ "${FORMAT}" == "pgen" ]]; then
    INFLAG="--pfile"
elif [[ "${FORMAT}" == "bed" ]]; then
    INFLAG="--bfile"
else
    echo "ERROR: FORMAT must be 'pgen' or 'bed'." >&2
    exit 1
fi

command -v plink2 >/dev/null 2>&1 || {
    echo "ERROR: plink2 not found on PATH." >&2; exit 1; }

# ---- 1. extract the selected individuals, per chromosome --------------------
for chr in ${CHROMS}; do
    echo "=== Subsetting chr${chr} ==="
    plink2 \
        ${INFLAG} "${RAW_DIR}/chr${chr}" \
        --keep "${KEEP_FILE}" \
        --make-bed \
        --threads "${THREADS}" \
        --out "${SUB_DIR}/by_chr/chr${chr}_sub"
done

# ---- 2. merge all chromosomes into one fileset ------------------------------
# Build a merge list containing every chromosome EXCEPT the first (the first is
# passed directly to --bfile and the rest are merged onto it).
FIRST_CHR=$(echo ${CHROMS} | awk '{print $1}')
MERGE_LIST="${SUB_DIR}/by_chr/merge_list.txt"
: > "${MERGE_LIST}"
for chr in ${CHROMS}; do
    [[ "${chr}" == "${FIRST_CHR}" ]] && continue
    echo "${SUB_DIR}/by_chr/chr${chr}_sub" >> "${MERGE_LIST}"
done

echo ""
echo "=== Merging chromosomes ${CHROMS} ==="
plink2 \
    --bfile "${SUB_DIR}/by_chr/chr${FIRST_CHR}_sub" \
    --pmerge-list "${MERGE_LIST}" bfile \
    --make-bed \
    --threads "${THREADS}" \
    --out "${OUT_PREFIX}"

echo ""
echo "Merged dataset:"
ls -lh "${OUT_PREFIX}".{bed,bim,fam}
echo ""
echo "Variants : $(wc -l < "${OUT_PREFIX}.bim")"
echo "Samples  : $(wc -l < "${OUT_PREFIX}.fam")"
echo ""
echo "Next step:  Rscript scripts/03_simulate_phenotype.R"
