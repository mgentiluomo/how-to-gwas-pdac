#!/bin/bash
# =============================================================================
# 01_download_hapnest.sh
# -----------------------------------------------------------------------------
# Downloads the per-chromosome HAPNEST pre-generated genotype files that the
# tutorial needs, plus the sample-description file that lists each synthetic
# individual's genetic-ancestry group.
#
# HAPNEST is hosted on EBI BioStudies under accession S-BSST936. The full
# dataset is enormous (1,008,000 individuals, all chromosomes); for teaching we
# only download SIX chromosomes — one per simulated causal locus — which keeps
# the download and all later steps light enough to run on a laptop.
#
# IMPORTANT
#   * File names and the exact download URL on EBI may change over time.
#     Always confirm the current paths on the S-BSST936 study page before a
#     fresh download, and adjust BASE_URL / file patterns below if needed.
#   * HAPNEST pre-generated data is distributed in PLINK2 format
#     (.pgen / .pvar / .psam). If your copy is in PLINK1 format
#     (.bed / .bim / .fam), set FORMAT=bed below; later scripts handle both.
# =============================================================================

set -euo pipefail

# ---- configuration (override by exporting the variable before running) ------
BASE_URL="${BASE_URL:-https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/936/S-BSST936/Files}"
RAW_DIR="${RAW_DIR:-./data/hapnest_raw}"
CHROMS="${CHROMS:-1 5 9 13 16 17}"        # one chromosome per causal locus
FORMAT="${FORMAT:-pgen}"                    # 'pgen' (PLINK2) or 'bed' (PLINK1)
SAMPLE_FILE="${SAMPLE_FILE:-hapnest_samples.tsv}"   # ancestry / metadata table

mkdir -p "${RAW_DIR}"

# ---- pick a download tool ----------------------------------------------------
if command -v wget >/dev/null 2>&1; then
    DL() { wget -c -q --show-progress -O "$2" "$1"; }
elif command -v curl >/dev/null 2>&1; then
    DL() { curl -fL --retry 3 -o "$2" "$1"; }
else
    echo "ERROR: neither wget nor curl is available." >&2
    exit 1
fi

echo "============================================================"
echo " HAPNEST download"
echo "   base URL : ${BASE_URL}"
echo "   format   : ${FORMAT}"
echo "   chroms   : ${CHROMS}"
echo "   target   : ${RAW_DIR}"
echo "============================================================"

# ---- choose the file extensions for the requested format --------------------
if [[ "${FORMAT}" == "pgen" ]]; then
    EXTS=(pgen pvar psam)
elif [[ "${FORMAT}" == "bed" ]]; then
    EXTS=(bed bim fam)
else
    echo "ERROR: FORMAT must be 'pgen' or 'bed' (got '${FORMAT}')." >&2
    exit 1
fi

# ---- download the sample / ancestry description file ------------------------
# This table maps each synthetic individual to a genetic-ancestry group; it is
# what 02_subset_individuals.R uses to build the 80/10/5/5 ancestry mix.
echo ""
echo ">>> Sample description file"
DL "${BASE_URL}/${SAMPLE_FILE}" "${RAW_DIR}/${SAMPLE_FILE}" \
    || echo "WARNING: could not fetch ${SAMPLE_FILE}; check its exact name on EBI."

# ---- download each chromosome ------------------------------------------------
for chr in ${CHROMS}; do
    echo ""
    echo ">>> Chromosome ${chr}"
    for ext in "${EXTS[@]}"; do
        url="${BASE_URL}/chr${chr}.${ext}"
        out="${RAW_DIR}/chr${chr}.${ext}"
        if [[ -s "${out}" ]]; then
            echo "    already present: ${out} (skipping)"
        else
            echo "    downloading chr${chr}.${ext} ..."
            DL "${url}" "${out}"
        fi
    done
done

echo ""
echo "============================================================"
echo " Download complete. Files in ${RAW_DIR}:"
ls -lh "${RAW_DIR}"
echo "============================================================"
echo ""
echo "Next step:  Rscript scripts/02_subset_individuals.R"
