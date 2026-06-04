#!/usr/bin/env bash
# =============================================================================
#  download_data.sh
#  Download the demonstration dataset for the PDAC GWAS tutorial from the
#  GitHub Release into demo_dataset/data/.
#
#  The genotype/phenotype files are large binaries and are NOT stored in git;
#  they live on the Release and are fetched here.
#
#  Usage:  bash demo_dataset/download_data.sh
# =============================================================================
set -euo pipefail

REPO="mgentiluomo/how-to-gwas-pdac"
TAG="v0.1-data"

# Resolve destination relative to this script, so it works from any CWD
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${SCRIPT_DIR}/data"

FILES=(
  "pdac_demo.bed"
  "pdac_demo.bim"
  "pdac_demo.fam"
  "phenotype.txt"
  "covariates.txt"
  "survival.txt"
  "sample_ancestry.tsv"
)

BASE_URL="https://github.com/${REPO}/releases/download/${TAG}"

mkdir -p "${DEST}"
echo ">>> Downloading demonstration dataset (release ${TAG}) into ${DEST}/"
echo ""

for f in "${FILES[@]}"; do
  out="${DEST}/${f}"
  if [ -f "${out}" ]; then
    echo "    [skip] ${f} already present"
    continue
  fi
  echo "    [get ] ${f}"
  if ! curl -fSL --retry 3 -o "${out}" "${BASE_URL}/${f}"; then
    echo ""
    echo "ERROR: could not download ${f}"
    echo "Check that release '${TAG}' exists at:"
    echo "  https://github.com/${REPO}/releases"
    exit 1
  fi
done

echo ""
echo ">>> Done. Files in ${DEST}/:"
ls -lh "${DEST}/"

echo ""
echo ">>> Quick sanity check (requires plink):"
if command -v plink >/dev/null 2>&1; then
  plink --bfile "${DEST}/pdac_demo" --freq --out /tmp/_dl_check 2>/dev/null \
    | grep -E "variants|people" || true
  rm -f /tmp/_dl_check.*
else
  echo "    (plink not found in PATH — skipping; see ../env/software_versions.md)"
fi

echo ""
echo "Data ready in demo_dataset/data/. Start the tutorial with sections/00_introduction/"
