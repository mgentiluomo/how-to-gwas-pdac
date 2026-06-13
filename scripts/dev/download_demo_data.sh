#!/usr/bin/env bash

################################################################################
# Download PDAC GWAS Demo Dataset
#
# Downloads demo dataset files to demo_data/ in your tutorial project folder.
# Works when called as: bash scripts/dev/download_demo_data.sh
#
# Usage: bash download_demo_data.sh
#
# Files downloaded:
#   - pdac_demo.bed, pdac_demo.bim, pdac_demo.fam (genotypes)
#   - phenotype.txt, covariates.txt, survival.txt (phenotypes)
#   - sample_ancestry.tsv (ancestry/ethnicity)
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

REPO='mgentiluomo/how-to-gwas-pdac'
TAG='v0.1-data'
BASE_URL="https://github.com/$REPO/releases/download/$TAG"

DEST_DIR='demo_data'

FILES=(
  'pdac_demo.bed'
  'pdac_demo.bim'
  'pdac_demo.fam'
  'phenotype.txt'
  'covariates.txt'
  'survival.txt'
  'sample_ancestry.tsv'
)

# Create destination directory
mkdir -p "$DEST_DIR"
echo ">>> Downloading dataset into $DEST_DIR/"
echo ""

# Download each file
for file in "${FILES[@]}"; do
    out_path="$DEST_DIR/$file"
    
    # Skip if already exists
    if [ -f "$out_path" ]; then
        echo "[skip] $file - already present"
        continue
    fi
    
    url="$BASE_URL/$file"
    echo -n "[GET] $file... "
    
    # Download with curl, follow redirects, show progress
    if curl -L -o "$out_path" --progress-bar "$url" 2>/dev/null; then
        size=$(du -h "$out_path" | cut -f1)
        echo " OK ($size)"
    else
        echo " ERROR"
        rm -f "$out_path"
        exit 1
    fi
done

echo ""
echo "Done. Files in $DEST_DIR :"
ls -lh "$DEST_DIR" | tail -n +2 | awk '{printf "  %-25s %10s\n", $9, $5}'
echo ""

################################################################################
# Verify SHA256 hashes
################################################################################

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Verifying SHA256 hashes..."
echo ""

# Define expected hashes (from GitHub Release v0.1-data).
# This Bash-3-compatible format also works with the default macOS shell.
HASHES=(
  "pdac_demo.bed d5c7a3c5816b70f3cb6c12b08b79940495aff5f32d765f5173c53dc46f64e984"
  "pdac_demo.bim 9399f4ab3a8929dd1feac345ad57a7c0b0a262574327a106b3371e32c50e5c0c"
  "pdac_demo.fam fb4f97191a32e108452d284054a033897d763e7acc138fb0c6a7f90a7fc129f3"
  "phenotype.txt 6eb084e67daf5e06df7212aad8f9fd8f1df0929f4641a36f95c248a5939822cc"
  "covariates.txt a0ee20ec8ca388bcf8986701c0843c236eb5da88b33d8540ad3cdf2e291561be"
  "survival.txt 6525a67c087dc5d89f528be5190a1a93f3ce6d635e861170d29ac5046eea69bb"
  "sample_ancestry.tsv 9d78bb59aa1bb16c83c457dcf3c376fd500bd2ca55913c5690c33dae25c06e23"
)

PASS=0
FAIL=0

for hash_entry in "${HASHES[@]}"; do
  filename="${hash_entry%% *}"
  expected_hash="${hash_entry#* }"
  filepath="${DEST_DIR}/${filename}"
  
  if [ ! -f "$filepath" ]; then
    echo "❌ MISSING: $filename"
    FAIL=$((FAIL + 1))
    continue
  fi
  
  # Compute hash (support both sha256sum and shasum)
  if command -v sha256sum &> /dev/null; then
    actual_hash=$(sha256sum "$filepath" | awk '{print $1}')
  else
    actual_hash=$(shasum -a 256 "$filepath" | awk '{print $1}')
  fi
  
  if [ "$actual_hash" == "$expected_hash" ]; then
    echo "✓ OK: $filename"
    PASS=$((PASS + 1))
  else
    echo "❌ HASH MISMATCH: $filename"
    echo "   Expected: $expected_hash"
    echo "   Got:      $actual_hash"
    FAIL=$((FAIL + 1))
  fi
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ $FAIL -eq 0 ]; then
  echo "✓ All files verified successfully! ($PASS/$PASS)"
  echo ""
  echo "Ready to run QC pipeline:"
  echo "  bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh"
  exit 0
else
  echo "❌ Verification failed: $FAIL file(s) with issues"
  echo ""
  echo "Solutions:"
  echo "  1. Delete mismatched files and re-download from:"
  echo "     https://github.com/mgentiluomo/how-to-gwas-pdac/releases/tag/v0.1-data"
  echo "  2. Or run this script again: bash scripts/dev/download_demo_data.sh"
  exit 1
fi
