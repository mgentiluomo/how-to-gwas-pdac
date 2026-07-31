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
TAG='v0.3-data'
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

# Define expected hashes (from GitHub Release v0.3-data). The three genotype
# files and sample_ancestry.tsv are unchanged from v0.1-data; the phenotype,
# covariates and survival files were regenerated and their hashes differ.
# This Bash-3-compatible format also works with the default macOS shell.
HASHES=(
  "pdac_demo.bed d5c7a3c5816b70f3cb6c12b08b79940495aff5f32d765f5173c53dc46f64e984"
  "pdac_demo.bim 9399f4ab3a8929dd1feac345ad57a7c0b0a262574327a106b3371e32c50e5c0c"
  "pdac_demo.fam fb4f97191a32e108452d284054a033897d763e7acc138fb0c6a7f90a7fc129f3"
  "phenotype.txt e1b54b0519d84200aa91fa354153ae789c23c187a12d23415bdb0e7c30f60364"
  "covariates.txt 70201ed41b39bf81ac6537d5b672262d914ed7a8d026caebc9947b721cadd735"
  "survival.txt d68493ee84ae153cf494f818d69bf425489a048bf35e94b6ecc8b05c97163685"
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
  echo "Please continue getting started:"
  echo "  Step 5 — Install dependencies"
  exit 0
else
  echo "❌ Verification failed: $FAIL file(s) with issues"
  echo ""
  echo "Solutions:"
  echo "  1. Delete mismatched files and re-download from:"
  echo "     https://github.com/mgentiluomo/how-to-gwas-pdac/releases/tag/v0.3-data"
  echo "  2. Or run this script again: bash scripts/dev/download_demo_data.sh"
  exit 1
fi
