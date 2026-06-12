#!/bin/bash

# GWAS Tutorial Setup Test
# Tests if folders are created, demo data downloaded, and tools are working

echo ""
echo "=================================="
echo "GWAS Tutorial Setup Test"
echo "=================================="
echo ""

# Test folders
echo "1. Checking folder structure..."
FOLDERS=("scripts" "demo_data" "tools/bin" "data_processed" "results")
FOLDER_OK=true

for folder in "${FOLDERS[@]}"; do
    if [ -d "$folder" ]; then
        echo "   ✓ $folder"
    else
        echo "   ✗ $folder (missing)"
        FOLDER_OK=false
    fi
done

# Test demo data files
echo ""
echo "2. Checking demo dataset files..."
DATA_FILES=("pdac_demo.bed" "pdac_demo.bim" "pdac_demo.fam" "phenotype.txt" "covariates.txt" "survival.txt" "sample_ancestry.tsv")
DATA_OK=true

for file in "${DATA_FILES[@]}"; do
    if [ -f "demo_data/$file" ]; then
        echo "   ✓ $file"
    else
        echo "   ✗ $file (missing)"
        DATA_OK=false
    fi
done

# Test tools
echo ""
echo "3. Checking installed tools..."
TOOLS_OK=true

# Test PLINK2
if command -v plink2 &> /dev/null; then
    PLINK2_VERSION=$(plink2 --version 2>&1 | head -1)
    echo "   ✓ plink2 ($PLINK2_VERSION)"
else
    echo "   ✗ plink2 (not found)"
    TOOLS_OK=false
fi

# Test PLINK 1.9
if command -v plink &> /dev/null; then
    PLINK_VERSION=$(plink --version 2>&1 | head -1)
    echo "   ✓ plink ($PLINK_VERSION)"
else
    echo "   ✗ plink (not found)"
    TOOLS_OK=false
fi

# Test R
if command -v R &> /dev/null; then
    R_VERSION=$(R --version 2>&1 | head -1)
    echo "   ✓ R ($R_VERSION)"
fi

# Test wget
if command -v wget &> /dev/null; then
    echo "   ✓ wget"
fi

# Test git
if command -v git &> /dev/null; then
    echo "   ✓ git"
fi

# Summary
echo ""
echo "=================================="
echo "Summary"
echo "=================================="

if [ "$FOLDER_OK" = true ] && [ "$DATA_OK" = true ] && [ "$TOOLS_OK" = true ]; then
    echo "✓ Initial test done! All folders, data, and tools are working."
    echo ""
    echo "Next command to run the QC pipeline:"
    echo "  bash scripts/01_initial_qc_stats.sh"
    echo ""
    exit 0
else
    echo "✗ Some issues found. Please check the output above and fix them."
    echo ""
    if [ "$FOLDER_OK" = false ]; then
        echo "  • Missing folders: Run Step 2 again"
    fi
    if [ "$DATA_OK" = false ]; then
        echo "  • Missing data files: Run Step 4 again"
    fi
    if [ "$TOOLS_OK" = false ]; then
        echo "  • Missing tools: Run Step 5 again"
    fi
    echo ""
    exit 1
fi
