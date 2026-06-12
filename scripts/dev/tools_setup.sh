#!/bin/bash

################################################################################
# WSL Setup Script for PDAC GWAS Tutorial
# 
# Run this ONCE in WSL to install all dependencies:
#   bash wsl_setup.sh
#
# You may be prompted for your WSL password (sudo access required)
#
################################################################################

set -e

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║   PDAC GWAS Tutorial — WSL Dependency Installation            ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "This script will install:"
echo "  • PLINK1.9 (genetic analysis — backward compatibility)"
echo "  • PLINK2 (genetic analysis — modern, faster)"
echo "  • R & R-dev (statistical computing)"
echo "  • Git, curl, wget (utilities)"
echo ""
echo "Requires: sudo access (you may be prompted for your password)"
echo ""

# ============================================================================
# Step 1: Update package manager
# ============================================================================
echo "Step 1: Updating package manager..."
sudo apt-get update -y >/dev/null 2>&1
echo "  ✓ Done"

# ============================================================================
# Step 2: Install PLINK2
# ============================================================================
echo "Step 2: Installing PLINK2..."
sudo apt-get install -y plink2 >/dev/null 2>&1
PLINK2_VERSION=$(plink2 --version 2>&1 | head -1)
echo "  ✓ PLINK2 installed: $PLINK2_VERSION"

# ============================================================================
# Step 3: Install PLINK1.9
# ============================================================================
echo "Step 3: Installing PLINK1.9..."
sudo apt-get install -y plink >/dev/null 2>&1
PLINK_VERSION=$(plink --version 2>&1 | head -1)
echo "  ✓ PLINK1.9 installed: $PLINK_VERSION"

# ============================================================================
# Step 4: Install R
# ============================================================================
echo "Step 4: Installing R..."
sudo apt-get install -y r-base r-base-dev >/dev/null 2>&1
R_VERSION=$(R --version | head -1)
echo "  ✓ R installed: $R_VERSION"

# ============================================================================
# Step 5: Install utilities
# ============================================================================
echo "Step 5: Installing utilities (git, curl, wget)..."
sudo apt-get install -y git curl wget >/dev/null 2>&1
echo "  ✓ Utilities installed"

# ============================================================================
# Verification
# ============================================================================
echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                    Installation Complete!                     ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Verification:"
echo "  PLINK2:"
plink2 --version | head -1 | sed 's/^/    /'
echo "  PLINK1.9:"
plink --version | head -1 | sed 's/^/    /'
echo "  R:"
R --version | head -1 | sed 's/^/    /'
echo "  curl: $(curl --version | head -1 | sed 's/^/    /')"
echo ""
echo "Next steps:"
echo "  1. Navigate to the repo:"
echo "     cd /mnt/s/Github/how-to-gwas-pdac"
echo ""
echo "  2. Download demo data:"
echo "     bash demo_dataset/download_data.sh"
echo ""
echo "  3. Run first QC script:"
echo "     cd sections/01B_genotyping_qc/scripts"
echo "     bash 01_initial_qc_stats.sh"
echo ""
