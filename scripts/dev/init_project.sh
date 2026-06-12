#!/bin/bash

################################################################################
# Initialize PDAC GWAS Project Structure (bash/macOS/Linux)
#
# This script creates the complete project structure with all necessary folders
#
# Usage: bash init_project.sh
#
# Creates:
#   - scripts/         (where you'll add QC scripts)
#   - demo_data/       (for downloaded dataset)
#   - tools/bin/       (for plink, plink2, R executables)
#   - data_processed/  (intermediate outputs)
#   - results/         (organized by analysis type)
#
################################################################################

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

PROJECT_ROOT="$(pwd)"

echo -e "${CYAN}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║  PDAC GWAS Project Initialization                           ║${NC}"
echo -e "${CYAN}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}Project Root: $PROJECT_ROOT${NC}"
echo ""

# Define directory structure
dirs=(
    'scripts'
    'demo_data'
    'tools/bin'
    'data_processed'
    'results/qc'
    'results/pop_structure'
    'results/imputation'
    'results/association'
    'results/finemapping'
    'results/meta_analysis'
)

echo -e "${YELLOW}Creating directory structure...${NC}"
echo ""

for dir in "${dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo -e "${GREEN}  ✓${NC} Created: $dir"
    else
        echo -e "  • Already exists: $dir"
    fi
done

echo ""
echo -e "${YELLOW}Directory structure:${NC}"
echo ""

# Show tree structure (compatible with macOS and Linux)
tree -L 2 -d 2>/dev/null || find . -maxdepth 2 -type d | sort | sed 's|^\./||' | awk '
    NR==1 { print $0; next }
    {
        depth = gsub(/\//, "/")
        indent = sprintf("%*s", depth*4, "")
        sub(/.*\//, "", $0)
        print indent "├── " $0
    }
'

echo ""
echo -e "${CYAN}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║  Next Steps:                                                 ║${NC}"
echo -e "${CYAN}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}1. Add QC scripts to scripts/ folder:${NC}"
echo -e "   ${NC}git clone https://github.com/mgentiluomo/how-to-gwas-pdac.git"
echo -e "   ${NC}cp how-to-gwas-pdac/sections/01B_genotyping_qc/scripts/*.sh scripts/"
echo ""
echo -e "${YELLOW}2. Download demo data:${NC}"
echo -e "   ${NC}bash scripts/dev/download_demo_data.sh"
echo ""
echo -e "${YELLOW}3. Install PLINK tools to tools/bin/ (see docs/helpers/SETUP_PROJECT_STRUCTURE.md):${NC}"
echo -e "   ${NC}Download PLINK2 and PLINK1.9 from official sources"
echo -e "   ${NC}Place executables in tools/bin/"
echo ""
echo -e "${YELLOW}4. Set PATH and run first QC script:${NC}"
echo -e "   ${NC}export PATH=\"\$(pwd)/tools/bin:\$PATH\""
echo -e "   ${NC}bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh"
echo ""
echo -e "${CYAN}For detailed instructions, see:${NC}"
echo -e "   ${NC}docs/helpers/SETUP_PROJECT_STRUCTURE.md"
echo ""
