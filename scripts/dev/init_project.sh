#!/usr/bin/env bash

################################################################################
# Initialize PDAC GWAS Project Structure (bash/macOS/Linux)
#
# This script creates the complete project structure with all necessary folders
#
# Usage: bash init_project.sh
#
# Creates:
#   - scripts/         (where section scripts are copied)
#   - scripts/dev/     (utility scripts)
#   - demo_data/       (for downloaded dataset)
#   - tools/bin/       (for PLINK tool links)
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
    'scripts/dev'
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
echo -e "${YELLOW}1. Download and organize all scripts by section:${NC}"
echo -e "   ${NC}git clone https://github.com/mgentiluomo/how-to-gwas-pdac.git"
echo -e "   ${NC}find how-to-gwas-pdac/scripts/dev -maxdepth 1 -type f -exec cp {} scripts/dev/ \\;"
echo -e "   ${NC}: > scripts/dev/section_manifest.txt"
echo -e "   ${NC}: > scripts/dev/script_manifest.txt"
echo -e "   ${NC}for section_dir in how-to-gwas-pdac/sections/*/; do"
echo -e "   ${NC}  section=\$(basename \"\$section_dir\")"
echo -e "   ${NC}  if [ -d \"\$section_dir/scripts\" ]; then"
echo -e "   ${NC}    mkdir -p \"scripts/\$section\""
echo -e "   ${NC}    echo \"scripts/\$section\" >> scripts/dev/section_manifest.txt"
echo -e "   ${NC}    find \"\$section_dir/scripts\" -maxdepth 1 -type f -exec cp {} \"scripts/\$section/\" \\;"
echo -e "   ${NC}    find \"scripts/\$section\" -maxdepth 1 -type f | sort >> scripts/dev/script_manifest.txt"
echo -e "   ${NC}  fi"
echo -e "   ${NC}done"
echo -e "   ${NC}find scripts/dev -maxdepth 1 -type f -name \"*.sh\" | sort >> scripts/dev/script_manifest.txt"
echo -e "   ${NC}sort -u scripts/dev/section_manifest.txt -o scripts/dev/section_manifest.txt"
echo -e "   ${NC}sort -u scripts/dev/script_manifest.txt -o scripts/dev/script_manifest.txt"
echo -e "   ${NC}rm -r how-to-gwas-pdac"
echo ""
echo -e "${YELLOW}2. Download and verify demo data:${NC}"
echo -e "   ${NC}bash scripts/dev/download_demo_data.sh"
echo ""
echo -e "${YELLOW}3. Install/check required tools and write the tool manifest:${NC}"
echo -e "   ${NC}bash scripts/dev/tools_setup.sh"
echo ""
echo -e "${YELLOW}4. Test your setup:${NC}"
echo -e "   ${NC}bash scripts/dev/test.sh"
echo ""
echo -e "${YELLOW}5. Run your first QC script:${NC}"
echo -e "   ${NC}bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh"
echo ""
echo -e "${CYAN}For detailed instructions, see:${NC}"
echo -e "   ${NC}docs/helpers/SETUP_PROJECT_STRUCTURE.md"
echo ""
