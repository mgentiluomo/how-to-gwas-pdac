#!/usr/bin/env bash

################################################################################
# Section 1A: Genotyping technologies — Step 04: pre-flight checks and formats
#
# PURPOSE:
#   Two things before any analysis begins.
#
#   First, the checks that catch the errors which are cheap to fix now and
#   expensive to fix later: unplaced variants, duplicated positions,
#   strand-ambiguous variants, and the genome build.
#
#   Second, the format conversions you will actually need: PLINK 2 native
#   format for speed, VCF for imputation servers and for sharing, text for
#   inspection.
#
# INPUT:
#   - results/raw/merged/DB_merged_final    (the converted example, Step 03)
#   - results/qc/pdac_demo_08_filt or demo_data/pdac_demo   (the demo dataset)
#
# OUTPUT:
#   - results/raw/pdac_demo_01A_preflight.txt
#   - results/raw/format_demo.*             (conversion examples)
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "$SCRIPT_DIR/../../scripts/dev" ]; then
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
elif [ -d "$SCRIPT_DIR/../../../scripts/dev" ]; then
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
else
  PROJECT_ROOT="$(pwd)"
fi
cd "$PROJECT_ROOT"

OUT_DIR="${1:-results/raw}"
EXAMPLE="${OUT_DIR}/merged/DB_merged_final"
DEMO="${2:-demo_data/pdac_demo}"
REPORT="${OUT_DIR}/pdac_demo_01A_preflight.txt"

mkdir -p "$OUT_DIR"

{
echo "Pre-flight checks"
echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# --- 1. unplaced variants ----------------------------------------------------
UNPLACED=$(awk '$1 == 0 || $4 == 0' "${EXAMPLE}.bim" | wc -l)
TOTAL=$(wc -l < "${EXAMPLE}.bim")
echo "1. UNPLACED VARIANTS"
echo "   converted example: ${UNPLACED} of ${TOTAL} have no chromosome or position"
echo "   These cannot be merged with another dataset, lifted between builds, or"
echo "   imputed. They are not junk: their coordinates live in the array"
echo "   manifest, which must be used during conversion."
echo ""

# --- 2. duplicated positions -------------------------------------------------
DUP_POS=$(awk '$1 != 0 && $4 != 0 {print $1":"$4}' "${EXAMPLE}.bim" | sort | uniq -d | wc -l)
DUP_ID=$(cut -f2 "${EXAMPLE}.bim" | sort | uniq -d | wc -l)
echo "2. DUPLICATED POSITIONS AND IDENTIFIERS"
echo "   duplicated chr:pos: ${DUP_POS}"
echo "   duplicated variant IDs: ${DUP_ID}"
echo "   Arrays deliberately place several probes at the same coordinate, for"
echo "   example when tiling a gene for copy-number calling. Merging or"
echo "   imputing on position without resolving these produces silent errors."
echo ""

# --- 3. strand-ambiguous variants --------------------------------------------
AMB=$(awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") ||
           ($5=="C" && $6=="G") || ($5=="G" && $6=="C")' "${EXAMPLE}.bim" | wc -l)
echo "3. STRAND-AMBIGUOUS VARIANTS"
echo "   A/T or C/G in the converted example: ${AMB} of ${TOTAL}"
echo "   For these the strand cannot be resolved by comparing alleles, because"
echo "   A complements T and C complements G. They are the classic source of"
echo "   sign errors when combining studies (Section 5)."
echo ""

# --- 4. the demonstration dataset, for comparison ----------------------------
if [ -f "${DEMO}.bim" ]; then
  D_TOTAL=$(wc -l < "${DEMO}.bim")
  D_UNPLACED=$(awk '$1 == 0 || $4 == 0' "${DEMO}.bim" | wc -l)
  D_DUP=$(awk '{print $1":"$4}' "${DEMO}.bim" | sort | uniq -d | wc -l)
  D_AMB=$(awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") ||
               ($5=="C" && $6=="G") || ($5=="G" && $6=="C")' "${DEMO}.bim" | wc -l)
  echo "4. THE SAME CHECKS ON THE DEMONSTRATION DATASET"
  echo "   variants: ${D_TOTAL}"
  echo "   unplaced: ${D_UNPLACED}"
  echo "   duplicated positions: ${D_DUP}"
  echo "   strand-ambiguous: ${D_AMB}"
  echo "   The demonstration data were prepared from a curated reference panel,"
  echo "   so they carry no unplaced variants. They are not free of the other"
  echo "   problems: 22 positions are still duplicated and 2,430 variants remain"
  echo "   strand-ambiguous. Even curated data need these checks."
  echo ""
fi

# --- 5. genome build ---------------------------------------------------------
echo "5. GENOME BUILD"
echo "   Nothing in a .bim file records the build. It has to be established"
echo "   from the manifest and then written down. A GRCh37 dataset analysed as"
echo "   GRCh38 produces coordinates that are wrong by megabases in places, and"
echo "   the analysis will run without complaint. Everything in this guide is"
echo "   GRCh38, with variant identifiers in chr:pos:ref:alt form."
} | tee "$REPORT"

echo ""
echo "=== Format conversions you will actually use ==="
echo ""

if [ -f "${DEMO}.bed" ]; then
  echo "--- PLINK 1 binary to PLINK 2 native (pgen), the fastest working format"
  plink2 --bfile "$DEMO" --chr 22 --make-pgen --out "${OUT_DIR}/format_demo" > /dev/null 2>&1
  echo "    plink2 --bfile ${DEMO} --make-pgen --out out"

  echo "--- PLINK to VCF, required by imputation servers and for sharing"
  plink2 --bfile "$DEMO" --chr 22 --export vcf bgz --out "${OUT_DIR}/format_demo" > /dev/null 2>&1
  echo "    plink2 --bfile ${DEMO} --export vcf bgz --out out"

  echo "--- PLINK binary to text, for looking at the data by eye"
  plink2 --bfile "$DEMO" --chr 22 --export ped --out "${OUT_DIR}/format_demo" > /dev/null 2>&1
  echo "    plink2 --bfile ${DEMO} --export ped --out out"

  echo ""
  echo "Sizes for chromosome 22 of the demonstration dataset:"
  ls -l "${OUT_DIR}"/format_demo.* 2>/dev/null | awk '{printf "  %10d  %s\n", $5, $9}'
fi

echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  Section 1B: bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh"
echo ""
