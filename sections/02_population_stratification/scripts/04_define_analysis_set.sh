#!/usr/bin/env bash

################################################################################
# Section 2: Population stratification — Step 04: defining the analysis set
#
# PURPOSE:
#   Decide which individuals enter the primary association analysis, and record
#   the decision with its consequences for case and control counts.
#
#   Three strategies are defensible, and the choice is a judgement, not a rule:
#
#     1. Restrict to the largest ancestry group and analyse it alone. Clean, and
#        the approach taken here, but it discards data.
#     2. Retain everyone and model ancestry group as a covariate. Preserves
#        sample size; assumes effect sizes are portable across populations.
#     3. Analyse each ancestry separately and meta-analyse. Neither discards
#        data nor assumes portability, and is preferable whenever each group is
#        large enough to stand alone.
#
#   For this demonstration, option 1 keeps the worked example simple and honest:
#   the non-European groups here are too small to form meaningful strata, which
#   is exactly the situation of most real rare-cancer cohorts.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.fam
#   - demo_data/sample_ancestry.tsv
#   - demo_data/phenotype.txt
#
# OUTPUT:
#   - results/pca/pdac_demo_02_eur_keep.txt          FID/IID of the analysis set
#   - results/pca/pdac_demo_02_analysis_set.tsv      counts, for reporting
#
# IMPORTANT:
#   In a real study the assignment comes from the PCA itself, or from projection
#   onto a reference panel such as 1000 Genomes or HGDP, with a threshold on the
#   estimated ancestry proportion. Here the shipped labels are used, and Step 03
#   is the check that they agree with the genetic data.
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

FAM="${1:-results/qc/pdac_demo_08_filt.fam}"
ANC="${2:-demo_data/sample_ancestry.tsv}"
PHENO="${3:-demo_data/phenotype.txt}"
OUT_DIR="${4:-results/pca}"
TARGET="eur"

mkdir -p "$OUT_DIR"

KEEP="${OUT_DIR}/pdac_demo_02_eur_keep.txt"
COUNTS="${OUT_DIR}/pdac_demo_02_analysis_set.tsv"

echo ""
echo "=== Defining the analysis set (target ancestry: ${TARGET}) ==="
echo ""

# Individuals of the target ancestry that survived QC.
awk -v target="$TARGET" '
  NR==FNR { if (tolower($2) == target) keep[$1] = 1; next }
  ($2 in keep) { print $1"\t"$2 }
' "$ANC" "$FAM" > "$KEEP"

# Composition of the whole QC-passed set, and of the analysis set.
{
  echo -e "quantity\tvalue"
  echo -e "qc_passed_total\t$(wc -l < "$FAM")"
  for g in eur afr eas; do
    n=$(awk -v g="$g" 'NR==FNR { if (tolower($2)==g) k[$1]=1; next } ($2 in k)' "$ANC" "$FAM" | wc -l)
    echo -e "qc_passed_${g}\t${n}"
  done
  echo -e "analysis_set_total\t$(wc -l < "$KEEP")"

  # Case and control counts within the analysis set. Phenotype coding is
  # 1 = control, 2 = case, the PLINK convention.
  awk 'NR==FNR { keep[$2]=1; next }
       FNR==1 { next }
       ($2 in keep) { if ($3==2) ca++; else if ($3==1) co++ }
       END { print "analysis_set_cases\t" ca+0;
             print "analysis_set_controls\t" co+0;
             if (ca+co > 0) printf "analysis_set_effective_n\t%.1f\n", 4*ca*co/(ca+co) }
  ' "$KEEP" "$PHENO"
} > "$COUNTS"

cat "$COUNTS"

echo ""
echo "The effective sample size, 4 x cases x controls / (cases + controls), is"
echo "the number that governs power, not the headline total. Note how much"
echo "smaller it is than the raw count when the design is unbalanced."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  bash scripts/02_population_stratification/05_pca_within_eur.sh"
echo ""
