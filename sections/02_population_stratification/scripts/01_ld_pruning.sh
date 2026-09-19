#!/usr/bin/env bash

################################################################################
# Section 2: Population stratification — Step 01: LD pruning
#
# PURPOSE:
#   Principal component analysis must be computed on variants that are
#   approximately independent of one another. Variants in linkage disequilibrium
#   (LD) carry overlapping information, and dense regions such as the major
#   histocompatibility complex on chromosome 6 would otherwise dominate the
#   components, so that the leading axes describe local genomic architecture
#   rather than ancestry.
#
#   This step keeps roughly one representative variant per LD block.
#
# WHY --indep-pairwise IS NOT ENOUGH ON ITS OWN
#   Pairwise pruning works inside a sliding window of 50 variants. The regions
#   that most distort principal components are far longer than that: the MHC
#   spans roughly 8 Mb, the chromosome 8 inversion 5 Mb, the lactase region
#   4 Mb. Pruning thins those regions but leaves a representative variant per
#   window across the whole block, and those survivors remain correlated with
#   one another. They then behave as a single very heavily weighted locus, and
#   a leading component can end up describing which haplotype an individual
#   carries at that locus rather than their ancestry.
#
#   The regions are therefore removed outright, before pruning, with
#   --exclude bed0. This is the step the Methods text refers to as excluding
#   known long-range LD regions.
#
# INPUT:
#   - results/qc/pdac_demo_08_filt.bed/bim/fam   (QC-passed data, Section 1B)
#   - data_processed/highLD_b38.bed              (long-range LD regions, GRCh38)
#
# OUTPUT:
#   - results/pca/pdac_demo_02_prune.prune.in    (variants to keep)
#   - results/pca/pdac_demo_02_prune.prune.out   (variants pruned away)
#
# PARAMETERS (--indep-pairwise 50 5 0.2):
#   50   window size, in variants
#   5    step size: the window advances five variants at a time
#   0.2  r-squared threshold; within a window, one of any pair above this is
#        dropped
#
#   The r-squared threshold is a genuine decision, not a standard. 0.1 is more
#   stringent (|r| about 0.32), 0.2 more permissive (|r| about 0.45). What is
#   required is to report the value used and to confirm that the leading
#   components do not change materially when it is varied.
#
#   --maf applies a frequency floor before pruning. Low-frequency variants
#   contribute little to ancestry axes and their LD is poorly estimated in a
#   sample of this size. Note the floor is evaluated on the pooled cohort here,
#   because this PCA is the one that assigns ancestry and therefore cannot be
#   stratified by it; the within-group PCA in Step 05 applies its own.
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

IN="${1:-results/qc/pdac_demo_08_filt}"
OUT_DIR="${2:-results/pca}"
LRLD_BED="${3:-data_processed/highLD_b38.bed}"
PRUNE_MAF="${4:-0.05}"
OUT_PREFIX="${OUT_DIR}/pdac_demo_02_prune"

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Long-range LD guard. Fails loudly rather than proceeding without the
# exclusion, because a silent no-op produces output indistinguishable from a
# correct run and would put the Methods text out of step with the code.
# ---------------------------------------------------------------------------
require_lrld() {
  local bed="$1" bim="$2" n regions
  if [ ! -s "$bed" ]; then
    echo "✗ Long-range LD region file not found: $bed" >&2
    echo "  Create it before running this step, or pass a path as argument 3." >&2
    exit 1
  fi
  n=$(awk '
    FNR == NR {
      if ($0 ~ /^#/ || NF < 3) next
      k++; c[k] = $1; s[k] = $2 + 0; e[k] = $3 + 0
      next
    }
    { for (i = 1; i <= k; i++) if ($1 == c[i] && $4 > s[i] && $4 <= e[i]) { hit++; break } }
    END { print hit + 0 }
  ' "$bed" "$bim")
  regions=$(awk '!/^#/ && NF >= 3' "$bed" | wc -l)
  if [ "$n" -eq 0 ]; then
    echo "✗ No variants fall inside the long-range LD regions listed in $bed" >&2
    echo "  The exclusion would be a silent no-op. Usual causes:" >&2
    echo "    - chromosome codes differ ('1' in the .bim vs 'chr1' in the BED)" >&2
    echo "    - the BED is on a different build than the data (GRCh38 expected)" >&2
    echo "    - bed0 vs bed1 coordinate convention mismatch" >&2
    exit 1
  fi
  echo "Long-range LD exclusion: ${regions} regions, ${n} variants removed before pruning"
}

echo ""
echo "=== LD pruning ==="
echo ""

require_lrld "$LRLD_BED" "${IN}.bim"
echo ""

plink2 \
  --bfile "$IN" \
  --autosome \
  --maf "$PRUNE_MAF" \
  --exclude bed0 "$LRLD_BED" \
  --indep-pairwise 50 5 0.2 \
  --out "$OUT_PREFIX"

KEPT=$(wc -l < "${OUT_PREFIX}.prune.in")
DROPPED=$(wc -l < "${OUT_PREFIX}.prune.out")

echo ""
echo "Variants retained for PCA:  $KEPT"
echo "Variants pruned away:       $DROPPED"
echo ""
echo "A substantial fraction is expected to be pruned: this is redundancy being"
echo "removed, not data being lost. The pruned set is used only for PCA and"
echo "relatedness, never for association testing."
echo ""
echo "Note that prune.in + prune.out is smaller than the input, because the"
echo "autosome restriction, the MAF floor and the long-range LD regions are"
echo "applied before pruning begins and those variants never enter either list."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  bash scripts/02_population_stratification/02_pca_all.sh"
echo ""
