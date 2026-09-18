#!/usr/bin/env bash

################################################################################
# Section 1B: Genotyping QC — Step 04: Heterozygosity outliers, WITHIN ANCESTRY
#
#   F is computed on an LD-pruned, common, autosomal
#   variant set, pruned separately within each ancestry group, with long-range
#   LD regions excluded. The method-of-moments F assumes the contributing
#   variants are approximately independent; computed on all 424,000 autosomal
#   variants it is dominated by whichever LD blocks happen to be densest on the
#   array, and the +/- 3 SD interval inherits that arbitrariness. Pruning within
#   the group rather than pooled is deliberate: which variants are correlated
#   depends on the population's allele frequencies, so a pooled pruned set is
#   not independent inside any single group.
#
#   The within-group MAF floor also removes the monomorphic variants that the
#   unpruned version reported (50,265 in AFR, 82,704 in EAS, 4,203 in EUR).
#   Those contribute nothing to observed heterozygosity but do enter the
#   expected-heterozygosity denominator, biasing F by an amount that differs
#   between groups precisely because the number of monomorphic sites differs.
#
# INPUT   pdac_demo_03_filt.bed/bim/fam, demo_data/sample_ancestry.tsv,
#         data_processed/highLD_b38.bed
# OUTPUT  pdac_demo_04_prune_<group>.prune.in, pdac_demo_04_het_<group>.het,
#         pdac_demo_04_het_outliers.txt, pdac_demo_04_het_summary.tsv,
#         pdac_demo_04_filt.bed/bim/fam
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

SEED="${1:-2026}"
DATASET_INPUT="${2:-results/qc}"
DATASET_NAME="pdac_demo"
ANCESTRY_FILE="${3:-demo_data/sample_ancestry.tsv}"
OUT_DIR="${4:-results/qc}"
SD_MULT="${5:-3}"
LRLD_BED="${6:-data_processed/highLD_b38.bed}"
PRUNE_MAF="${7:-0.05}"

mkdir -p "$OUT_DIR"

if [ ! -s "$ANCESTRY_FILE" ]; then
  echo "ERROR: ancestry file not found: $ANCESTRY_FILE" >&2
  echo "  Step 04 is ancestry-aware. In a real study the groups come from a" >&2
  echo "  preliminary PCA; here they ship with the demonstration data." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Long-range LD guard
#
#   --indep-pairwise cannot remove long-range LD regions on its own: a
#   50-variant window is far shorter than the MHC or the chromosome 8
#   inversion, so representatives survive from across the whole region and
#   continue to carry correlated information into whatever is computed next.
#
#   The guard fails loudly rather than proceeding without the exclusion,
#   because a silent no-op produces output that is indistinguishable from a
#   correct run and would put the Methods text out of step with the code.
# ---------------------------------------------------------------------------
require_lrld() {
  local bed="$1" bim="$2" n regions
  if [ ! -s "$bed" ]; then
    echo "✗ Long-range LD region file not found: $bed" >&2
    echo "  Create it before running this step, or pass a path as argument 6." >&2
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
  echo "Long-range LD exclusion: ${regions} regions, ${n} variants in scope"
}

require_lrld "$LRLD_BED" "${DATASET_INPUT}/${DATASET_NAME}_03_filt.bim"

# NOTE: do not name this variable GROUPS. In bash, GROUPS is a reserved array
# holding the current user's Unix groups; assigning to it is ignored and returns
# a non-zero status, which under 'set -e' terminates the script silently.
ANC_GROUPS=$(awk '{print $2}' "$ANCESTRY_FILE" | sort -u | grep -v '^$')
echo "=== Heterozygosity within ancestry groups: $(echo "$ANC_GROUPS" | tr '\n' ' ')"

: > "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt"
printf "group\tn\tmean_F\tsd_F\tlower\tupper\toutliers\tn_variants\n" \
  > "${OUT_DIR}/${DATASET_NAME}_04_het_summary.tsv"

for G in $ANC_GROUPS; do
  awk -v g="$G" '$2 == g { print $1"\t"$1 }' "$ANCESTRY_FILE" \
    > "${OUT_DIR}/${DATASET_NAME}_04_keep_${G}.txt"

  # Prune inside the group. --keep is applied before --maf and before
  # --indep-pairwise, so both the frequency floor and the LD estimates use that
  # group's own genotypes. --autosome makes explicit what PLINK was previously
  # doing silently when it excluded the 6,000 chrX variants from --het.
  plink2 \
    --bfile "${DATASET_INPUT}/${DATASET_NAME}_03_filt" \
    --keep "${OUT_DIR}/${DATASET_NAME}_04_keep_${G}.txt" \
    --autosome \
    --maf "$PRUNE_MAF" \
    --exclude bed0 "$LRLD_BED" \
    --indep-pairwise 50 5 0.2 \
    --out "${OUT_DIR}/${DATASET_NAME}_04_prune_${G}"

  NVAR=$(wc -l < "${OUT_DIR}/${DATASET_NAME}_04_prune_${G}.prune.in")
  echo "  ${G}: ${NVAR} pruned autosomal variants for F"

  # --het inside the group, so BOTH observed and expected heterozygosity use
  # that group's own allele frequencies. Using pooled frequencies would leave
  # the Wahlund inflation in F even after stratifying the threshold.
  plink2 \
    --bfile "${DATASET_INPUT}/${DATASET_NAME}_03_filt" \
    --keep "${OUT_DIR}/${DATASET_NAME}_04_keep_${G}.txt" \
    --extract "${OUT_DIR}/${DATASET_NAME}_04_prune_${G}.prune.in" \
    --het \
    --out "${OUT_DIR}/${DATASET_NAME}_04_het_${G}"

  awk -v g="$G" -v k="$SD_MULT" -v nv="$NVAR" -v out="${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt" '
    NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "F") fi = i; if ($i == "IID") ii = i } ; next }
    { n++; id[n] = $ii; f[n] = $fi + 0; s += $fi; ss += ($fi)^2 }
    END {
      m = s / n
      sd = sqrt((ss - n * m * m) / (n - 1))
      lo = m - k * sd; hi = m + k * sd
      c = 0
      for (i = 1; i <= n; i++)
        if (f[i] < lo || f[i] > hi) { print id[i]"\t"id[i] >> out; c++ }
      printf "%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n", g, n, m, sd, lo, hi, c, nv
    }
  ' "${OUT_DIR}/${DATASET_NAME}_04_het_${G}.het" \
    >> "${OUT_DIR}/${DATASET_NAME}_04_het_summary.tsv"
done

column -t "${OUT_DIR}/${DATASET_NAME}_04_het_summary.tsv"
N_OUT=$(sort -u "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt" | grep -c . || true)
echo ""
echo "Heterozygosity outliers, union across groups: ${N_OUT}"
echo "A mean F close to zero in every group is the check that this worked."
echo "n_variants differs between groups by design: each set is pruned on that"
echo "group's own LD structure, so the counts are not expected to match."

sort -u "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.txt" \
  > "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.uniq.txt"

plink2 \
  --bfile "${DATASET_INPUT}/${DATASET_NAME}_03_filt" \
  --remove "${OUT_DIR}/${DATASET_NAME}_04_het_outliers.uniq.txt" \
  --make-bed \
  --out "${OUT_DIR}/${DATASET_NAME}_04_filt"

echo ""
echo "[NEXT] bash scripts/01B_genotyping_qc/05_variant_callrate.sh"
