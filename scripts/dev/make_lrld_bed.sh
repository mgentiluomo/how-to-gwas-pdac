#!/usr/bin/env bash

###############################################################################
# make_lrld_bed.sh
#
# Writes data_processed/highLD_b38.bed, the list of long-range linkage
# disequilibrium regions excluded before LD pruning, principal component
# analysis and relatedness estimation.
#
# WHY THIS SCRIPT EXISTS
#   Five steps of this pipeline require this file, and data_processed/ is not
#   under version control. Before this script existed the file was made by hand
#   on one machine, so a reader who cloned the repository and followed the
#   instructions stopped at section 2 step 1, and the canonical run could not be
#   reproduced anywhere else. The coordinates now live in version control and
#   the file is derived from them deterministically.
#
# WHY THE REGIONS ARE EXCLUDED AT ALL
#   Pairwise pruning works inside a sliding window of a few dozen variants. The
#   regions that most distort principal components are far longer: the major
#   histocompatibility complex spans about 8 Mb, the chromosome 8 inversion
#   about 4 Mb. Pruning thins such a block but leaves one survivor per window,
#   and the survivors stay correlated with one another, so they behave as a
#   single heavily weighted locus and a leading component can describe which
#   haplotype an individual carries rather than their ancestry.
#   Price AL, Weale ME, Patterson N, Myers SR, Need AC, Shianna KV, Ge D,
#   Rotter JI, Torres E, Taylor KD, Goldstein DB, Reich D. Long-range LD can
#   confound genome scans in admixed populations. Am J Hum Genet
#   2008;83(1):132-5. PMID 18606306. doi:10.1016/j.ajhg.2008.06.005
#
# PROVENANCE OF THE COORDINATES
#   GRCh38 coordinates as distributed with the plinkQC R package, file
#   inst/extdata/high-LD-regions-hg38-GRCh38.txt
#   (https://github.com/meyer-lab-cshl/plinkQC), which carries the long-range
#   LD region list lifted to GRCh38.
#
#   Two documented edits were made to that file, both inspectable here:
#     1. Entries shorter than 100 kb were dropped. The distributed file
#        contains 43 entries, 20 of them shorter than 1 kb and one of 1.9 kb.
#        Those are liftOver fragments, not long-range LD regions: an interval
#        of 22 bp cannot confound a principal component. The retained regions
#        span 248 kb to 14.6 Mb, with no entry between 1.9 kb and 248 kb, so
#        the cut falls in an empty part of the distribution. None of the
#        dropped entries contains a variant in the demonstration dataset, so
#        on these data the edit changes nothing; it is made so that the file a
#        reader inspects contains only regions that are what they claim to be.
#     2. Three regions nested inside larger ones were merged into them, so that
#        each interval appears once. Merging cannot change which variants are
#        excluded, and it was verified not to.
#
#   Coordinates are BED0: half-open, zero-based start, as PLINK 2 reads them
#   with --exclude bed0. Chromosome codes are bare integers to match the .bim
#   of this dataset; a file using chr-prefixed codes would exclude nothing and
#   the guards in the consuming scripts would stop the run.
#
# USAGE
#   bash scripts/dev/make_lrld_bed.sh [OUT_BED] [BIM_TO_CHECK]
#
#   OUT_BED        default data_processed/highLD_b38.bed
#   BIM_TO_CHECK   optional .bim; if given, reports how many of its variants
#                  fall inside the regions, which is the check that the file
#                  matches the build and the chromosome coding of the data.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

OUT="${1:-data_processed/highLD_b38.bed}"
BIM="${2:-}"

mkdir -p "$(dirname "$OUT")"

cat > "$OUT" <<'REGIONS'
1	47761740	51761740
2	85919365	100517106
2	182427027	189427029
3	47483505	49987563
3	83368158	86868160
5	44464140	51168409
5	129636407	132636409
6	25391792	33424245
6	57788603	58453888
6	61109122	61357029
6	139637169	142137170
7	54964812	66897578
8	8105067	12105082
8	43025699	48924888
8	110918594	113918595
10	36671065	43184546
11	88127183	91127184
12	32955798	41319931
20	33948532	36438183
REGIONS

N_REGIONS=$(awk 'NF>=3' "$OUT" | wc -l)
COVERAGE=$(awk '{t+=$3-$2} END{printf "%.1f", t/1e6}' "$OUT")

echo ""
echo "=== Long-range LD regions ==="
echo ""
echo "  wrote      $OUT"
echo "  regions    $N_REGIONS"
echo "  coverage   ${COVERAGE} Mb"
echo "  build      GRCh38, BED0, bare chromosome codes"

# Self-checks. A malformed region file is worse than a missing one, because the
# consuming steps would run and produce plausible output.
awk 'NF<3 { print "  MALFORMED line " NR; bad=1 } $2+0>=$3+0 { print "  NON-POSITIVE width, line " NR; bad=1 } END{ exit bad+0 }' "$OUT" \
  || { echo "  the written file is malformed" >&2; exit 1; }

if [ -n "$BIM" ]; then
  if [ ! -s "$BIM" ]; then
    echo "  cannot check against $BIM: not found" >&2
    exit 1
  fi
  awk -v bim="$BIM" '
    FNR==NR { if ($0 ~ /^#/ || NF < 3) next; k++; c[k]=$1; s[k]=$2+0; e[k]=$3+0; next }
    { tot++; for (i=1;i<=k;i++) if ($1==c[i] && $4>s[i] && $4<=e[i]) { hit++; break } }
    END {
      printf "  checked    %s\n", bim
      printf "  variants   %d of %d fall inside the regions (%.2f%%)\n", hit+0, tot, 100*(hit+0)/tot
      if (hit+0 == 0) {
        print "  the exclusion would be a silent no-op: check build and chromosome coding" > "/dev/stderr"
        exit 1
      }
    }
  ' "$OUT" "$BIM"
fi

echo ""
