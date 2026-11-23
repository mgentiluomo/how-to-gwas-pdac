#!/usr/bin/env bash

###############################################################################
# make_lrld_bed.sh
#
# Writes data_processed/highLD_b38.bed, the list of long-range linkage
# disequilibrium regions excluded before LD pruning and principal component
# analysis. run_all.sh calls it before the first step, so the file is always
# derived from this script and never kept by hand.
#
# WHY THE REGIONS ARE EXCLUDED AT ALL
#   Pairwise pruning works inside a sliding window of a few dozen variants. The
#   regions that most distort principal components are far longer: the major
#   histocompatibility complex spans about 8 Mb, the chromosome 8 inversion
#   about 4 Mb. Pruning thins such a block but leaves one survivor per window,
#   and the survivors stay correlated with one another, so a leading component
#   can describe which haplotype an individual carries rather than ancestry.
#
# SOURCE: ALL 24 REGIONS OF PRICE ET AL. (2008)
#   Price AL, Weale ME, Patterson N, Myers SR, Need AC, Shianna KV, Ge D,
#   Rotter JI, Torres E, Taylor KD, Goldstein DB, Reich D. Long-range LD can
#   confound genome scans in admixed populations. Am J Hum Genet
#   2008;83(1):132-5. PMID 18606306. doi:10.1016/j.ajhg.2008.06.005
#
#   The letter lists 24 autosomal regions of long-range LD found by PCA in
#   European samples: 3 in its Table 1 (6: 25.5-33.5, 8: 112-115,
#   11: 46-57 Mb) and 21 in the text under the table. Coordinates are given in
#   Mb at 0.5 Mb resolution. They are copied below exactly as printed
#   (PRICE_2008_MB), so any reader can check them against the article.
#
# SOURCE BUILD: NCBI36 (hg18), AN INFERENCE, NOT A STATEMENT OF THE ARTICLE
#   The letter does not name its genome build. NCBI36 is adopted because:
#     1. GRCh37 did not exist in 2008; the candidates are NCBI35 and NCBI36.
#     2. The three SNPs of Table 1 (rs169679 at 29.0 Mb, rs896760 at 113.5 Mb,
#        rs637249 at 56.0 Mb) sit at 28.96, 113.59 and 56.00 Mb in both
#        NCBI35 and NCBI36 (Ensembl assembly mapping from their GRCh38
#        positions), consistent with the printed values at 0.5 Mb
#        resolution; the SNPs do not distinguish the two builds.
#     3. Read as NCBI36, the 18 regions that the plinkQC GRCh38 list retains
#        map to exactly the coordinates of that list (see COMPARISON below).
#        plinkQC used UCSC liftOver, an independent tool, so this is a check
#        of the mapping as well as evidence of how the list is usually read.
#   Read as NCBI35 instead, six regions (1: 48-52; 2: 86-100.5, 134.5-138,
#   183-190; 7: 55-66; 12: 109.5-112) would move by 0.02 to 0.19 Mb, less
#   than the 0.5 Mb resolution of the source.
#
# LIFTOVER TO GRCh38
#   Tool: Ensembl REST assembly mapper, endpoint
#     https://rest.ensembl.org/map/human/NCBI36/{chr}:{pos}..{pos}:1/GRCh38
#   Ensembl release 116 (GET /info/software), queried 23 November 2026.
#   No UCSC chain file is used. Each region boundary is mapped as a single
#   base, not the region as one interval: interval liftover of multi-Mb
#   regions fails or fragments where the assemblies differ inside the region
#   (around centromeres in particular), which is how six of the 24 regions
#   were lost from the plinkQC list. Every boundary mapped to a single
#   position on the same chromosome, and every region kept start < end.
#   Region lengths change where GRCh38 resized a centromeric gap (for
#   example 3: 89-97.5 becomes 7.43 Mb, 12: 33-40 becomes 8.36 Mb).
#   Conversion: Price start S Mb -> 1-based position S x 1e6 in NCBI36 ->
#   mapped GRCh38 position P -> BED0 start P - 1; Price end E Mb -> E x 1e6
#   -> mapped position Q -> BED0 end Q.
#   The mapping is embedded below (GRCh38 table), so the canonical run needs
#   no network. `bash scripts/dev/make_lrld_bed.sh --derive` repeats the
#   queries and fails if any boundary no longer maps to the embedded value.
#
# COMPARISON WITH THE plinkQC GRCh38 LIST
#   plinkQC (https://github.com/meyer-lab-cshl/plinkQC), package version
#   1.1.0, file inst/extdata/high-LD-regions-hg38-GRCh38.txt as on branch
#   master at commit 2e64f66a473d22efd95c5e8ae9cea77329fa5567 (13 October
#   2026); the file itself last changed in commit
#   43a07bdd330c42533ced961df6e25590dd15f45e (9 July 2019). 43 entries, md5
#   f41e67afac5200e4e77681630c3e9387. Consulted 23 November 2026. This
#   pipeline used that list, filtered to 19 regions, until November 2026.
#     - 17 of the 24 Price regions appear there with exactly the coordinates
#       below; 6: 57-64 appears as two pieces (57.79-58.45 and 61.11-61.36
#       Mb) inside the 57.03-63.23 Mb region below.
#     - 6 Price regions are absent: 2: 134.5-138, 3: 89-97.5, 5: 98-100.5,
#       5: 135.5-138.5, 11: 46-57 (one of the three of Table 1) and
#       12: 109.5-112. In their place the list carries fragments of 22 to
#       557 bp, several on other chromosomes, which are what remains of
#       intervals that liftOver could not map whole.
#     - The earlier version of this script dropped entries shorter than
#       100 kb as liftOver fragments. That was correct of the fragments but
#       missed what they meant: six regions had been lost upstream.
#   On the demonstration data (results/qc/pdac_demo_08_filt.bim) the 19
#   regions held 13,044 variants and the 24 below hold 15,200.
#
#   Coordinates are BED0: half-open, zero-based start, as PLINK 2 reads them
#   with --exclude bed0. Chromosome codes are bare integers to match the .bim
#   of this dataset; a file using chr-prefixed codes would exclude nothing and
#   the guards in the consuming scripts would stop the run.
#
# USAGE
#   bash scripts/dev/make_lrld_bed.sh [OUT_BED] [BIM_TO_CHECK]
#   bash scripts/dev/make_lrld_bed.sh --derive     (needs network, python3)
#
#   OUT_BED        default data_processed/highLD_b38.bed
#   BIM_TO_CHECK   optional .bim; if given, reports how many of its variants
#                  fall inside the regions, the check that the file matches
#                  the build and chromosome coding of the data.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Price et al. 2008, as printed: chromosome, start Mb, end Mb (NCBI36, see above).
PRICE_2008_MB='1	48	52
2	86	100.5
2	134.5	138
2	183	190
3	47.5	50
3	83.5	87
3	89	97.5
5	44.5	50.5
5	98	100.5
5	129	132
5	135.5	138.5
6	25.5	33.5
6	57	64
6	140	142.5
7	55	66
8	8	12
8	43	50
8	112	115
10	37	43
11	46	57
11	87.5	90.5
12	33	40
12	109.5	112
20	32	34.5'

# The same 24 regions on GRCh38, BED0, from the Ensembl mapping above. Same order.
GRCH38_BED0='1	47761740	51761740
2	85919365	100517106
2	134025958	137525960
2	182427027	189427029
3	47483505	49987563
3	83368158	86868160
3	88868159	96298466
5	44464140	51168409
5	98636395	101136397
5	129636407	132636409
5	136136411	139136412
6	25391792	33424245
6	57027242	63232136
6	139637169	142137170
7	54964812	66897578
8	8105067	12105082
8	43025699	48924888
8	110918594	113918595
10	36671065	43184546
11	46021872	57475951
11	88127183	91127184
12	32955798	41319931
12	110577811	113077812
20	33948532	36438183'

if [ "${1:-}" = "--derive" ]; then
  command -v python3 >/dev/null || { echo "python3 is needed for --derive" >&2; exit 1; }
  echo "Re-deriving the GRCh38 coordinates from Price et al. 2008 (NCBI36) via Ensembl REST"
  PRICE="$PRICE_2008_MB" EMB="$GRCH38_BED0" python3 -I - <<'PY'
import os, json, urllib.request, time, sys
def get(url):
    for i in range(4):
        try:
            return json.load(urllib.request.urlopen(url, timeout=30))
        except Exception:
            time.sleep(2 * (i + 1))
    sys.exit("Ensembl REST unreachable: " + url)
rel = get("https://rest.ensembl.org/info/software?content-type=application/json")["release"]
print(f"  Ensembl release {rel}")
def mp(c, p):
    d = get(f"https://rest.ensembl.org/map/human/NCBI36/{c}:{p}..{p}:1/GRCh38?content-type=application/json")
    m = [x["mapped"] for x in d["mappings"]]
    if len(m) != 1 or m[0]["seq_region_name"] != c:
        sys.exit(f"  {c}:{p} does not map to a single position on chromosome {c}: {m}")
    return m[0]["start"]
bad = 0
for src, emb in zip(os.environ["PRICE"].splitlines(), os.environ["EMB"].splitlines()):
    c, s, e = src.split(); s1 = round(float(s) * 1e6); e1 = round(float(e) * 1e6)
    got = (c, mp(c, s1) - 1, mp(c, e1))
    exp = tuple([emb.split()[0]] + [int(x) for x in emb.split()[1:]])
    ok = got == exp; bad += not ok
    print(f"  {c}:{s}-{e} Mb -> {got[1]}-{got[2]}  {'ok' if ok else 'DIFFERS from embedded ' + str(exp[1:])}")
if bad: sys.exit(f"  {bad} region(s) differ from the embedded table")
print("  all 24 regions reproduce the embedded GRCh38 table")
PY
  exit 0
fi

OUT="${1:-data_processed/highLD_b38.bed}"
BIM="${2:-}"

mkdir -p "$(dirname "$OUT")"
printf '%s\n' "$GRCH38_BED0" > "$OUT"

N_REGIONS=$(awk 'NF>=3' "$OUT" | wc -l)
COVERAGE=$(awk '{t+=$3-$2} END{printf "%.1f", t/1e6}' "$OUT")

echo ""
echo "=== Long-range LD regions ==="
echo ""
echo "  wrote      $OUT"
echo "  regions    $N_REGIONS (Price et al. 2008, NCBI36 -> GRCh38)"
echo "  coverage   ${COVERAGE} Mb"
echo "  build      GRCh38, BED0, bare chromosome codes"

# Self-checks. A malformed region file is worse than a missing one, because the
# consuming steps would run and produce plausible output.
[ "$N_REGIONS" -eq 24 ] || { echo "  expected 24 regions, wrote $N_REGIONS" >&2; exit 1; }
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
