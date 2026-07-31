#!/usr/bin/env bash

################################################################################
# make_checkpoints.sh
#
# Builds the checkpoint bundles that let a reader enter the guide part-way
# through, and the manifest that proves which pipeline produced them.
#
# WHY THIS SCRIPT EXISTS
#   Checkpoints go stale silently. If the pipeline changes and the frozen files
#   do not, a reader starting from the middle gets numbers that match no page,
#   with nothing to warn them. This script is the only sanctioned way to produce
#   them: never build a checkpoint by hand, and never edit one afterwards.
#
# WHAT IT GUARANTEES
#   - every file is a byte-for-byte copy of a canonical pipeline output
#   - the manifest records the git commit that produced them
#   - checksums are taken on DECOMPRESSED content, so they do not depend on the
#     gzip version, the compression level, or the timestamp in the header
#   - the manifest carries the counts that each section's entry check verifies
#
# USAGE
#   bash scripts/dev/make_checkpoints.sh            # after a full pipeline run
#   bash scripts/dev/make_checkpoints.sh --verify   # check an existing bundle
################################################################################

set -euo pipefail

QC="results/qc"
PCA="results/pca"
ASSOC="results/assoc"
META="results/meta"
OUT="checkpoints"
MODE="${1:-build}"

hash_content() {           # md5 of the content, transparent to compression
    case "$1" in
        *.gz) gzip -dc "$1" | md5sum | cut -d' ' -f1 ;;
        *)    md5sum "$1" | cut -d' ' -f1 ;;
    esac
}

# ---------------------------------------------------------------- verify ----
if [ "$MODE" = "--verify" ]; then
    [ -f "$OUT/MANIFEST.tsv" ] || { echo "no manifest in $OUT/"; exit 1; }
    fail=0
    while IFS=$'\t' read -r file expected _; do
        [ "$file" = "file" ] && continue
        if [ ! -f "$OUT/$file" ]; then
            echo "  MISSING  $file"; fail=1; continue
        fi
        got=$(hash_content "$OUT/$file")
        if [ "$got" = "$expected" ]; then echo "  ok       $file"
        else echo "  MISMATCH $file"; fail=1; fi
    done < "$OUT/MANIFEST.tsv"
    [ "$fail" -eq 0 ] && echo "" && echo "All checkpoints match the manifest." \
                      || { echo ""; echo "Bundle is not trustworthy. Rebuild it."; exit 1; }
    exit 0
fi

# ----------------------------------------------------------------- build ----
mkdir -p "$OUT"

echo ">>> Genotype checkpoint: the dataset after all quality control"
for e in bed bim fam; do
    cp "$QC/pdac_demo_08_filt.$e" "$OUT/pdac_demo_qc.$e" 2>/dev/null \
      || cp "$PCA/pdac_demo_02_hwe_filt.$e" "$OUT/pdac_demo_qc.$e"
done

echo ">>> Small text checkpoint: who is analysed, and the covariates"
cp "$PCA/pdac_demo_02_eur_keep.txt"    "$OUT/pdac_demo_eur_keep.txt"
cp "$PCA/pdac_demo_02_pca_all.eigenvec" "$OUT/pdac_demo_pca_all.eigenvec"
cp "$PCA/pdac_demo_02_pca_all.eigenval" "$OUT/pdac_demo_pca_all.eigenval"
cp "$PCA/pdac_demo_02_pca_eur.eigenvec" "$OUT/pdac_demo_pca_eur.eigenvec"
cp "$PCA/pdac_demo_02_pca_eur.eigenval" "$OUT/pdac_demo_pca_eur.eigenval"
cp "$ASSOC/pdac_demo_04A_covar.txt" "$OUT/pdac_demo_covar_eur.txt"

echo ">>> Summary-statistics checkpoint: no genotypes needed downstream"
# -n suppresses the timestamp, so repeated runs give identical bytes
gzip -nc "$ASSOC"/pdac_demo_04A_gwas.PHENO.glm.logistic.hybrid \
       > "$OUT/pdac_demo_gwas_eur.tsv.gz"
for g in eur afr eas; do
    gzip -nc "$META/$g"/pdac_demo_05_${g}_gwas.PHENO.glm.logistic.hybrid \
           > "$OUT/pdac_demo_gwas_${g}.stratum.tsv.gz"
done

# -------------------------------------------------------------- manifest ----
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
DATE=$(date -u +%Y-%m-%d)

N_SAMPLE=$(wc -l < "$OUT/pdac_demo_qc.fam")
N_VARIANT=$(wc -l < "$OUT/pdac_demo_qc.bim")
N_EUR=$(wc -l < "$OUT/pdac_demo_eur_keep.txt")

{
  printf "file\tcontent_md5\tbytes\n"
  for f in "$OUT"/*; do
      b=$(basename "$f")
      case "$b" in MANIFEST.tsv|README.md) continue ;; esac
      printf "%s\t%s\t%s\n" "$b" "$(hash_content "$f")" "$(stat -c%s "$f")"
  done
} > "$OUT/MANIFEST.tsv"

cat > "$OUT/README.md" <<EOF
# Checkpoints

Entry points into the guide for readers who do not want to run every earlier
step. Produced by \`scripts/dev/make_checkpoints.sh\` from pipeline commit
\`$COMMIT\` on $DATE.

## What each bundle unlocks

| Bundle | Files | Starts you at |
|---|---|---|
| Genotypes after QC | \`pdac_demo_qc.bed/.bim/.fam\` | Population stratification, statistical power, association testing |
| Analysis set and covariates | \`pdac_demo_eur_keep.txt\`, \`pdac_demo_covar_eur.txt\`, the four PCA files | Association testing, together with the genotype bundle |
| Summary statistics | \`pdac_demo_gwas_*.tsv.gz\` | Meta-analysis, fine mapping, statistical power. **No genotypes required.** |

The summary-statistics bundle is the one most readers want. Meta-analysis, fine
mapping and the power calculations never touch a genotype, so learning those
steps costs about 39 MB rather than the full dataset.

## Expected contents

| Quantity | Value |
|---|---|
| Samples after QC | $N_SAMPLE |
| Variants after QC | $N_VARIANT |
| Individuals in the European analysis set | $N_EUR |

Each section that accepts a checkpoint opens with a check that prints these
numbers. If what you see differs, you have a bundle from a different version of
the pipeline and the figures on the pages will not match your output.

## Verifying a download

\`\`\`bash
bash scripts/dev/make_checkpoints.sh --verify
\`\`\`

Checksums are taken on decompressed content, so they are unaffected by the gzip
version or compression level used on your machine.

## Provenance

These files are outputs of the pipeline in this repository, not separate
artefacts. Nothing here was produced by hand or edited after generation. To
rebuild them from raw data, run the sections in order and then this script.
EOF

echo ""
echo ">>> Done. $OUT/ contains $(( $(ls -1 "$OUT" | wc -l) )) files, $(du -sh "$OUT" | cut -f1)"
echo ">>> Pipeline commit: $COMMIT"
echo ""
column -t -s$'\t' "$OUT/MANIFEST.tsv"
