#!/usr/bin/env bash

################################################################################
# Section 1A: Genotyping technologies — Step 03: to binary, merged and repaired
#
# PURPOSE:
#   Take the per-sample long-format files from Step 02 and produce one merged
#   PLINK binary dataset, then repair the allele coding that merging leaves
#   incomplete.
#
#       per-sample .lgen/.map/.fam
#           --(plink --lfile --recode)-->      per-sample .ped/.map
#           --(plink --merge-list --recode)--> merged .ped/.map
#           --(fix alleles, --make-bed)-->     merged .bed/.bim/.fam
#
# WHY PER SAMPLE THEN MERGE:
#   Final Reports usually arrive one file per sample, and genotyping runs arrive
#   months or years apart. Converting each sample independently and merging at
#   the end means a new run is added by converting it and re-running the merge,
#   rather than reprocessing everything. It also stops one bad sample from
#   taking down the whole batch.
#
# THE REPAIR THAT IS EASY TO MISS:
#   A variant that is monomorphic across the samples in hand shows only one
#   allele, so PLINK writes "0" in the other allele column of the .bim. That
#   zero is not missing data, it is an unobserved allele, and it breaks later
#   merges, strand checks and imputation. The manifest knows what the second
#   allele should be, so this script fills it in: for each variant with a zero,
#   it takes the observed allele, looks up the design pair, and writes the other
#   one. This is what makes the dataset combinable with anything else.
#
# INPUT:
#   - results/raw/lgen/*.lgen|.map|.fam                 (Step 02)
#   - results/raw/pdac_demo_01A_reference_alleles.csv   (Step 02)
#
# OUTPUT:
#   - results/raw/binary/*.ped|.map                     per sample
#   - results/raw/merged/DB_merged_final.bed|.bim|.fam
#   - results/raw/pdac_demo_01A_format_sizes.tsv
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
LGEN_DIR="${OUT_DIR}/lgen"
BIN_DIR="${OUT_DIR}/binary"
MRG_DIR="${OUT_DIR}/merged"
REF="${OUT_DIR}/pdac_demo_01A_reference_alleles.csv"

mkdir -p "$BIN_DIR" "$MRG_DIR"

SAMPLES=$(ls "${LGEN_DIR}"/*.lgen 2>/dev/null | xargs -n1 basename | sed 's/\.lgen$//' || true)
if [ -z "$SAMPLES" ]; then
  echo "No .lgen files in ${LGEN_DIR}. Run Step 02 first."
  exit 1
fi
N=$(echo "$SAMPLES" | wc -l)

echo ""
echo "=== Converting ${N} sample(s) to text format ==="
echo ""
for S in $SAMPLES; do
  plink --lfile "${LGEN_DIR}/${S}" --recode --allow-no-sex \
        --out "${BIN_DIR}/${S}" > /dev/null 2>&1
  echo "  ${S}"
done
rm -f "${BIN_DIR}"/*.nosex

echo ""
echo "=== Merging ==="
echo ""
FIRST=$(echo "$SAMPLES" | head -1)

for S in $SAMPLES; do
  plink --file "${BIN_DIR}/${S}" --make-bed --allow-no-sex \
        --out "${BIN_DIR}/${S}" > /dev/null 2>&1
done

if [ "$N" -gt 1 ]; then
  # Samples are added one at a time. Merging incrementally is slower than a
  # single merge-list, but when a conflict appears you know exactly which sample
  # caused it, and you can correct that sample alone. Flipping every dataset at
  # once resolves one conflict and creates others.
  cp "${BIN_DIR}/${FIRST}.bed" "${MRG_DIR}/acc.bed"
  cp "${BIN_DIR}/${FIRST}.bim" "${MRG_DIR}/acc.bim"
  cp "${BIN_DIR}/${FIRST}.fam" "${MRG_DIR}/acc.fam"
  echo "  base: ${FIRST}"

  for S in $(echo "$SAMPLES" | tail -n +2); do
    set +e
    plink --bfile "${MRG_DIR}/acc" --bmerge "${BIN_DIR}/${S}" \
          --make-bed --allow-no-sex --out "${MRG_DIR}/acc_new" > /dev/null 2>&1
    OK=$?
    set -e

    if [ "$OK" -ne 0 ] && [ -f "${MRG_DIR}/acc_new-merge.missnp" ]; then
      NC=$(wc -l < "${MRG_DIR}/acc_new-merge.missnp")
      echo "  ${S}: merge conflict on ${NC} variants"
      echo "    A variant with three or more alleles once two samples are put"
      echo "    together does not exist biologically. It means the same variant"
      echo "    was coded differently in the two files, and the usual cause is a"
      echo "    strand difference. Excluding these variants hides the problem;"
      echo "    flipping the strand of the offending sample tests it."
      plink --bfile "${BIN_DIR}/${S}" --flip "${MRG_DIR}/acc_new-merge.missnp" \
            --make-bed --allow-no-sex --out "${BIN_DIR}/${S}_flip" > /dev/null 2>&1
      set +e
      plink --bfile "${MRG_DIR}/acc" --bmerge "${BIN_DIR}/${S}_flip" \
            --make-bed --allow-no-sex --out "${MRG_DIR}/acc_new" > /dev/null 2>&1
      OK=$?
      set -e
      if [ "$OK" -eq 0 ]; then
        echo "    resolved by flipping ${NC} variants in ${S}: it was a strand difference"
      else
        echo "    NOT resolved by flipping. The samples were converted against"
        echo "    different references. Fix that upstream rather than here."
        exit 1
      fi
    else
      echo "  ${S}: merged cleanly"
    fi

    mv "${MRG_DIR}/acc_new.bed" "${MRG_DIR}/acc.bed"
    mv "${MRG_DIR}/acc_new.bim" "${MRG_DIR}/acc.bim"
    mv "${MRG_DIR}/acc_new.fam" "${MRG_DIR}/acc.fam"
  done

  mv "${MRG_DIR}/acc.bed" "${MRG_DIR}/DB_merged.bed"
  mv "${MRG_DIR}/acc.bim" "${MRG_DIR}/DB_merged.bim"
  mv "${MRG_DIR}/acc.fam" "${MRG_DIR}/DB_merged.fam"
else
  cp "${BIN_DIR}/${FIRST}.bed" "${MRG_DIR}/DB_merged.bed"
  cp "${BIN_DIR}/${FIRST}.bim" "${MRG_DIR}/DB_merged.bim"
  cp "${BIN_DIR}/${FIRST}.fam" "${MRG_DIR}/DB_merged.fam"
  echo "  single sample: nothing to merge"
fi

echo ""
echo "=== Repairing unobserved alleles from the manifest ==="
echo ""
Rscript --vanilla -e "
bim <- read.delim('${MRG_DIR}/DB_merged.bim', header = FALSE, stringsAsFactors = FALSE)
ref <- read.csv('${REF}', stringsAsFactors = FALSE)
before <- sum(bim\$V5 == '0')
m <- ref[match(bim\$V2, ref\$Name), c('A1', 'A2')]
fix <- bim\$V5 == '0' & !is.na(m\$A1)
bim\$V5[fix] <- ifelse(bim\$V6[fix] == m\$A1[fix], m\$A2[fix],
                ifelse(bim\$V6[fix] == m\$A2[fix], m\$A1[fix], bim\$V5[fix]))
after <- sum(bim\$V5 == '0')
write.table(bim, '${MRG_DIR}/DB_merged.bim', col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = '\t')
cat('  variants with an unobserved allele:', before, '\n')
cat('  repaired from the manifest:        ', before - after, '\n')
cat('  still unresolved:                  ', after, '\n')
"

plink --bfile "${MRG_DIR}/DB_merged" --make-bed --freq --allow-no-sex \
      --out "${MRG_DIR}/DB_merged_final" > /dev/null 2>&1

echo ""
echo "=== Size of the same data in each format ==="
echo ""
{
  echo -e "format\tfiles\tbytes"
  RAW=$(cat example_data/*_final_report.csv 2>/dev/null | wc -c)
  LG=$(cat "${LGEN_DIR}"/*.lgen "${LGEN_DIR}"/*.map "${LGEN_DIR}"/*.fam | wc -c)
  PD=$(cat "${BIN_DIR}"/*.ped "${BIN_DIR}"/*.map 2>/dev/null | wc -c)
  BN=$(( $(stat -c%s "${MRG_DIR}/DB_merged_final.bed") + $(stat -c%s "${MRG_DIR}/DB_merged_final.bim") + $(stat -c%s "${MRG_DIR}/DB_merged_final.fam") ))
  echo -e "final report\tcsv\t${RAW}"
  echo -e "plink long\tlgen+map+fam\t${LG}"
  echo -e "plink text\tped+map\t${PD}"
  echo -e "plink binary\tbed+bim+fam\t${BN}"
  echo -e "genotypes only\tbed\t$(stat -c%s "${MRG_DIR}/DB_merged_final.bed")"
} | tee "${OUT_DIR}/pdac_demo_01A_format_sizes.tsv"

echo ""
echo "Compare the last row with the first. The genotypes occupy a small fraction"
echo "of the text export; the rest was identifiers, separators and scores"
echo "repeated on every line."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  bash scripts/01A_study_design/04_preflight_and_formats.sh"
echo ""
