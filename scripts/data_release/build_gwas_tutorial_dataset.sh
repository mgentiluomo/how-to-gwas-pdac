#!/usr/bin/env bash
################################################################################
# build_gwas_tutorial_dataset.sh
#
# First step of the genotype layer of the demo dataset: which individuals of
# the HGDP and 1000 Genomes panels are in it, which ancestry label each one
# carries, and which variants. It writes hgdp1kg_gsa_tutorial.bed/bim/fam and
# sample_ancestry.tsv with the SOURCE sample identifiers.
#
# THIS IS A RECORD, NOT A STEP OF THE PIPELINE. It is the script that built
# the dataset on 4 June 2026 (date of the script and of its log on the build
# server). It was copied here on
# 24 November 2026 and has NOT been run again; scripts/dev/run_all.sh does
# not call it. The commands below are those that were run. Only the five
# path assignments at the top were changed: they were fixed paths on the
# build server and now come from the environment (MD5 of the script as run:
# 6c9a8f2ac0becfdb53e388552e6ff893). The intermediate files in /tmp are as
# in the original.
#
# WHAT IS STILL NOT IN THE REPOSITORY
#   The steps after this one: the anonymisation of the identifiers (MG####;
#   the mapping is withheld on purpose), the injected missingness, the
#   simulated chromosome X. The phenotype is simulated by
#   demo_dataset/scripts/04_simulate_phenotype.R.
#
# INPUT, all public
#   1. Phased genotypes: "HGDP+1KG Phased Reference Panel for Local Ancestry
#      Inference (hg38)", Zenodo, version v1.0, published 2026-01-05,
#      https://doi.org/10.5281/zenodo.18156285
#      (files hgdp1kgp_chr<N>.shapeit5_phased.filter1_SNP_maf005.rechr.vcf.gz,
#      saved here as hgdp1kgp_chr<N>.vcf.gz in VCF_DIR; checked on
#      24 November 2026: the sizes of chromosomes 4, 5 and 9 on the build
#      server equal those listed by Zenodo). The panel holds 4,091 samples.
#      It derives from the gnomAD v3.1 HGDP+1KG callset (Koenig et al. 2024).
#   2. Sample metadata (META): gnomAD v3.1 HGDP+1KG, gnomad_meta_v1.tsv,
#      https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/metadata_and_qc/gnomad_meta_v1.tsv
#      (MD5 c048408444fa4a4d3c91422d2f314f6e, 4,150 samples, 184 columns).
#      The copy used on the build server (MD5
#      35844f5782df8f1c5c8ea5d18dbb966e) is that file with 14 further columns
#      appended, of an origin that is not recorded; this script does not read
#      them. Checked on 24 November 2026: the 184 columns of the public file
#      are identical in the copy for all 4,150 samples.
#   3. GSA_BED: positions of the array, GRCh38, two columns read (chromosome,
#      1-based position). MD5 of the file used:
#      8496bb41513a96db1556f602ef20a5d8. How it was derived from the
#      manufacturer's manifest is not recorded.
#
# THE ANCESTRY LABELS
#   The label of an individual is the genetic-ancestry group gnomAD inferred
#   for the sample, column population_inference.pop of the metadata: nfe and
#   fin are written as eur, afr as afr, eas as eas; every other group is left
#   out. Only samples with high_quality equal to true are considered. The
#   label is NOT the population or the region the sample was collected under
#   (hgdp_tgp_meta.Population, hgdp_tgp_meta.Genetic.region): one individual
#   of population PUR and one of population Mozabite carry the label afr
#   because gnomAD infers them so. sections/08_replicates/scripts/
#   09_set_composition.R checks the rule on the released labels (it holds
#   for 1,461 of 1,461 individuals).
#
# THE DRAW OF 350 AFRICAN AND 350 EAST ASIAN INDIVIDUALS
#   All eligible European samples are taken (767). Of the 896 eligible
#   African and 801 eligible East Asian samples, 350 each are drawn by a
#   shuffle written in awk, srand(2026) and srand(2027). The sequence of
#   rand() after srand(seed) is NOT guaranteed to be the same in every awk
#   implementation (they do not share one generator), and the input order is
#   that of the metadata file, so the same draw is to be expected only with
#   the same awk and the same file. On the build server (bio-plink3) awk is
#   GNU Awk 5.0.1 (/usr/bin/gawk) and bcftools is 1.15.1, as found on
#   24 November 2026; the versions at the time of the build, and the PLINK 2
#   version used for the last command, were not recorded.
#
# WHY 349 AFRICAN INDIVIDUALS AND NOT 350, AND 762 EUROPEAN AND NOT 767
#   The 1,467 samples selected (767 + 350 + 350) are intersected with the
#   samples of the chromosome 22 file of the phased panel. Six are not in
#   the panel: five European and one of the 350 African samples drawn. 1,461
#   remain: 762 eur, 349 afr, 350 eas. Verified on 24 November 2026 from the
#   intermediate lists the build left on the server (samples_keep_raw.txt,
#   1,467 lines; samples_keep.txt, 1,461), without rerunning anything; the
#   counts printed by the build ("EUR: 767, AFR: 896, EAS: 801", "Final
#   samples: 1461") are in its log.
#
# VARIANTS
#   Per chromosome: positions of GSA_BED, biallelic SNPs, minor allele count
#   at least 1 among the selected samples; then PLINK 2 with --maf 0.01
#   --max-alleles 2 --snps-only just-acgt.
#
# USAGE (not needed to follow the tutorial: the release holds the result)
#   GWAS_BUILD_BASE=/path/to/workdir \
#   GWAS_BUILD_VCF_DIR=/path/to/zenodo_vcfs \
#   GWAS_BUILD_META=/path/to/gnomad_meta_v1.tsv \
#   GWAS_BUILD_GSA_BED=/path/to/gsa_positions_grch38.bed \
#     bash scripts/data_release/build_gwas_tutorial_dataset.sh
#   Needs bcftools, plink2 and awk on PATH. Output in
#   ${GWAS_BUILD_BASE}/processed/final_dataset (GWAS_BUILD_FINAL to change it).
################################################################################
set -euo pipefail
BASE="${GWAS_BUILD_BASE:?set GWAS_BUILD_BASE to a working directory}"
VCF_DIR="${GWAS_BUILD_VCF_DIR:-${BASE}/raw_data/zenodo_hgdp1kg}"
META="${GWAS_BUILD_META:-${BASE}/sample_meta.tsv}"
GSA_BED="${GWAS_BUILD_GSA_BED:-${BASE}/gsa_positions_grch38.bed}"
FINAL="${GWAS_BUILD_FINAL:-${BASE}/processed/final_dataset}"
TMP="${GWAS_BUILD_TMP:-${BASE}/tmp}"
mkdir -p "${BASE}/processed" "${FINAL}" "${TMP}"
echo ">>> [1/5] Building sample lists..."
rm -f /tmp/eur_hq.txt /tmp/afr_hq.txt /tmp/eas_hq.txt
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; next} $(col["high_quality"])!="true"{next} {pop=$(col["population_inference.pop"]); sid=$1; if(pop=="nfe"||pop=="fin") print sid>>"/tmp/eur_hq.txt"; else if(pop=="afr") print sid>>"/tmp/afr_hq.txt"; else if(pop=="eas") print sid>>"/tmp/eas_hq.txt"}' "${META}"
echo "  EUR: $(wc -l < /tmp/eur_hq.txt), AFR: $(wc -l < /tmp/afr_hq.txt), EAS: $(wc -l < /tmp/eas_hq.txt)"
awk -v seed=2026 -v n=350 'BEGIN{srand(seed)} {a[NR]=$0} END{for(i=NR;i>1;i--){j=int(rand()*(i-1))+1;t=a[i];a[i]=a[j];a[j]=t}; for(i=1;i<=n&&i<=NR;i++) print a[i]}' /tmp/afr_hq.txt > /tmp/afr_keep.txt
awk -v seed=2027 -v n=350 'BEGIN{srand(seed)} {a[NR]=$0} END{for(i=NR;i>1;i--){j=int(rand()*(i-1))+1;t=a[i];a[i]=a[j];a[j]=t}; for(i=1;i<=n&&i<=NR;i++) print a[i]}' /tmp/eas_hq.txt > /tmp/eas_keep.txt
cat /tmp/eur_hq.txt /tmp/afr_keep.txt /tmp/eas_keep.txt | awk '{print $1, $1}' > "${TMP}/samples_keep_raw.txt"
bcftools query -l "${VCF_DIR}/hgdp1kgp_chr22.vcf.gz" > /tmp/vcf_sample_list.txt
awk 'NR==FNR{v[$1]=1;next} v[$1]{print}' /tmp/vcf_sample_list.txt "${TMP}/samples_keep_raw.txt" > "${TMP}/samples_keep.txt"
awk '{print $1}' "${TMP}/samples_keep.txt" > "${TMP}/samples_onecol.txt"
echo "  Final samples: $(wc -l < ${TMP}/samples_keep.txt)"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; next} {pop=$(col["population_inference.pop"]); if(pop=="nfe"||pop=="fin") pop="eur"; print $1"\t"pop}' "${META}" | awk 'NR==FNR{keep[$1]=1;next} keep[$1]{print}' "${TMP}/samples_onecol.txt" - > "${FINAL}/sample_ancestry.tsv"
echo ">>> [2/5] Indexing VCFs..."
for chr in $(seq 1 22); do
  vcf="${VCF_DIR}/hgdp1kgp_chr${chr}.vcf.gz"
  if [ ! -f "${vcf}.csi" ]; then bcftools index "${vcf}"; fi
done
echo "  Done."
echo ">>> [3/5] Filtering per chromosome (-T, 1-based)..."
> "${TMP}/vcf_merge_list.txt"
for chr in $(seq 1 22); do
  echo "  chr${chr}..."
  awk -v c="${chr}" -F'\t' '$1==c{print $1"\t"$2}' "${GSA_BED}" > "${TMP}/gsa_chr${chr}.tsv"
  bcftools view -T "${TMP}/gsa_chr${chr}.tsv" -S "${TMP}/samples_onecol.txt" --min-ac 1 --max-alleles 2 --type snps -O z -o "${TMP}/chr${chr}_filtered.vcf.gz" "${VCF_DIR}/hgdp1kgp_chr${chr}.vcf.gz"
  bcftools index "${TMP}/chr${chr}_filtered.vcf.gz"
  echo "${TMP}/chr${chr}_filtered.vcf.gz" >> "${TMP}/vcf_merge_list.txt"
done
echo ">>> [4/5] Merging and converting..."
bcftools concat --file-list "${TMP}/vcf_merge_list.txt" -O z -o "${TMP}/all_merged.vcf.gz"
bcftools index "${TMP}/all_merged.vcf.gz"
bcftools stats "${TMP}/all_merged.vcf.gz" | grep "^SN" | grep -E "samples|records"
plink2 --vcf "${TMP}/all_merged.vcf.gz" --keep "${TMP}/samples_keep.txt" --double-id --maf 0.01 --max-alleles 2 --snps-only just-acgt --make-bed --out "${FINAL}/hgdp1kg_gsa_tutorial" --threads 4
echo ">>> [5/5] Done."
echo "  Variants: $(wc -l < ${FINAL}/hgdp1kg_gsa_tutorial.bim)"
echo "  Samples:  $(wc -l < ${FINAL}/hgdp1kg_gsa_tutorial.fam)"
ls -lh "${FINAL}/hgdp1kg_gsa_tutorial."*
