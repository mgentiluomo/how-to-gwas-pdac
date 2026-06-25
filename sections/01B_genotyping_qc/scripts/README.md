# Genotyping QC Scripts

This folder contains annotated scripts for the **full sample and variant-level QC pipeline**.

## Scripts (in order)

| # | Script | Purpose | Input | Output | Key Threshold |
|---|--------|---------|-------|--------|---------------|
| 1 | `01_initial_qc_stats.sh` | Baseline statistics and plots | `demo_data/pdac_demo.bed/bim/fam` | afreq, het, smiss, vmiss, PNG plots | — |
| 2 | `02_sample_callrate.sh` | Remove high-missingness samples | Raw demo data | sample-filtered dataset | --mind 0.02 |
| 3 | `03_sex_check.sh` | Identify sex discordance | Dataset from 02 | sexcheck, discordant list | F stat anomaly |
| 4 | `04_heterozygosity.sh` | Remove contamination/outliers | Dataset from 03 | het outliers list, plots | F ± 3 SD |
| 5 | `05_variant_callrate.sh` | Remove low-genotype variants | Dataset from 04 | filtered variants | --geno 0.05 |
| 6 | `06_hardy_weinberg.sh` | Remove HWE-deviant variants | Dataset from 05 | HWE filtered (controls) | p > 1e-6 |
| 7 | `07_relatedness.sh` | Prune related samples (KING) | Dataset from 06 | phenotype-aware pruned samples list | KING kinship > 0.1875 |
| 8 | `08_maf_filter.sh` | Apply MAF threshold (context-dependent) | Dataset from 07 | final variants | --maf 0.01 (PDAC) |
| 9 | `09_qc_summary.sh` + `09_qc_report.R` | Summarize QC impact | All datasets | combined PDF report, count table, text summary | — |

## System Requirements & Setup

### ⚠️ Important for Windows Users

**This pipeline requires PLINK2, R, and bash — all work best in a Linux environment.**

If you're on Windows:
- **❌ NOT recommended:** Running directly in PowerShell/CMD (fragile, requires manual binary downloads)
- **✅ STRONGLY RECOMMENDED:** Use WSL2 (Windows Subsystem for Linux)

#### Quick WSL Setup

```bash
# In PowerShell (Administrator):
wsl --install

# Then, in WSL terminal:
cd /mnt/s/Github/how-to-gwas-pdac
bash scripts/dev/tools_setup.sh    # Install PLINK tools (one-time)
```

See [`../../WSL_SETUP.md`](../../WSL_SETUP.md) for detailed Windows/WSL instructions.

### Linux/Mac Users

Run `bash scripts/dev/tools_setup.sh` from the tutorial project root to install/check PLINK2, PLINK1.9, METAL, REGENIE, and R.

## Running the Full Pipeline

Each script is **self-contained**: it outputs the command for the next step at the end.
Steps 06 and 07 use `demo_data/phenotype.txt`, which must contain `FID IID PHENO`
with `PHENO` coded as `1 = control` and `2 = case`.

```bash
bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh     # ~1 min (includes R visualization)
bash scripts/01B_genotyping_qc/02_sample_callrate.sh      # ~1 min
bash scripts/01B_genotyping_qc/03_sex_check.sh            # ~1 min
bash scripts/01B_genotyping_qc/04_heterozygosity.sh       # ~2 min (includes R visualization)
bash scripts/01B_genotyping_qc/05_variant_callrate.sh     # ~1 min
bash scripts/01B_genotyping_qc/06_hardy_weinberg.sh       # ~2 min
bash scripts/01B_genotyping_qc/07_relatedness.sh          # ~5 min (kinship computation)
bash scripts/01B_genotyping_qc/08_maf_filter.sh           # ~1 min
bash scripts/01B_genotyping_qc/09_qc_summary.sh           # ~1 min (combined QC report)
```

## Reproducibility & Customization

All scripts expose a `--seed` argument:
```bash
bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh 2026
```

And accept dataset paths as arguments (for chaining or reuse):
```bash
bash scripts/01B_genotyping_qc/02_sample_callrate.sh 2026 /path/to/raw_data
```

## Output Files

After running the full pipeline, you will have:

```
results/qc/pdac_demo_01_qc.afreq       — allele frequencies
results/qc/pdac_demo_01_qc.smiss       — per-sample missing rates
results/qc/pdac_demo_01_qc.vmiss       — per-variant missing rates
results/qc/pdac_demo_01_qc_variant_density_by_chromosome.tsv
results/qc/pdac_demo_01_qc_variant_density_by_chromosome.png
results/qc/pdac_demo_01_qc_sample_missingness_histogram.png
results/qc/pdac_demo_01_qc_variant_missingness_histogram.png
results/qc/pdac_demo_01_qc_allele_frequency_distribution.png
results/qc/pdac_demo_01_qc_heterozygote_rate_distribution.png
results/qc/pdac_demo_01_qc_sample_missingness_vs_heterozygosity.png
results/qc/pdac_demo_04_het.het        — heterozygosity (F statistic)
results/qc/pdac_demo_03_sexcheck.sexcheck — sex concordance check
results/qc/pdac_demo_04_het_outliers.txt — list of het outliers to remove
results/qc/pdac_demo_06_hwe.hardy      — Hardy-Weinberg p-values
results/qc/pdac_demo_06_hwe_pvalue_distribution.png — HWE p-value histogram
results/qc/pdac_demo_07_king_cutoff.king.cutoff.in.id — PLINK2 --king-cutoff 0.1875 keep list
results/qc/pdac_demo_07_king_cutoff.king.cutoff.out.id — PLINK2 --king-cutoff 0.1875 removal list
results/qc/pdac_demo_07_ibd_pi_hat.genome — PLINK1.9 --genome PI_HAT report
results/qc/pdac_demo_07_ibd_pi_hat_annotated.tsv — PI_HAT pairs with category labels
results/qc/pdac_demo_07_ibd_pi_hat_summary.tsv — PI_HAT relationship category counts
results/qc/pdac_demo_07_related_pairs.kin0 — KING table for third-degree-or-closer pairs
results/qc/pdac_demo_07_related_pairs_annotated.tsv — related pairs with category labels
results/qc/pdac_demo_07_relatedness_summary.tsv — duplicate/twin, degree, and unrelated sample counts
results/qc/pdac_demo_07_king_vs_pihat_comparison.tsv — pair-level KING/PI_HAT comparison
results/qc/pdac_demo_07_king_vs_pihat_summary.tsv — KING/PI_HAT agreement summary
results/qc/pdac_demo_07_pruning_components.tsv — connected components among pairs above the pruning threshold
results/qc/pdac_demo_07_relatedness_pruning_diagnostics.tsv — headline diagnostics for sample loss
results/qc/pdac_demo_07_removal_decisions.tsv — phenotype-aware relatedness pruning decisions
results/qc/pdac_demo_07_removal_summary.tsv — case/control removal counts
results/qc/pdac_demo_07_removal_phenotype_counts.tsv — final case/control counts from removal list
results/qc/pdac_demo_07_removal_list.txt — samples removed by phenotype-aware pruning
results/qc/pdac_demo_08_afreq.afreq    — final allele frequencies (post-QC)
results/qc/pdac_demo_09_qc_counts.tsv  — final sample/variant count table
results/qc/pdac_demo_09_qc_summary.txt — plain-text QC summary
results/qc/pdac_demo_09_qc_report.pdf  — combined professional QC report
```

## Context: PDAC and QC Choices

For **pancreatic cancer GWAS**:
- Small sample size (~250 EUR cases) makes each sample/variant valuable
- Rare variants contribute to heritability in complex disease; don't over-filter MAF
- Use **low MAF threshold** (0.01) to retain signal from rare variants
- Separate HWE p-values for cases and controls (cases may deviate due to disease association)
- Prefer removing controls over cases during relatedness pruning when a related pair is case-control
- Use KING kinship >0.1875 for final removal, while also reporting PLINK1.9 PI_HAT for teaching and review
- Remember that PI_HAT and KING kinship are related but not interchangeable scales
- Document thresholds and filtering rationale for reproducibility and transparency
