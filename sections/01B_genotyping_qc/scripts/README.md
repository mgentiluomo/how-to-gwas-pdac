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
| 7 | `07_relatedness.sh` | Prune related samples (IBD) | Dataset from 06 | pruned samples list | PI_HAT > 0.1875 |
| 8 | `08_maf_filter.sh` | Apply MAF threshold (context-dependent) | Dataset from 07 | final variants | --maf 0.01 (PDAC) |
| 9 | `09_qc_summary.sh` | Summarize QC impact | All datasets | summary table, decision tree | — |

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

```bash
bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh     # ~1 min (includes R visualization)
bash scripts/01B_genotyping_qc/02_sample_callrate.sh      # ~1 min
bash scripts/01B_genotyping_qc/03_sex_check.sh            # ~1 min
bash scripts/01B_genotyping_qc/04_heterozygosity.sh       # ~2 min (includes R visualization)
bash scripts/01B_genotyping_qc/05_variant_callrate.sh     # ~1 min
bash scripts/01B_genotyping_qc/06_hardy_weinberg.sh       # ~2 min
bash scripts/01B_genotyping_qc/07_relatedness.sh          # ~5 min (kinship computation)
bash scripts/01B_genotyping_qc/08_maf_filter.sh           # ~1 min
bash scripts/01B_genotyping_qc/09_qc_summary.sh           # ~1 min (report final numbers)
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
results/qc/pdac_demo_04_het.het        — heterozygosity (F statistic)
results/qc/pdac_demo_03_sexcheck.sexcheck — sex concordance check
results/qc/pdac_demo_04_het_outliers.txt — list of het outliers to remove
results/qc/pdac_demo_06_hwe.hardy      — Hardy-Weinberg p-values
results/qc/pdac_demo_07_kinship.*      — kinship output files
results/qc/pdac_demo_08_afreq.afreq    — final allele frequencies (post-QC)
results/qc/pdac_demo_qc_summary.txt    — final report
results/qc/pdac_demo_qc_decision_tree.pdf — decision tree visualization
```

## Context: PDAC and QC Choices

For **pancreatic cancer GWAS**:
- Small sample size (~250 EUR cases) makes each sample/variant valuable
- Rare variants contribute to heritability in complex disease; don't over-filter MAF
- Use **low MAF threshold** (0.01) to retain signal from rare variants
- Separate HWE p-values for cases and controls (cases may deviate due to disease association)
- Document thresholds and filtering rationale for reproducibility and transparency
