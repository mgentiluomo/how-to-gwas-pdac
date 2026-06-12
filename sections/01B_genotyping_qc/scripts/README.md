# Genotyping QC Scripts

This folder contains annotated scripts for the **full sample and variant-level QC pipeline**.

## Scripts (in order)

| # | Script | Purpose | Input | Output | Key Threshold |
|---|--------|---------|-------|--------|---------------|
| 1 | `01_initial_qc_stats.sh` | Baseline statistics | pdac_demo.bed/bim/fam | afreq, het, imiss, lmiss | — |
| 2 | `02_sample_callrate.sh` | Remove high-missingness samples | Stats from 01 | imiss filtered | --mind 0.02 |
| 3 | `03_sex_check.sh` | Identify sex discordance | Dataset from 02 | sexcheck, discordant list | F stat anomaly |
| 4 | `04_heterozygosity.sh` | Remove contamination/outliers | Dataset from 02 | het outliers list, plots | F ± 3 SD |
| 5 | `05_variant_callrate.sh` | Remove low-genotype variants | Dataset from 02 | filtered variants | --geno 0.05 |
| 6 | `06_hardy_weinberg.sh` | Remove HWE-deviant variants | Dataset from 05 | HWE filtered (controls) | p > 1e-6 |
| 7 | `07_relatedness.sh` | Prune related samples (IBD) | Dataset from 05 | pruned samples list | PI_HAT > 0.1875 |
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
bash wsl_setup.sh    # Install dependencies (one-time)
```

See [`../../WSL_SETUP.md`](../../WSL_SETUP.md) for detailed Windows/WSL instructions.

### Linux/Mac Users

Ensure you have PLINK2 and R installed:

```bash
# Ubuntu/Debian
sudo apt-get install plink2 r-base

# macOS (with Homebrew)
brew install plink2
```

## Running the Full Pipeline

Each script is **self-contained**: it outputs the command for the next step at the end.

```bash
bash 01_initial_qc_stats.sh     # ~1 min
bash 02_sample_callrate.sh      # ~1 min
bash 03_sex_check.sh            # ~1 min
bash 04_heterozygosity.sh       # ~2 min (includes R visualization)
bash 05_variant_callrate.sh     # ~1 min
bash 06_hardy_weinberg.sh       # ~2 min (separate by case/control)
bash 07_relatedness.sh          # ~5 min (kinship computation)
bash 08_maf_filter.sh           # ~1 min
bash 09_qc_summary.sh           # ~1 min (report final numbers)
```

## Reproducibility & Customization

All scripts expose a `--seed` argument:
```bash
bash 01_initial_qc_stats.sh 2026
```

And accept dataset paths as arguments (for chaining or reuse):
```bash
bash 02_sample_callrate.sh /path/to/dataset
```

## Output Files

After running the full pipeline, you will have:

```
pdac_demo_qc.*.afreq       — allele frequencies
pdac_demo_qc.*.imiss       — per-sample missing rates
pdac_demo_qc.*.lmiss       — per-variant missing rates
pdac_demo_qc.*.het         — heterozygosity (F statistic)
pdac_demo_qc.*.sexcheck    — sex concordance check
pdac_demo_qc.*.het.outliers.txt — list of het outliers to remove
pdac_demo_qc.*.hwe         — Hardy-Weinberg p-values
pdac_demo_qc.*.king.gz     — kinship coefficients (binary format)
pdac_demo_qc.*.king_removal.txt — list of related samples to remove
pdac_demo_qc.*.afreq.filtered   — final allele frequencies (post-QC)
pdac_demo_qc_summary.txt   — final report
pdac_demo_qc_decision_tree.png — decision tree visualization
```

## Context: PDAC and QC Choices

For **pancreatic cancer GWAS**:
- Small sample size (~250 EUR cases) makes each sample/variant valuable
- Rare variants contribute to heritability in complex disease; don't over-filter MAF
- Use **low MAF threshold** (0.01) to retain signal from rare variants
- Separate HWE p-values for cases and controls (cases may deviate due to disease association)
- Document thresholds and filtering rationale for reproducibility and transparency
